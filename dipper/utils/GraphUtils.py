import logging
import hashlib
from rdflib import URIRef, ConjunctiveGraph
from rdflib import util as rdflib_util
from rdflib.namespace import DC, RDF, OWL
from xml.sax import SAXParseException

from dipper.utils.CurieUtil import CurieUtil

__author__ = 'nlw'

logger = logging.getLogger(__name__)


class GraphUtils:

    def __init__(self, curie_map):
        self.curie_map = curie_map
        self.cu = CurieUtil(curie_map)

        return

    def write(self, graph, fileformat=None, file=None):
        """
        A basic graph writer (to stdout) for any of the sources.
        this will write raw triples in rdfxml, unless specified.
        to write turtle, specify format='turtle'
        an optional file can be supplied instead of stdout
        :return: None

        """

        filewriter = None
        if fileformat is None:
            fileformat = 'rdfxml'
        if file is not None:
            filewriter = open(file, 'wb')

            logger.info("Writing triples in %s to %s", fileformat, file)
            graph.serialize(filewriter, format=fileformat)
            filewriter.close()
        else:
            print(graph.serialize(format=fileformat).decode())
        return

    @staticmethod
    def get_properties_from_graph(graph):
        """
        Wrapper for RDFLib.graph.predicates() that returns a unique set
        :param graph: RDFLib.graph
        :return: set, set of properties
        """
        # collapse to single list
        property_set = set()
        for row in graph.predicates():
            property_set.add(row)

        return property_set

    @staticmethod
    def add_property_axioms(graph, properties):
        ontology_graph = ConjunctiveGraph()
        GH = 'https://raw.githubusercontent.com'
        MI = '/monarch-initiative'
        ontologies = [
            GH + MI + '/SEPIO-ontology/master/src/ontology/sepio.owl',
            GH + MI + '/GENO-ontology/develop/src/ontology/geno.owl',
            GH + '/oborel/obo-relations/master/ro.owl',
            'http://purl.obolibrary.org/obo/iao.owl',
            'http://purl.obolibrary.org/obo/ero.owl',
            GH + '/jamesmalone/OBAN/master/ontology/oban_core.ttl',
            'http://purl.obolibrary.org/obo/pco.owl',
            'http://purl.obolibrary.org/obo/xco.owl'
        ]

        # random timeouts can waste hours. (too many redirects?)
        # there is a timeout param in urllib.request,
        # but it is not exposed by rdflib.parsing
        # so retry once on URLError
        for ontology in ontologies:
            logger.info("parsing: " + ontology)
            try:
                ontology_graph.parse(
                    ontology, format=rdflib_util.guess_format(ontology))
            except SAXParseException as e:
                logger.error(e)
                logger.error('Retrying as turtle: ' + ontology)
                ontology_graph.parse(ontology, format="turtle")
            except OSError as e:  # URLError:
                # simple retry
                logger.error(e)
                logger.error('Retrying: ' + ontology)
                ontology_graph.parse(
                    ontology, format=rdflib_util.guess_format(ontology))

        # Get object properties
        graph = GraphUtils.add_property_to_graph(
            ontology_graph.subjects(RDF['type'], OWL['ObjectProperty']),
            graph, OWL['ObjectProperty'], properties)

        # Get annotation properties
        graph = GraphUtils.add_property_to_graph(
            ontology_graph.subjects(RDF['type'], OWL['AnnotationProperty']),
            graph, OWL['AnnotationProperty'], properties)

        # Get data properties
        graph = GraphUtils.add_property_to_graph(
            ontology_graph.subjects(RDF['type'], OWL['DatatypeProperty']),
            graph, OWL['DatatypeProperty'], properties)

        for row in graph.predicates(DC['source'], OWL['AnnotationProperty']):
            if row == RDF['type']:
                graph.remove(
                    (DC['source'], RDF['type'], OWL['AnnotationProperty']))
        graph.add((DC['source'], RDF['type'], OWL['ObjectProperty']))

        # Hardcoded properties
        graph.add((
            URIRef('https://monarchinitiative.org/MONARCH_cliqueLeader'),
            RDF['type'], OWL['AnnotationProperty']))

        graph.add((URIRef('https://monarchinitiative.org/MONARCH_anonymous'),
                  RDF['type'], OWL['AnnotationProperty']))

        return graph

    @staticmethod
    def add_property_to_graph(results, graph, property_type, property_list):

        for row in results:
            if row in property_list:
                graph.add((row, RDF['type'], property_type))
        return graph

    def digest_id(wordage):
        '''
        Form a deterministic digest of input
        Leading 'b' is an experiment forcing the first char to be non numeric
        but valid hex
        Not required for RDF but some other contexts do not want the leading
        char to be a digit

        : param str wordage arbitrary string
        : return str
        '''
        return 'b' + hashlib.sha1(wordage.encode('utf-8')).hexdigest()[1:20]
