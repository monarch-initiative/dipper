import logging
import hashlib

from xml.sax import SAXParseException
from collections import defaultdict
from rdflib import URIRef, ConjunctiveGraph, util as rdflib_util
from rdflib.namespace import DC, RDF, OWL

from dipper.utils.CurieUtil import CurieUtil

__author__ = 'nlw'

LOG = logging.getLogger(__name__)


class GraphUtils:

    def __init__(self, curie_map):
        self.curie_map = curie_map
        self.cu = CurieUtil(curie_map)

        return

    @staticmethod
    def write(graph, fileformat=None, filename=None):
        """
        A basic graph writer (to stdout) for any of the sources.
        this will write raw triples in rdfxml, unless specified.
        to write turtle, specify format='turtle'
        an optional file can be supplied instead of stdout
        :return: None

        """

        filewriter = None
        if fileformat is None:
            fileformat = 'turtle'
        if filename is not None:

            with open(filename, 'wb') as filewriter:
                LOG.info("Writing triples in %s to %s", fileformat, filename)
                # rdflib serialize
                graph.serialize(filewriter, format=fileformat)
        else:
            print(graph.serialize(fileformat).decode())
        return

    @staticmethod
    def get_properties_from_graph(graph):
        """
        Wrapper for RDFLib.graph.predicates() that returns a unique set
        :param graph: RDFLib.graph
        :return: set, set of properties
        """
        # collapse to single list
        property_set = list()
        for row in graph.predicates():
            property_set.append(row)

        return set(property_set)

    @staticmethod
    def add_property_axioms(graph, properties):
        ontology_graph = ConjunctiveGraph()
        GH = 'https://raw.githubusercontent.com'
        OBO = 'http://purl.obolibrary.org/obo'
        ontologies = [
            OBO + '/sepio.owl',
            OBO + '/geno.owl',
            OBO + '/iao.owl',
            OBO + '/ero.owl',
            OBO + '/pco.owl',
            OBO + '/xco.owl',
            OBO + '/ro.owl',
            GH + '/jamesmalone/OBAN/master/ontology/oban_core.ttl',
        ]

        # random timeouts can waste hours. (too many redirects?)
        # there is a timeout param in urllib.request,
        # but it is not exposed by rdflib.parsing
        # so retry once on URLError
        for ontology in ontologies:
            LOG.info("parsing: " + ontology)
            try:
                ontology_graph.parse(
                    ontology, format=rdflib_util.guess_format(ontology))
            except SAXParseException as e:
                LOG.error(e)
                LOG.error('Retrying as turtle: ' + ontology)
                ontology_graph.parse(ontology, format="turtle")
            except OSError as e:  # URLError:
                # simple retry
                LOG.error(e)
                LOG.error('Retrying: ' + ontology)
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
            URIRef('https://monarchinitiative.org/MONARCH_cliqueLeader'), RDF['type'],
            OWL['AnnotationProperty']))

        graph.add((
            URIRef('https://monarchinitiative.org/MONARCH_anonymous'), RDF['type'],
            OWL['AnnotationProperty']))

        return graph

    @staticmethod
    def add_property_to_graph(results, graph, property_type, property_list):

        for row in results:
            if row in property_list:
                graph.add((row, RDF['type'], property_type))
        return graph

    @staticmethod
    def digest_id(wordage):   # same as source/Source.hash_id(wordage)
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

    @staticmethod
    def compare_graph_predicates(graph1, graph2):
        '''
        From rdf graphs, count predicates in each and return a list of
        : param graph1 graph, hopefully RDFlib-like
        : param graph2 graph, ditto
        : return dict with count of predicates in each graph:
        : e.g.:
        :         {
        :         "has_a_property": {
        :                 "graph1": 1234,
        :                 "graph2": 1023},
        :         "has_another_property": {
        :                 "graph1": 94,
        :                 "graph2": 51}
        :         }
        '''
        # dict of dicts that acts sensibly when a key that doesn't
        # exist is accessed
        counts = defaultdict(lambda: defaultdict(int))
        for this_g in [graph1, graph2]:
            for this_p in this_g.predicates():
                counts[this_p][str(this_g.identifier)] = \
                    counts[this_p][str(this_g.identifier)] + 1
        return counts

    @staticmethod
    def count_predicates(graph):
        '''
        From rdf graphs, count predicates in each and return a list of
        : param graph
        : return dict with count of predicates in each graph:
        : e.g.:
        :         {
        :         "has_a_property": 1234,
        :         "has_another_property": 482
        :         }
        '''
        # dict of dicts that acts sensibly when a key that doesn't
        # exist is accessed
        counts = defaultdict(int)
        for this_p in graph.predicates():
            counts[this_p] = counts[this_p] + 1
        return counts
