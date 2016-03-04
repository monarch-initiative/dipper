import re
import logging
import hashlib

from rdflib import Namespace, URIRef, Literal
from rdflib.namespace import RDF, OWL, RDFS, XSD

from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

__author__ = 'nlw'

logger = logging.getLogger(__name__)


class Assoc:
    """
    An abstract class for OBAN (Monarch)-style associations,
    to enable attribution of source and evidence
    on statements.
    """

    assoc_types = {
        'association': 'OBAN:association'
    }

    annotation_properties = {
        'replaced_by': 'IAO:0100001',
        'consider': 'OIO:consider',
        'hasExactSynonym': 'OIO:hasExactSynonym',
        'hasRelatedSynonym': 'OIO:hasRelatedSynonym',
        'definition': 'IAO:0000115',
        'has_xref': 'OIO:hasDbXref',
    }

    object_properties = {
        'has_disposition': 'GENO:0000208',
        'has_phenotype': 'RO:0002200',
        'in_taxon': 'RO:0002162',
        'has_quality': 'RO:0000086',
        'towards': 'RO:0002503',
        'has_subject': 'OBAN:association_has_subject',
        'has_object': 'OBAN:association_has_object',
        'has_predicate': 'OBAN:association_has_object_property',
        'is_about': 'IAO:00000136',
        'has_evidence': 'RO:0002558',
        'has_source': 'dc:source',
        'has_provenance': 'OBAN:has_provenance'
    }

    datatype_properties = {
        'position': 'faldo:position',
        'has_measurement': 'IAO:0000004'
    }

    properties = annotation_properties.copy()
    properties.update(object_properties)
    properties.update(datatype_properties)

    OWLCLASS = OWL['Class']
    OWLIND = OWL['NamedIndividual']
    OBJECTPROP = OWL['ObjectProperty']
    ANNOTPROP = OWL['AnnotationProperty']
    DATAPROP = OWL['DatatypeProperty']

    SUBCLASS = RDFS['subClassOf']
    BASE = Namespace(curie_map.get()[''])

    def __init__(self, definedby):
        self.cu = CurieUtil(curie_map.get())
        self.gu = GraphUtils(curie_map.get())

        # core parts of the association
        self.definedby = definedby
        self.sub = self.obj = self.rel = None
        self.assoc_id = None

        self.description = None
        self.source = []
        self.evidence = []
        # this is going to be used for the refactored evidence/provenance
        self.provenance = []

        self.score = None
        self.score_type = None
        self.score_unit = None

        return

    def get_properties(self):
        return self.properties

    def _is_valid(self):

        # check if sub/obj/rel are none...throw error
        if self.sub is None:
            raise ValueError('No subject set for this association')
        if self.obj is None:
            raise ValueError('No object set for this association')
        if self.rel is None:
            raise ValueError('No relation set for this association')

        return True

    def _add_basic_association_to_graph(self, g):

        if not self._is_valid():
            return

        # first, add the direct triple
        # anonymous (blank) nodes are indicated with underscore
        s = self.gu.getNode(self.sub)
        o = self.gu.getNode(self.obj)
        p = self.gu.getNode(self.rel)

        g.add((s, p, o))

        if self.assoc_id is None:
            self.set_association_id()

        node = self.gu.getNode(self.assoc_id)
        g.add((node, RDF['type'],
               self.gu.getNode(self.assoc_types['association'])))

        self.gu.addTriple(g, self.assoc_id,
                          self.object_properties['has_subject'], self.sub)
        self.gu.addTriple(g, self.assoc_id,
                          self.object_properties['has_object'], self.obj)
        self.gu.addTriple(g, self.assoc_id,
                          self.object_properties['has_predicate'], self.rel)

        if self.description is not None:
            self.gu.addDescription(g, self.assoc_id, self.description)

        if self.evidence is not None and len(self.evidence) > 0:
            for e in self.evidence:
                self.gu.addTriple(g, self.assoc_id,
                                  self.object_properties['has_evidence'], e)

        if self.source is not None and len(self.source) > 0:
            for s in self.source:
                if re.match('http', s):
                    # TODO assume that the source is a publication?
                    # use Reference class here
                    self.gu.addTriple(g, self.assoc_id,
                                      self.object_properties['has_source'], s,
                                      True)
                else:
                    self.gu.addTriple(g, self.assoc_id,
                                      self.object_properties['has_source'], s)

        if self.provenance is not None and len(self.provenance) > 0:
            for p in self.provenance:
                self.gu.addTriple(g, self.assoc_id,
                                  self.object_properties['has_provenance'], p)

        if self.score is not None:
            self.gu.addTriple(g,
                              self.assoc_id, self.properties['has_measurement'],
                              Literal(self.score, datatype=XSD['float']), True)
            # TODO
            # update with some kind of instance of scoring object
            # that has a unit and type

        return

    def add_association_to_graph(self, g):

        self._add_basic_association_to_graph(g)

        return

    def set_subject(self, identifier):
        self.sub = identifier
        return

    def set_object(self, identifier):
        self.obj = identifier
        return

    def set_relationship(self, identifier):
        self.rel = identifier
        return

    def set_association_id(self, assoc_id=None):
        """
        This will set the association ID based on the internal parts
            of the association.
        To be used in cases where an external association identifier
            should be used.

        :param assoc_id:
        :return:
        """
        if assoc_id is None:
            self.assoc_id = self.make_association_id(self.definedby, self.sub,
                                                     self.rel, self.obj)
        else:
            self.assoc_id = assoc_id

        return

    def get_association_id(self):

        return self.assoc_id

    def set_description(self, description):
        self.description = description

        return

    def set_score(self, score, unit=None, score_type=None):

        self.score = score
        self.score_unit = unit
        self.score_type = score_type

        return

    def add_evidence(self, identifier):
        """
        Add an evidence code to the association object (maintained as a list)
        :param identifier:
        :return:
        """

        if identifier is not None and identifier.strip() != '':
            self.evidence += [identifier]

        return

    def add_source(self, identifier):
        """
        Add a source identifier (such as publication id)
        to the association object (maintained as a list)
        TODO we need to greatly expand this function!

        :param identifier:
        :return:
        """

        if identifier is not None and identifier.strip() != '':
            self.source += [identifier]

        return

    def add_provenance(self, identifier):

        if identifier is not None and identifier.strip() != '':
            self.provenance += [identifier]

        return

    def load_all_properties(self, g):
        props = {
            self.OBJECTPROP: self.object_properties,
            self.ANNOTPROP: self.annotation_properties,
            self.DATAPROP: self.datatype_properties
        }

        for p in props:
            self.gu.loadProperties(g, props[p], p)

        return

    def _get_source_uri(self, pub_id):
        """
        Given some kind of pub_id (which might be a CURIE or url),
        convert it into a proper node.

        :param pub_id:
        :return: source: Well-formed URI for the given identifier (or url)
        """

        source = None
        if re.compile('http').match(pub_id):
            source = URIRef(pub_id)
        else:
            u = self.gu.getNode(pub_id)
            if u is not None:
                source = URIRef(u)
            else:
                logger.error("An id we don't know how to deal with: %s", pub_id)

        return source

    @staticmethod
    def make_association_id(definedby, subject, predicate, object,
                            attributes=None):
        """
        A method to create unique identifiers for OBAN-style associations,
        based on all the parts of the association
        If any of the items is empty or None, it will convert it to blank.
        It effectively md5 hashes the (+)-joined string from the values.
        Subclasses of Assoc can submit an additional array of attributes
        that will be added to the ID.

        :param definedby: The (data) resource that provided the annotation
        :param subject:
        :param predicate:
        :param object:
        :param attributes:
        :return:
        """

        # note others available:
        #   md5(), sha1(), sha224(), sha256(), sha384(), and sha512()
        # TEC: at our scale, md5 is in danger of having collisions.
        # putting definedby first,
        # as this will usually be the datasource providing the annotation
        # this will end up making the first few parts of the id
        # be the same for all annotations in that resource
        items_to_hash = [definedby, subject, predicate, object]
        if attributes is not None:
            items_to_hash += attributes

        for i, val in enumerate(items_to_hash):
            if val is None:
                items_to_hash[i] = ''

        byte_string = '+'.join(items_to_hash).encode("utf-8")

        # TODO put this in a util?
        return ':'.join(('MONARCH', hashlib.md5(byte_string).hexdigest()))
