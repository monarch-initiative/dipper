__author__ = 'nlw'

import re
import logging
import hashlib

from rdflib import Namespace, URIRef, Literal
from rdflib.namespace import RDF, DC, OWL, RDFS, XSD

from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


class Assoc:
    """
    An abstract class for OBAN (Monarch)-style associations, to enable attribution of source and evidence
    on statements.
    """

    assoc_types = {
        'association' : 'Annotation:'  # FIXME 'OBAN:association'
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
        'has_subject': ':hasSubject',  # FIXME 'OBAN::association_has_subject'
        'has_object': ':hasObject',  # FIXME 'OBAN:association_has_object'
        'has_predicate': ':hasPredicate',
        'is_about': 'IAO:00000136',
        'has_evidence': 'RO:0002558'
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
        self.definedby = definedby
        self.sub = self.obj = self.rel = None
        self.assoc_id = None
        return

    def get_namespaces(self):
        if self.namespaces:
            return self.namespaces

        return None

    def get_properties(self):
        return self.properties

    def createAssociationNode(self, g):
        #TODO make a general association object following our pattern

        return

    def set_assoc_id(self, assoc_id):
        """
        To be used in cases where an external association identifier should be used.
        :param assoc_id:
        :return:
        """
        self.annot_id = assoc_id

        return

    def addAssociationToGraph(self, g):
        gu = self.gu
        # first, add the direct triple
        # anonymous nodes are indicated with underscore
        s = gu.getNode(self.sub)
        o = gu.getNode(self.obj)
        p = gu.getNode(self.rel)

        g.add((s, p, o))

        # now, create the reified relationship with our annotation pattern
        node = gu.getNode(self.assoc_id)
        g.add((node, RDF['type'], gu.getNode(self.assoc_types['association'])))
        g.add((node, gu.getNode(self.properties['has_subject']), s))
        g.add((node, gu.getNode(self.properties['has_object']), o))
        g.add((node, gu.getNode(self.properties['has_predicate']), p))

        # this is handling the occasional messy pubs that are sometimes literals
        if self.pub_id is not None:
            source = self._get_source(self.pub_id)

        if self.pub_id is not None and self.pub_id.strip() != '':
            self._add_source_node(g, node, source, self.pub_id)

        if self.evidence is None or self.evidence.strip() == '':
            pass
        else:
            self.addEvidence(g,self.evidence,self.annot_id)

        # Check if publications are in list form
        if self.pub_list is not None:
            self._process_pub_list(self.pub_list, g, node)

        return

    def setSubject(self, identifier):
        self.sub = identifier
        return

    def setObject(self, identifier):
        self.obj = identifier
        return

    def setRelationship(self, identifier):
        self.rel = identifier
        return

    def addEvidence(self, g, evidence_identifier, annot_id=None):
        """
        Add the evidence to the annotation object; if one is supplied, add it to that.
        :param g:
        :param evidence_identifier:
        :return:
        """
        evidence = self.gu.getNode(evidence_identifier)
        if (annot_id is None and self.annot_id is not None):
            annot_id = self.annot_id
        if annot_id is not None:
            node = self.gu.getNode(annot_id)
            g.add((node, self.object_properties['has_evidence'], evidence))

        return

    def addSource(self, g, assoc, source_id):
        """
        TODO we need to greatly expand this function!
        :param g:
        :param assoc:
        :param source_identifier:
        :return:
        """
        source = self.gu.getNode(source_id)
        self.gu.addIndividualToGraph(g, source_id, None)
        node = self.gu.getNode(assoc)
        g.add((node, DC['source'], source))

        return

    def addDescription(self, g, assoc, description):
        node = self.gu.getNode(assoc)
        g.add((node, DC['description'], Literal(description)))
        return

    def loadAllProperties(self, g):
        props = { self.OBJECTPROP: self.object_properties,
              self.ANNOTPROP: self.annotation_properties,
              self.DATAPROP: self.datatype_properties}

        for p in props:
            self.gu.loadProperties(g, props[p], p)

        return

    def addScore(self, g, assoc, score, score_type=None):
        node = self.gu.getNode(assoc)
        has_measurement=self.gu.getNode(self.properties['has_measurement'])
        g.add((node, has_measurement, Literal(score, datatype=XSD['float'])))
        return

    def _process_pub_list(self, pub_list, g, node):
        """
        :param pub_list: list of references with proper prefix
        :return None
        """
        for pub in pub_list:
            source = self._get_source(pub)

            if pub.strip() != '':
                self._add_source_node(g, node, source, pub)

        return

    def _get_source(self, pub_id):
        """
        :param pub_id: Prefixed publication ID
        :return: source: Source IRI containing publication ID
        """
        source = None
        cu = self.cu
        if re.compile('http').match(pub_id):
            source = URIRef(pub_id)
        else:
            u = cu.get_uri(pub_id)
            if u is not None:
                source = URIRef(u)

        return source

    def _add_source_node(self, g, node, source, pub_id):
        """
        :param g: graph
        :param node:
        :param source:
        :param pub_id:
        :return: None
        """
        if source is not None and source != URIRef('[]'):
            g.add((node, DC['source'], source))
            g.add((source, RDF['type'], self.OWLIND))
        else:
            logger.warn("source as a literal -- skipping: %s", pub_id)
            # g.add((node, DC['source'], Literal(pub_id)))

        return

    def make_association_id(self, definedby, subject, predicate, object, attributes=None):
        """
        A method to create unique identifiers for OBAN-style associations, based on all the parts of the assn
        If any of the items is empty or None, it will convert it to blank.
        It effectively md5 hashes the (+)-joined string from the values.
        Subclasses of Assoc can submit an additional array of attributes that will be added to the ID.

        :param definedby: The (data) resource that provided the annotation
        :param subject:
        :param predicate:
        :param object:
        :param evidence:
        :param source:
        :param attributes:
        :return:
        """
        # note others available: md5(), sha1(), sha224(), sha256(), sha384(), and sha512()

        # putting definedby first, as this will usually be the datasource providing the annotation
        # this will end up making the first few parts of the id be the same for all annotations in that resource
        items_to_hash = [definedby, subject, predicate, object]
        if attributes is not None:
            items_to_hash += attributes

        for i, val in enumerate(items_to_hash):
            if val is None:
                items_to_hash[i] = ''

        byte_string = '+'.join(items_to_hash).encode("utf-8")

        return ':'.join(('MONARCH', hashlib.md5(byte_string).hexdigest()))
