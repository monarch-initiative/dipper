__author__ = 'nlw'

import re
import logging

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
        'has_subject': ':hasSubject',
        'has_object': ':hasObject',
        'has_predicate': ':hasPredicate',
        'is_about': 'IAO:00000136',
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


    def __init__(self):
        self.cu = CurieUtil(curie_map.get())
        self.gu = GraphUtils(curie_map.get())
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

    def addAssociationToGraph(self, g):
        cu = self.cu
        gu = self.gu
        # first, add the direct triple
        # anonymous nodes are indicated with underscore
        s = gu.getNode(self.sub)
        o = gu.getNode(self.obj)
        p = gu.getNode(self.rel)

        g.add((s, p, o))

        # now, create the reified relationship with our annotation pattern
        node = gu.getNode(self.annot_id)
        g.add((node, RDF['type'], URIRef(cu.get_uri('Annotation:'))))
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
            g.add((node, DC['evidence'], evidence))

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
            logger.warn("source as a literal -- is this ok?")
            g.add((node, DC['source'], Literal(pub_id)))
            # else:
            #   print("WARN:",self.entity_id,'+',self.phenotype_id,'has no source information for the association (',self.evidence,')')

        return