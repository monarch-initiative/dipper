__author__ = 'nicole'

from rdflib import Namespace, URIRef, Literal
from rdflib.namespace import RDF,DC,OWL,RDFS
from utils.CurieUtil import CurieUtil
import re
import curie_map


class Assoc:
    '''
    An abstract class for Monarch-style associations, to enable attribution of source and evidence
    on statements.
    '''

    relationships = {
        'has_disposition':'GENO:0000208',
        'has_phenotype':'RO:0002200',
        'replaced_by' : 'IAO:0100001',
        'consider' : 'OIO:consider',
        'hasExactSynonym' : 'OIO:hasExactSynonym',
        'hasRelatedSynonym' : 'OIO:hasRelatedSynonym',
        'definition' : 'IAO:0000115',
        'in_taxon' : 'RO:0002162'
    }

    OWLCLASS=OWL['Class']
    OWLIND=OWL['NamedIndividual']
    OWLPROP=OWL['ObjectProperty']
    SUBCLASS=RDFS['subClassOf']
    BASE=Namespace(curie_map.get()[''])


    def __init__(self):
        self.cu = CurieUtil(curie_map.get())
        return

    def get_namespaces(self):
        if (self.namespaces):
            return self.namespaces

        return None

    def get_relationships(self):
        return self.relationships

    def createAssociationNode(self, g):
        #TODO make a general association object following our pattern

        return

    def addAssociationToGraph(self,g):
        cu = self.cu

        #first, add the direct triple
        s = URIRef(cu.get_uri(self.sub))
        p = URIRef(cu.get_uri(self.rel))
        o = URIRef(cu.get_uri(self.obj))

        g.add((s, RDF['type'], self.OWLCLASS))
        g.add((o, RDF['type'], self.OWLCLASS))
        g.add((s, p, o))

        #now, create the reified relationship with our annotation pattern
        node = URIRef(cu.get_uri(self.annot_id))
        g.add((node, RDF['type'], URIRef(cu.get_uri('Annotation:'))))
        g.add((node, self.BASE['hasSubject'], s))
        g.add((node, self.BASE['hasObject'], o))

        #this is handling the occasional messy pubs that are sometimes literals
        if (self.pub_id is not None):
            if (re.compile('http').match(self.pub_id)):
                source = URIRef(self.pub_id)
            else:
                u = cu.get_uri(self.pub_id)
                if (u is not None):
                    source = URIRef(u)
                else:
                    source = None

        evidence = URIRef(cu.get_uri(self.evidence))
        if (self.pub_id is not None and self.pub_id.strip() != ''):
            if (source is not None and source != URIRef('[]')):
                g.add((node, DC['source'], source))
                g.add((source, RDF['type'], self.OWLIND))
            else:
                print("WARN: source as a literal -- is this ok?")
                g.add((node, DC['source'], Literal(self.pub_id)))
            #            else:
            #                print("WARN:",self.entity_id,'+',self.phenotype_id,'has no source information for the association (',self.evidence,')')

        if (self.evidence is None or self.evidence.strip() == ''):
            print("WARN:", self.sub, '+', self.obj, 'has no evidence code')
        else:
            g.add((node, DC['evidence'], evidence))

        return

    def setSubject(self,identifier):
        self.sub = identifier
        return

    def setObject(self,identifier):
        self.obj = identifier
        return

    def setRelationship(self,identifier):
        self.rel = identifier
        return

    def loadObjectProperties(self,g):
        '''
        Given a graph, it will load any of the items in relationships dictionary
        as owl['ObjectProperty'] types
        :param g: a graph
        :return: None
        '''
        cu = self.cu
        for k in self.relationships:
            g.add((URIRef(cu.get_uri(self.relationships[k])),RDF['type'],self.OWLPROP))
        return