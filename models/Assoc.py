__author__ = 'nicole'

from rdflib import Namespace, URIRef, Literal
from rdflib.namespace import RDF,DC,OWL
from utils.CurieUtil import CurieUtil
import re


class Assoc:
    '''
    An abstract class for Monarch-style associations, to enable attribution of source and evidence
    on statements.
    '''

    # TODO turn this into a dictionary from the context file
    curie_map = {
        'MONARCH': 'http://www.monarchinitiative.org/MONARCH_',
        'HP': 'http://purl.obolibrary.org/obo/HP_',
        'ECO': 'http://purl.obolibrary.org/obo/ECO_',
        'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/',
        'ISBN-10' : 'http://www.monarchinitiative.org/ISBN_',
        'ISBN-13' : 'http://www.monarchinitiative.org/ISBN_',
        'ISBN-15' : 'http://www.monarchinitiative.org/ISBN_',
        'ISBN' : 'http://www.monarchinitiative.org/ISBN_',
        'RO' : 'http://purl.obolibrary.org/obo/RO_',
        'GENO' : 'http://purl.obolibrary.org/obo/GENO_',
        'OBO' : 'http://purl.obolibrary.org/obo/',
        'Annotation' : 'http://www.w3.org/ns/oa#Annotation',
        '' : 'http://www.monarchinitiative.org/'  #base
    }

    relationships = {
        'has_disposition':'GENO:0000208',
        'has_phenotype':'RO:0002200'
    }

    OWLCLASS=OWL['Class']
    OWLIND=OWL['NamedIndividual']
    OWLPROP=OWL['ObjectProperty']
    OWLSUBCLASS=OWL['subclassOf']
    BASE=Namespace(curie_map[''])


    def __init__(self):
        return

    def get_namespaces(self):
        return self.curie_map

    def get_relationships(self):
        return self.relationships

    def createAssociationNode(self, g):
        #TODO make a general association object following our pattern

        return

    def addAssociationToGraph(self,g):
        namespaces = self.curie_map

        #first, add the direct triple
        s = URIRef(self.cu.get_uri(self.sub))
        p = URIRef(self.cu.get_uri(self.rel))
        o = URIRef(self.cu.get_uri(self.obj))

        g.add((s, RDF['type'], self.OWLCLASS))
        g.add((o, RDF['type'], self.OWLCLASS))
        g.add((s, p, o))

        #now, create the reified relationship with our annotation pattern
        node = URIRef(self.cu.get_uri(self.annot_id))
        g.add((node, RDF['type'], URIRef(self.cu.get_uri('Annotation:'))))
        g.add((node, self.BASE['hasSubject'], s))
        g.add((node, self.BASE['hasObject'], o))

        #this is handling the occasional messy pubs that are sometimes literals
        if (self.pub_id is not None):
            if (re.compile('http').match(self.pub_id)):
                source = URIRef(self.pub_id)
            else:
                source = URIRef(self.cu.get_uri(self.pub_id))

        evidence = URIRef(self.cu.get_uri(self.evidence))
        if (self.pub_id is not None and self.pub_id.strip() != ''):
            if (source != URIRef('[]')):
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
        cu = CurieUtil(self.curie_map)
        for k in self.relationships:
            g.add((URIRef(cu.get_uri(self.relationships[k])),RDF['type'],self.OWLPROP))
        return