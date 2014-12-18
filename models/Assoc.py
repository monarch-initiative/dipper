__author__ = 'nicole'

from rdflib import Namespace, OWL, RDF, URIRef
from utils.CurieUtil import CurieUtil


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