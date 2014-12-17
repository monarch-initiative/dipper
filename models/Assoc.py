__author__ = 'nicole'

from rdflib import Namespace, OWL, RDF, URIRef


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
        'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/PMID_',
        'RO' : 'http://purl.obolibrary.org/obo/RO_',
        'OBO' : 'http://purl.obolibrary.org/obo/',
        'Annotation' : 'http://www.w3.org/ns/oa#Annotation',
        '' : 'http://www.monarchinitiative.org/'  #base
    }

    relationships = {
        'has_disposition':'http://purl.obolibrary.org/obo/GENO_0000208',
        'has_phenotype':'http://purl.obolibrary.org/obo/RO_0002200'
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
        for k in self.relationships:
            g.add((URIRef(self.relationships[k]),RDF['type'],self.OWLPROP))
        return