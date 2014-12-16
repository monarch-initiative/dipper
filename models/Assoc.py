__author__ = 'nicole'

from rdflib import Namespace, OWL, RDF, URIRef


class Assoc:
    ###### Namespaces ######
    # TODO turn this into a dictionary from the context file
    curie_map = {
        'MONARCH': 'http://www.monarchinitiative.org/MONARCH_',
        'HP': 'http://purl.obolibrary.org/obo/HP_',
        'ECO': 'http://purl.obolibrary.org/obo/ECO_',
        'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/PMID_',
        'RO' : 'http://purl.obolibrary.org/obo/RO_',
        'OBO' : 'http://purl.obolibrary.org/obo/'
    }

    relationships = {
        'has_disposition':'http://purl.obolibrary.org/obo/GENO_0000208',
        'has_phenotype':'http://purl.obolibrary.org/obo/RO_0002200'
    }


    OWLCLASS=OWL['Class']
    OWLIND=OWL['NamedIndividual']
    OWLANNOT=OWL['Annotation']
    OWLPROP=OWL['ObjectProperty']
    OWLSUBCLASS=OWL['subclassOf']


    def __init__(self):
        return

    def get_namespaces(self):
        return self.curie_map

    def get_relationships(self):
        return self.relationships

    def createAssociationNode(self, g):
        #todo make a general association object following our pattern

        return

    def loadObjectProperties(self,g):
        for k in self.relationships:
            g.add((URIRef(self.relationships[k]),RDF['type'],self.OWLPROP))
        return