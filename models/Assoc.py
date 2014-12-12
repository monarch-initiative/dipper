__author__ = 'nicole'

from rdflib import Namespace


class Assoc:
    ###### Namespaces ######
    # TODO turn this into a dictionary from the context file
    namespaces = {
        'MONARCH': 'http://www.monarchinitiative.org/',
        'HP': 'http://purl.obofoundry.org/obo/HP_',
        'ECO': 'http://purl.obofoundry.org/obo/ECO_',
        'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/'
    }

    def __init__(self):
        return

    def get_namespaces(self):
        return self.namespaces

    def createAssociationNode(self, g):
        #todo make a general association object following our pattern

        return