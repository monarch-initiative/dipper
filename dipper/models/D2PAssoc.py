__author__ = 'nicole'

from rdflib import Namespace, Literal

from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper.curie import curie_map


# This one is specific for making a disease-to-phenotype
class D2PAssoc(Assoc):
    '''
    A specific association class for defining Disease-to-Phenotype relationships
    This assumes that a graph is created outside of this class, and nodes get added.
    By default, an association will assume the "has_phenotype" relationship, unless
    otherwise specified.


    '''

    def __init__(self, assoc_id, entity_id, phenotype_id, onset, frequency, pub, evidence_code):
        self.annot_id = assoc_id
        self.entity_id = entity_id
        self.phenotype_id = phenotype_id
        self.onset = onset
        self.frequency = frequency
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.relationships['has_phenotype']
        self.cu = CurieUtil(curie_map.get())
        self.pub_list = None

        self.setSubject(entity_id)
        self.setObject(phenotype_id)

        return

    def set_relationship(self, rel):
        self.rel = rel
        return

    def addAssociationNodeToGraph(self, g):
        '''
        The reified relationship between a disease and a phenotype is decorated with some provenance information.
        This makes the assumption that both the disease and phenotype are classes.

        currently hardcoded to map the annotation to the monarch namespace
        :param g:
        :return:
        '''
        namespaces = curie_map.get()

        #add the basic association nodes
        self.addAssociationToGraph(g)

        #add the specific attributes for this association type
        n = Namespace(namespaces['MONARCH'])
        node = n[self.annot_id]
        if (self.frequency is not None and self.frequency != ''):
            #FIXME what is the real predicate here?
            g.add((node, self.BASE['frequencyOfPhenotype'], Literal(self.frequency)))
        if (self.onset is not None and self.onset != ''):
            #FIXME what is the real predicate here?
            g.add((node, self.BASE['onset'], n[self.onset]))

        return g

