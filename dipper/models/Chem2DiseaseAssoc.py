from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper.curie import curie_map


class Chem2DiseaseAssoc(Assoc):

    def __init__(self, assoc_id, chem_id, phenotype_id, pub_list, rel, evidence):
        super().__init__()
        self.annot_id = assoc_id
        self.chem_id = chem_id
        self.phenotype_id = phenotype_id
        self.pub_list = pub_list
        self.pub_id = None
        self.evidence = evidence
        self.rel = rel
        self.cu = CurieUtil(curie_map.get())

        self.setSubject(chem_id)
        self.setObject(phenotype_id)
        return

    def addAssociationNodeToGraph(self, g):
        '''
        The reified relationship between a genotype (or any genotype part) and a phenotype
        is decorated with some provenance information.
        This makes the assumption that both the genotype and phenotype are classes.

        currently hardcoded to map the annotation to the monarch namespace
        :param g:
        :return:
        '''

        self.addAssociationToGraph(g)

        return g