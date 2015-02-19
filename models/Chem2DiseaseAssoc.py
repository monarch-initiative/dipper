from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil
from conf import curie_map


class Chem2DiseaseAssoc(Assoc):
    """
    A specific association class for defining Chemical-to-Phenotype relationships
    This assumes that a graph is created outside of this class, and nodes get added.
    By default, an association will assume the "has_phenotype" relationship, unless
    otherwise specified.
    Note that genotypes are expected to be created and defined outside of this association,
    most likely by calling methods in the Genotype() class.
    """

    def __init__(self, assoc_id, chem_id, phenotype_id, pub_list, evidence_code):
        self.annot_id = assoc_id
        self.chem_id = chem_id
        self.phenotype_id = phenotype_id
        self.pub_list = pub_list
        self.pub_id = None
        self.evidence = evidence_code
        self.rel = self.relationships['has_phenotype']  # default to has_phenotype
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