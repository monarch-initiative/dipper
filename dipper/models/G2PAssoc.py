__author__ = 'nicole'

from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper.curie import curie_map


class G2PAssoc(Assoc):
    '''
    A specific association class for defining Genotype-to-Phenotype relationships
    This assumes that a graph is created outside of this class, and nodes get added.
    By default, an association will assume the "has_phenotype" relationship, unless
    otherwise specified.
    Note that genotypes are expected to be created and defined outside of this association,
    most likely by calling methods in the Genotype() class.
    '''

    def __init__(self, assoc_id, entity_id, phenotype_id, pub, evidence_code):
        super().__init__()
        self.annot_id = assoc_id
        self.entity_id = entity_id
        self.phenotype_id = phenotype_id
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.relationships['has_phenotype']  # default to has_phenotype
        self.cu = CurieUtil(curie_map.get())
#        self.gu = GraphUtils(curie_map.get())
        self.pub_list = None

        self.setSubject(entity_id)
        self.setObject(phenotype_id)
        return

    def set_relationship(self, rel):
        self.rel = rel
        return

    def set_stage(self, start_stage_id, end_stage_id):
        self.start_stage_id = start_stage_id
        self.end_stage_id = end_stage_id
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

        #TODO add staging information here
        #if (self.start_stage_id is not None):
        #    g.add((node, self.BASE['during'], URIRef(self.cu.get_uri(self.start_stage_id))))

        return g

