__author__ = 'nlw'


from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Assoc import Assoc
from dipper import curie_map


class OrthologyAssoc(Assoc):

    ortho_rel = {
        'orthologous': 'RO:HOM0000017',  # in orthology relationship with
        'least_diverged_orthologous': 'RO:HOM0000020',  # in 1 to 1 orthology relationship with
        'homologous': 'RO:HOM0000019',  # in 1 to 1 homology relationship with
        'paralogous': 'RO:HOM0000011',  # in paralogy relationship with (generic)
        'in_paralogous': 'RO:HOM0000023',  # in in-paralogy relationship with
        'ohnologous': 'RO:HOM0000022',  # in ohnology relationship with
        'xenologous': 'RO:HOM0000018',  # in xenology relationship with
        'has_member': 'RO:0002351'
        }

    def __init__(self, assoc_id, gene1, gene2, pub, evidence_code):
        super().__init__()
        self.cu = CurieUtil(curie_map.get())
        self.gu = GraphUtils(curie_map.get())
        self.object_properties.update(self.ortho_rel)
        self.properties.update(self.ortho_rel)
        self.annot_id = assoc_id
        self.gene1 = gene1
        self.gene2 = gene2
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.properties['orthologous']  # default
        self.pub_list = None

        self.setSubject(gene1)
        self.setObject(gene2)

        return

    def addOrthologyAssociationToGraph(self, g):

        self.addAssociationToGraph(g)

        return

    def addGeneFamilyToGraph(self, g, family_id):
        """
        Make an association between a group of genes and some grouping class.
        We make the assumption that the genes in the association are part of the supplied
        family_id, and that the genes have already been declared as classes elsewhere.
        The family_id is added as an individual of type DATA:gene_family.

        Triples:
        <family_id> a DATA:gene_family
        <family_id> RO:has_member <gene1>
        <family_id> RO:has_member <gene2>

        :param family_id:
        :param g: the graph to modify
        :return:
        """
        gene_family = 'DATA:3148'  # http://edamontology.org/data_3148

        # make the assumption that the genes have already been added as classes previously
        self.gu.addIndividualToGraph(g, family_id, None, gene_family)

        # add each gene to the family
        self.gu.addMember(g, family_id, self.gene1)
        self.gu.addMember(g, family_id, self.gene2)

        return
