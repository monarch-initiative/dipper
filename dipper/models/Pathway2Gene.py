# Right now just a copy of InteractionAssoc
from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Assoc import Assoc
from dipper import curie_map


class Pathway2Gene(Assoc):

    relationship_map = {
        'has_member': 'RO:0002351'
    }

    def __init__(self,assoc_id, pathway_id, gene_id, evidence_code, pathway_label=None, gene_label=None):
        super().__init__()
        self.cu = CurieUtil(curie_map.get())
        self.annot_id = assoc_id
        self.subj = pathway_id
        self.obj = gene_id
        self.pub_id = None
        self.evidence = evidence_code
        self.rel = self.relationship_map['has_member']  # default
        self.cu = CurieUtil(curie_map.get())
        self.pub_list = None
        self.pathway_label = pathway_label
        self.gene_label = gene_label

        self.setSubject(pathway_id)
        self.setObject(gene_id)

        return

    def addInteractionAssociationToGraph(self,graph):

        self.addAssociationToGraph(graph)

        #todo add some other stuff related to the experimental methods?

        return

    def add_gene_as_class(self, graph):
        gu = GraphUtils(curie_map.get())
        gu.addClassToGraph(graph, self.obj, self.gene_label)

        return

    def add_pathway_as_class(self, graph):
        gu = GraphUtils(curie_map.get())
        gu.addClassToGraph(graph, self.subj, self.pathway_label)

        return