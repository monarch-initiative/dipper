# Right now just a copy of InteractionAssoc
from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
from models.Assoc import Assoc
import curie_map

class Gene2Pathway(Assoc):

    relationship_map = {
        'genetically_interacts_with' : 'RO:0002435',
        'interacts_with' : 'RO:0002434',  #use this for directly interacts with.  better choice? psi-mi:"MI:0407"(direct interaction)
        'molecularly_interacts_with' : 'RO:0002436',  #should we use this instead for direct interaction?
        'colocalizes_with' : 'RO:0002325', #psi-mi:"MI:0403"(colocalization)
        'ubiquitinates' : 'RO:0002480'
    }

    def __init__(self,assoc_id, pathway_id, gene_id, evidence_code, pathway_label=None, gene_label=None):
        self.cu = CurieUtil(curie_map.get())
        self.annot_id = assoc_id
        self.subj = pathway_id
        self.obj = gene_id
        self.pub_id = None
        self.evidence = evidence_code
        self.rel = self.relationship_map['interacts_with']  # default
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