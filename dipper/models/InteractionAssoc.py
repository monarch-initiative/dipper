__author__ = 'nicole'

from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Assoc import Assoc
from dipper import curie_map


class InteractionAssoc(Assoc):

    rel = {
        'genetically_interacts_with' : 'RO:0002435',
        'interacts_with' : 'RO:0002434',  #use this for directly interacts with.  better choice? psi-mi:"MI:0407"(direct interaction)
        'molecularly_interacts_with' : 'RO:0002436',  #should we use this instead for direct interaction?
        'colocalizes_with' : 'RO:0002325', #psi-mi:"MI:0403"(colocalization)
        'ubiquitinates' : 'RO:0002480'
    }

    def __init__(self,assoc_id, subj, obj, pub, evidence_code):
        super()
        self.object_properties.update(self.rel)
        self.properties.update(self.rel)
        self.cu = CurieUtil(curie_map.get())
        self.gu = GraphUtils(curie_map.get())
        self.annot_id = assoc_id
        self.subj = subj
        self.obj = obj
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.rel['interacts_with']  # default
        self.cu = CurieUtil(curie_map.get())
        self.pub_list = None

        self.setSubject(subj)
        self.setObject(obj)


        return

    def addInteractionAssociationToGraph(self,g):

        self.addAssociationToGraph(g)

        #todo add some other stuff related to the experimental methods?

        return

