__author__ = 'nicole'

from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Assoc import Assoc
from dipper import curie_map


class InteractionAssoc(Assoc):

    rel = {
        'genetically_interacts_with': 'RO:0002435',
        'interacts_with': 'RO:0002434',  # using for directly interacts with.  better choice? psi-mi:"MI:0407"
        'molecularly_interacts_with': 'RO:0002436',  # should we use this instead for direct interaction?
        'colocalizes_with': 'RO:0002325',  # psi-mi:"MI:0403"(colocalization)
        'ubiquitinates': 'RO:0002480'
    }

    def __init__(self, definedby, subj, obj, rel=None):
        super().__init__(definedby)
        self.pub_list = None
        self.setSubject(subj)
        self.setObject(obj)
        if rel is None:
            rel = self.rel['interacts_with']
        self.setRelationship(rel)

        return

    def loadAllProperties(self, g):

        super().loadAllProperties(g)
        self.gu.loadObjectProperties(g, self.rel)

        return