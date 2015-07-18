__author__ = 'nlw'

from dipper.models.assoc.Association import Assoc


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

        self.set_subject(subj)
        self.set_object(obj)
        if rel is None:
            rel = self.rel['interacts_with']
        self.set_relationship(rel)

        return

    def load_all_properties(self, g):

        super().load_all_properties(g)
        self.gu.loadObjectProperties(g, self.rel)

        return