from dipper.models.assoc.Association import Assoc

__author__ = 'nlw'


class InteractionAssoc(Assoc):

    def __init__(self, graph, definedby, subj, obj, rel=None):
        super().__init__(graph, definedby)

        self.set_subject(subj)
        self.set_object(obj)
        if rel is None:
            rel = self.globaltt['interacts with']
        self.set_relationship(rel)
