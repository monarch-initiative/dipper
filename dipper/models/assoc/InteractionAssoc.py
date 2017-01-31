from dipper.models.assoc.Association import Assoc

__author__ = 'nlw'


class InteractionAssoc(Assoc):

    interaction_object_properties = {
        'genetically_interacts_with': 'RO:0002435',
        # using for directly interacts with.  better choice? psi-mi:"MI:0407"
        'interacts_with': 'RO:0002434',
        # should we use this instead for direct interaction?
        'molecularly_interacts_with': 'RO:0002436',
        # psi-mi:"MI:0403"(colocalization)
        'colocalizes_with': 'RO:0002325',
        'ubiquitinates': 'RO:0002480',
        'regulates': 'RO:0002448',
        'positively_regulates': 'RO:0003003',
        'negatively_regulates': 'RO:0003002',
    }

    def __init__(self, graph, definedby, subj, obj, rel=None):
        super().__init__(graph, definedby)

        self.set_subject(subj)
        self.set_object(obj)
        if rel is None:
            rel = self.rel['interacts_with']
        self.set_relationship(rel)

        return
