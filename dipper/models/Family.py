from dipper.graph.Graph import Graph
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

class Family():
    """
    Model mereological/part whole relationships

    Although these relations are more abstract, we often
    use them to model family relationships (proteins, humans, etc.)
    The naming of this class may change in the future to better
    reflect the meaning of the relations it is modeling
    """

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".format(graph))
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map

    def addMember(self, group_id, member_id,
                  group_category=None, member_category=None):
        self.graph.addTriple(group_id, self.globaltt['has member'], member_id,
                             subject_category=group_category,
                             object_category=member_category)
        return

    def addMemberOf(self, member_id, group_id,
                    group_category=None, member_category=None):
        self.graph.addTriple(member_id, self.globaltt['member of'], group_id,
                             subject_category=group_category,
                             object_category=member_category)
        return
