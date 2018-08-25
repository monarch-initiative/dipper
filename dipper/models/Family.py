from dipper.graph.Graph import Graph


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
            raise ValueError("{} is not a graph".graph)
        self.model = Model(graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map

    def addMember(self, group_id, member_id):
        self.graph.addTriple(group_id, self.globaltt['has member'], member_id)
        return

    def addMemberOf(self, member_id, group_id):
        self.graph.addTriple(member_id, self.globaltt['member of'], group_id)
        return
