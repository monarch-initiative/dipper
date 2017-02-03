from dipper.graph.Graph import Graph


class Family():
    """
    Model mereological/part whole relationships

    Although these relations are more abstract, we often
    use them to model family relationships (proteins, humans, etc.)
    The naming of this class may change in the future to better
    reflect the meaning of the relations it is modeling
    """

    object_properties = {
        # has member is a mereological relation
        # between a collection and an item
        'has_member': 'RO:0002351',

        # is member of is a mereological relation
        # between a item and a collection
        'member_of': 'RO:0002350'
    }

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)

    def addMember(self, group_id, member_id):
        self.graph.addTriple(
            group_id, self.object_properties['has_member'], member_id)
        return

    def addMemberOf(self, member_id, group_id):
        self.graph.addTriple(
            member_id, self.object_properties['member_of'], group_id)
        return