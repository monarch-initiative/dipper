import logging
from dipper.models.Model import Model
from dipper.graph.Graph import Graph

__author__ = 'nlw'

logger = logging.getLogger(__name__)


class Environment():
    """
    These methods provide convenient methods
    to add items related to an experimental environment
    and it's parts to a supplied graph.

    This is a stub ready for expansion.
    """


    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        return

    def addEnvironment(
            self, env_id, env_label, env_type=None, env_description=None):
        if env_type is None:
            env_type = self.globaltt['environmental_system']

        self.model.addIndividualToGraph(
            env_id, env_label, env_type, env_description)

        return

    def addEnvironmentalCondition(
            self, cond_id, cond_label, cond_type=None, cond_description=None):
        if cond_type is None:
            cond_type = self.globaltt['environmental_condition']

        self.model.addIndividualToGraph(
            cond_id, cond_label, cond_type, cond_description)

        return

    def addComponentToEnvironment(self, env_id, component_id):

        self.graph.addTriple(env_id, self.globaltt['has_part'], component_id)

        return

    def addComponentAttributes(self, component_id, entity_id, value=None, unit=None):

        self.graph.addTriple(
            component_id, self.globaltt['has_part'], entity_id)
        # TODO add value and units

        return
