__author__ = 'nlw'

import logging

from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map


logger = logging.getLogger(__name__)


class Environment():
    """
    These methods provide convenient methods to add items related to an experimental environment
    and it's parts to a supplied graph.

    This is a stub ready for expansion.
    """

    # special genotype parts mapped to their GENO and SO classes that we explicitly reference here
    environment_parts = {
        'environmental_system': 'ENVO:01000254',
        'environmental_condition' : 'XCO:0000000',
        'morpholio_reagent': 'REO:0000042',
        'talen_reagent': 'REO:0001022',
        'crispr_reagent': 'REO:crispr_TBD'
    }

    object_properties = {
        'has_part': 'BFO:0000051',
    }

    annotation_properties = {
    }

    properties = object_properties.copy()
    properties.update(annotation_properties)

    def __init__(self, graph):

        self.gu = GraphUtils(curie_map.get())

        self.graph = graph

        self.gu.loadProperties(self.graph, self.object_properties, self.gu.OBJPROP)

        return

    def addEnvironment(self, env_id, env_label, env_type=None, env_description=None):
        if env_type is None:
            env_type = self.environment_parts['environmental_system']

        self.gu.addIndividualToGraph(self.graph, env_id, env_label, env_type, env_description)

        return

    def addEnvironmentalCondition(self, cond_id, cond_label, cond_type=None, cond_description=None):
        if cond_type is None:
            cond_type = self.environment_parts['environmental_condition']

        self.gu.addIndividualToGraph(self.graph, cond_id, cond_label, cond_type, cond_description)

        return

    def addComponentToEnvironment(self, env_id, component_id):

        self.gu.addTriple(self.graph, env_id, self.gu.object_properties['has_part'], component_id)

        return

    def addComponentAttributes(self, component_id, entity_id, value=None, unit=None):

        self.gu.addTriple(self.graph, component_id, self.gu.object_properties['has_part'], entity_id)
        # TODO add value and units

        return
