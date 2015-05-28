__author__ = 'nlw'

from rdflib import RDF
import logging

from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
import re

logger = logging.getLogger(__name__)


class Pathway():
    """
    This provides convenience methods to deal with gene and protein collections in the context of pathways.
    """

    pathway_parts = {
        'signal_transduction' : 'GO:0007165'
    }

    object_properties = {
        'involved_in': 'RO:0002331'
    }

    properties = object_properties.copy()

    def __init__(self, graph):

        self.gu = GraphUtils(curie_map.get())

        self.graph = graph

        self.gu.loadProperties(self.graph, self.object_properties, self.gu.OBJPROP)

        return

    def addPathway(self, pathway_id, pathway_label, pathway_type=None, pathway_description=None):
        """
        Adds a pathway as a class.  If no specific type is specified, it will
        default to a subclass of "signal transduction".
        :param pathway_id:
        :param pathway_label:
        :param pathway_type:
        :param pathway_description:
        :return:
        """
        if pathway_type is None:
            pathway_type = self.pathway_parts['signal_transduction']
        self.gu.addClassToGraph(self.graph, pathway_id, pathway_label, pathway_type, pathway_description)
        return

    def addGeneToPathway(self, pathway_id, gene_id):
        """
        gene_id involved_in pathway_id

        :param pathway_id:
        :param gene_id:
        :return:
        """
        # TODO consider converting this and creating a BNode of a "gene product of gene_id"
        self.gu.addTriple(self.graph, gene_id, self.object_properties['involved_in'], pathway_id)
        return

    def addProteinToPathway(self, pathway_id, protein_id):
        self.gu.addTriple(self.graph, protein_id, self.object_properties['involved_in'], pathway_id)

        return