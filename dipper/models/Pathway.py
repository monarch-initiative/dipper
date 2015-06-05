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
        'signal_transduction' : 'GO:0007165',
        'cellular_process' : 'GO:0009987',
        'pathway' : 'PW:0000001'
    }

    object_properties = {
        'involved_in': 'RO:0002331',
        'gene_product_of' : 'RO:0002204',
        'has_gene_product' : 'RO:0002205'
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
        default to a subclass of "GO:cellular_process" and "PW:pathway".
        :param pathway_id:
        :param pathway_label:
        :param pathway_type:
        :param pathway_description:
        :return:
        """
        if pathway_type is None:
            pathway_type = self.pathway_parts['cellular_process']
        self.gu.addClassToGraph(self.graph, pathway_id, pathway_label, pathway_type, pathway_description)
        self.gu.addSubclass(self.graph, self.pathway_parts['pathway'], pathway_id)

        return

    def addGeneToPathway(self, pathway_id, gene_id):
        """
        gene_id involved_in pathway_id

        :param pathway_id:
        :param gene_id:
        :return:
        """
        gene_product = '_'+gene_id+'product'
        # FIXME figure out what the type of the gene product is (not necessarily a protein)
        self.gu.addClassToGraph(self.graph, gene_product, None)
        #self.gu.addTriple(self.graph, gene_product, self.object_properties['gene_product_of'], gene_id)
        self.gu.addTriple(self.graph, gene_id, self.object_properties['has_gene_product'], gene_product)
        self.addComponentToPathway(pathway_id, gene_product)

        return

    def addComponentToPathway(self, pathway_id, component_id):
        self.gu.addTriple(self.graph, component_id, self.object_properties['involved_in'], pathway_id)

        return