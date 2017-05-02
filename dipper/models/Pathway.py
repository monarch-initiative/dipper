import logging
import re
from dipper.models.Model import Model
from dipper.graph.Graph import Graph

__author__ = 'nlw'

logger = logging.getLogger(__name__)


class Pathway():
    """
    This provides convenience methods to deal with gene and protein collections
    in the context of pathways.
    """

    pathway_parts = {
        'signal_transduction': 'GO:0007165',
        'cellular_process': 'GO:0009987',
        'pathway': 'PW:0000001',
        'gene_product': 'CHEBI:33695'  # bioinformation molecule
    }

    object_properties = {
        'involved_in': 'RO:0002331',
        'gene_product_of': 'RO:0002204',
        'has_gene_product': 'RO:0002205'
    }

    properties = object_properties.copy()

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)

        return

    def addPathway(
            self, pathway_id, pathway_label, pathway_type=None,
            pathway_description=None):
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
        self.model.addClassToGraph(
            pathway_id, pathway_label, pathway_type, pathway_description)
        self.model.addSubClass(pathway_id, self.pathway_parts['pathway'])

        return

    def addGeneToPathway(self, gene_id, pathway_id):
        """
        When adding a gene to a pathway, we create an intermediate
        'gene product' that is involved in
        the pathway, through a blank node.

        gene_id RO:has_gene_product _gene_product
        _gene_product RO:involved_in pathway_id

        :param pathway_id:
        :param gene_id:
        :return:
        """

        gene_product = '_:'+re.sub(r':', '', gene_id)+'product'
        self.model.addIndividualToGraph(
            gene_product, None, self.pathway_parts['gene_product'])
        self.graph.addTriple(
            gene_id, self.object_properties['has_gene_product'], gene_product)
        self.addComponentToPathway(pathway_id, gene_product)

        return

    def addComponentToPathway(self, component_id, pathway_id):
        """
        This can be used directly when the component is directly involved in
        the pathway.  If a transforming event is performed on the component
        first, then the addGeneToPathway should be used instead.

        :param pathway_id:
        :param component_id:
        :return:
        """
        self.graph.addTriple(
            component_id, self.object_properties['involved_in'], pathway_id)

        return
