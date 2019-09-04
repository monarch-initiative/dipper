import logging
import re
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv
from dipper.graph.Graph import Graph

__author__ = 'nlw'

LOG = logging.getLogger(__name__)
# note: currently no log issued


class Pathway():
    """
    This provides convenience methods to deal with gene and protein collections
    in the context of pathways.
    """

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".format(graph))
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
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
            pathway_type = self.globaltt['cellular_process']
        self.model.addClassToGraph(
            pathway_id, pathway_label, pathway_type, pathway_description,
            class_category=blv.terms.Pathway.value)
        self.model.addSubClass(pathway_id, self.globaltt['pathway'],
                               child_category=blv.terms.Pathway.value)

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

        gene_product = '_:'+re.sub(r':', '', gene_id) + 'product'
        self.model.addIndividualToGraph(
            gene_product, None, self.globaltt['gene_product'],
            ind_category=blv.terms.GeneProduct.value)
        self.graph.addTriple(
            gene_id, self.globaltt['has gene product'], gene_product,
            subject_category=blv.terms.Gene.value,
            object_category=blv.terms.GeneProduct.value)
        self.addComponentToPathway(gene_product, pathway_id)

        return

    def addComponentToPathway(self, component_id, pathway_id,
                              component_category=None,
                              pathway_category=blv.terms.Pathway.value):
        """
        This can be used directly when the component is directly involved in
        the pathway.  If a transforming event is performed on the component
        first, then the addGeneToPathway should be used instead.

        :param pathway_id:
        :param component_id:
        :param component_category: biolink category for component_id
        :param pathway_category: biolink category for pathway_id
        :return:
        """
        self.graph.addTriple(component_id, self.globaltt['involved in'], pathway_id,
                             subject_category=component_category,
                             object_category=pathway_category)

        return
