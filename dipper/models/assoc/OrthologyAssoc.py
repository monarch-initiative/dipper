from dipper.models.assoc.Association import Assoc
from dipper.models.Family import Family
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

__author__ = 'nlw'


class OrthologyAssoc(Assoc):

    def __init__(
            self,
            graph,
            definedby,
            gene1,
            gene2,
            rel=None,
            subject_category=None,
            object_category=None
    ):
        super().__init__(graph, definedby)
        self.globaltt = graph.globaltt
        self.globaltcid = graph.globaltcid
        self.curie_map = graph.curie_map

        if rel is None:
            rel = self.globaltt['in orthology relationship with']  # default

        self.set_subject(gene1)
        self.set_object(gene2)
        self.subject_category = subject_category
        self.object_category = object_category
        self.set_relationship(rel)

    def add_gene_family_to_graph(self, family_id):
        """
        Make an association between a group of genes and some grouping class.
        We make the assumption that the genes in the association
        are part of the supplied family_id, and that the genes have
        already been declared as classes elsewhere.
        The family_id is added as an individual of type DATA:gene_family.

        Triples:
        <family_id> a EDAM-DATA:gene_family
        <family_id> RO:has_member <gene1>
        <family_id> RO:has_member <gene2>
        <gene1> biolink:category <subject_category>
        <gene2> biolink:category <object_category>
        :param family_id:
        :param g: the graph to modify
        :return:
        """
        family = Family(self.graph)
        gene_family = self.globaltt['gene_family']

        # make the assumption that the genes
        # have already been added as classes previously
        self.model.addIndividualToGraph(
            family_id, None, gene_family, ind_category=blv.terms['GeneFamily']
        )

        # add each gene to the family
        family.addMember(
            family_id,
            self.sub,
            group_category=blv.terms['GeneFamily'],
            member_category=self.subject_category
        )

        family.addMember(
            family_id,
            self.obj,
            group_category=blv.terms['GeneFamily'],
            member_category=self.object_category
        )
