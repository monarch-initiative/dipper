import logging

from dipper.models.assoc.Association import Assoc

__author__ = 'nlw'

LOG = logging.getLogger(__name__)
# note: currently no log issued


class G2PAssoc(Assoc):
    """
    A specific association class for defining Genotype-to-Phenotype
    relationships. This assumes that a graph is created outside of this class,
    and nodes get added.
    By default, an association will assume the "has_phenotype" relationship,
    unless otherwise specified.
    Note that genotypes are expected to be
    created and defined outside of this association,
    most likely by calling methods in the Genotype() class.

    """

    def __init__(
            self,
            graph,
            definedby,
            entity_id,
            phenotype_id,
            rel=None,
            entity_category=None,
            phenotype_category=None
    ):
        super().__init__(graph, definedby)
        self.entity_id = entity_id
        self.phenotype_id = phenotype_id

        if rel is None:
            rel = self.globaltt['has phenotype']

        self.start_stage_id = None
        self.end_stage_id = None
        self.environment_id = None
        self.stage_process_id = None

        self.set_subject(entity_id)
        self.set_object(phenotype_id)
        self.set_relationship(rel)

        self.subject_category = entity_category
        self.object_category = phenotype_category

        return

    def set_stage(self, start_stage_id, end_stage_id):
        if start_stage_id is not None and start_stage_id.strip() != '':
            self.start_stage_id = start_stage_id
        if end_stage_id is not None and end_stage_id.strip() != '':
            self.end_stage_id = end_stage_id

    def set_environment(self, environment_id):
        if environment_id is not None and environment_id.strip() != '':
            self.environment_id = environment_id

    def set_association_id(self, assoc_id=None):

        if assoc_id is None:
            self.assoc_id = self.make_g2p_id()
        else:
            self.assoc_id = assoc_id

    def add_association_to_graph(self, entity_category=None, phenotype_category=None):
        """
        Overrides  Association by including bnode support

        The reified relationship between a genotype (or any genotype part)
        and a phenotype is decorated with some provenance information.
        This makes the assumption that
        both the genotype and phenotype are classes.

        currently hardcoded to map the annotation to the monarch namespace
        :param g:
        :param entity_category: a biolink category CURIE for self.sub
        :param phenotype_category: a biolink category CURIE for self.obj
        :return:
        """
        # is this kosher?
        Assoc.add_association_to_graph(self)

        # make a blank stage bnode
        if self.start_stage_id or self.end_stage_id is not None:
            stage_process_str = '-'.join((
                str(self.start_stage_id), str(self.end_stage_id)))
            stage_process_id = self.make_id(stage_process_str.replace(':', ''), '_')
            self.model.addIndividualToGraph(
                stage_process_id, None, self.globaltt['developmental_process']
            )
            self.graph.addTriple(
                stage_process_id, self.globaltt['label'], stage_process_str
            )

            self.graph.addTriple(
                stage_process_id, self.globaltt['starts during'], self.start_stage_id
            )

            self.graph.addTriple(
                stage_process_id, self.globaltt['ends during'], self.end_stage_id
            )

            self.stage_process_id = stage_process_id
            self.graph.addTriple(
                self.assoc_id, self.globaltt['has_qualifier'], self.stage_process_id
            )

        if self.environment_id is not None:
            self.graph.addTriple(
                self.assoc_id, self.globaltt['has_qualifier'], self.environment_id
            )

    def make_g2p_id(self):
        """
        Make an association id for phenotypic associations that is defined by:
        source of association +
        (Annot subject) +
        relationship +
        phenotype/disease +
        environment +
        start stage +
        end stage

        :return:

        """

        attributes = [self.environment_id, self.start_stage_id, self.end_stage_id]
        assoc_id = self.make_association_id(
            self.definedby, self.entity_id, self.rel, self.phenotype_id, attributes)

        return assoc_id
