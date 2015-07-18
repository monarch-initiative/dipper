__author__ = 'nicole'

from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map


class G2PAssoc(Assoc):
    """
    A specific association class for defining Genotype-to-Phenotype relationships
    This assumes that a graph is created outside of this class, and nodes get added.
    By default, an association will assume the "has_phenotype" relationship, unless
    otherwise specified.
    Note that genotypes are expected to be created and defined outside of this association,
    most likely by calling methods in the Genotype() class.
    """

    def __init__(self, definedby, entity_id, phenotype_id, rel=None):
        super().__init__(definedby)
        self.entity_id = entity_id
        self.phenotype_id = phenotype_id
        self.pub_id = None
        self.evidence = None
        if rel is None:
            self.rel = self.properties['has_phenotype']  # default to has_phenotype
        else:
            self.rel = rel
        self.start_stage_id = None
        self.end_stage_id = None
        self.environment_id = None

        self.setSubject(entity_id)
        self.setObject(phenotype_id)
        return

    def set_relationship(self, rel):
        self.rel = rel
        return

    def set_stage(self, start_stage_id, end_stage_id):
        if start_stage_id is not None and start_stage_id.strip() != '':
            self.start_stage_id = start_stage_id
        if end_stage_id is not None and end_stage_id.strip() != '':
            self.end_stage_id = end_stage_id
        return

    def set_environment(self, environment_id):
        if environment_id is not None and environment_id.strip() != '':
            self.environment_id = environment_id

        return

    def addAssociationToGraph(self, g):
        """
        The reified relationship between a genotype (or any genotype part) and a phenotype
        is decorated with some provenance information.
        This makes the assumption that both the genotype and phenotype are classes.

        currently hardcoded to map the annotation to the monarch namespace
        :param g:
        :return:
        """
        self.annot_id = self.make_g2p_id(self.definedby)

        self.addAssociationToGraph(g)

        if self.start_stage_id is not None:
            self.gu.addTriple(g, self.annot_id,
                              self.gu.object_properties['has_begin_stage_qualifier'],
                              self.start_stage_id)
        if self.end_stage_id is not None:
            self.gu.addTriple(g, self.annot_id,
                              self.gu.object_properties['has_end_stage_qualifier'],
                              self.end_stage_id)

        if self.environment_id is not None:
            self.gu.addTriple(g, self.annot_id,
                              self.gu.object_properties['has_environment_qualifier'],
                              self.environment_id)
        return g


    def make_g2p_id(self, definedby):
        """
        Make an association id for phenotypic associations that is defined by:
            source of association + (Annot subject) + relationship + phenotype/disease
                + environment + start stage + end stage
        :param definedby:
        :return:
        """

        attributes = [self.environment_id, self.start_stage_id, self.end_stage_id]
        assoc_id = self.make_association_id(definedby, self.entity_id, self.rel, self.phenotype_id, attributes)

        return assoc_id