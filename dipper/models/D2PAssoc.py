__author__ = 'nlw'

from rdflib import Literal

from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils


class D2PAssoc(Assoc):
    """
    A specific association class for defining Disease-to-Phenotype relationships
    This assumes that a graph is created outside of this class, and nodes get added.
    By default, an association will assume the "has_phenotype" relationship, unless
    otherwise specified.
    """

    object_properties = {
        'onset': ':onset',
        'frequency': ':frequencyOfPhenotype'
    }

    def __init__(self, definedby, disease_id, phenotype_id, onset=None, frequency=None, rel=None):
        super().__init__(definedby)
        self.disease_id = disease_id
        self.phenotype_id = phenotype_id
        self.onset = onset
        self.frequency = frequency
        if rel is None:
            rel = self.properties['has_phenotype']

        self.setRelationship(rel)
        self.setSubject(disease_id)
        self.setObject(phenotype_id)

        return

    def addAssociationNodeToGraph(self, g):
        """
        The reified relationship between a disease and a phenotype is decorated with some provenance information.
        This makes the assumption that both the disease and phenotype are classes.

        :param g:
        :return:
        """

        # add the basic association nodes
        self.assoc_id = self.make_d2p_id(self.definedby)
        self.addAssociationToGraph(g)

        if self.frequency is not None and self.frequency != '':
            # FIXME what is the real predicate here?
            self.gu.addTriple(g, self.assoc_id, self.object_properties['frequency'], self.frequency, True)
        if self.onset is not None and self.onset != '':
            # FIXME what is the real predicate here?
            self.gu.addTriple(g, self.assoc_id, self.object_properties['onset'], self.onset)

        return g

    def make_d2p_id(self, definedby):
        """
        Make an association id for phenotypic associations with disease that is defined by:
            source of association + disease + relationship + phenotype
                + onset + frequency
        :param definedby:
        :return:
        """

        attributes = [self.onset, self.frequency]
        assoc_id = self.make_association_id(definedby, self.disease_id, self.rel, self.phenotype_id, attributes)

        return assoc_id

    def loadAllProperties(self, g):

        super().loadAllProperties(g)
        self.gu.loadObjectProperties(g, self.object_properties)

        return