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

    def __init__(self, assoc_id, entity_id, phenotype_id, onset, frequency, pub, evidence_code):
        super().__init__()
        self.annot_id = assoc_id
        self.entity_id = entity_id
        self.phenotype_id = phenotype_id
        self.onset = onset
        self.frequency = frequency
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.properties['has_phenotype']
        self.cu = CurieUtil(curie_map.get())

        self.pub_list = None

        self.setSubject(entity_id)
        self.setObject(phenotype_id)

        return

    def set_relationship(self, rel):
        self.rel = rel
        return

    def addAssociationNodeToGraph(self, g):
        """
        The reified relationship between a disease and a phenotype is decorated with some provenance information.
        This makes the assumption that both the disease and phenotype are classes.

        :param g:
        :return:
        """

        # add the basic association nodes
        self.addAssociationToGraph(g)

        # add the specific attributes for this association type
        gu = GraphUtils(curie_map.get())
        node = gu.getNode(self.annot_id)

        if self.frequency is not None and self.frequency != '':
            # FIXME what is the real predicate here?
            g.add((node, self.BASE['frequencyOfPhenotype'], Literal(self.frequency)))
        if self.onset is not None and self.onset != '':
            # FIXME what is the real predicate here?
            g.add((node, self.BASE['onset'], gu.getNode(self.onset)))

        return g