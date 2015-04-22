__author__ = 'nlw'

from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map


class DispositionAssoc(Assoc):
    """
    A specific Association model for Heritability annotations.  These are to be used between diseases and a
    heritability disposition.
    """

    def __init__(self, annot_id, entity_id, heritability_id, pub, evidence_code):
        super().__init__()
        self.annot_id = annot_id
        self.entity_id = entity_id
        self.heritability_id = heritability_id
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.object_properties['has_disposition']  # default to has_disposition
        self.cu = CurieUtil(curie_map.get())
        self.gu = GraphUtils(curie_map.get())
        self.pub_list = None

        self.setSubject(entity_id)
        self.setObject(heritability_id)

        return

    def set_relationship(self, rel):
        self.rel = rel
        return

    def addAssociationNodeToGraph(self, g):
        """
        The reified relationship between a disease and the heritability is decorated with some provenance information.
        This makes the assumption that both the disease and heritability are classes.

        :param g:
        :return:
        """

        # add the basic association nodes
        self.addAssociationToGraph(g)

        return g