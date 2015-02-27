from dipper.models.Assoc import Assoc
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map


class Chem2DiseaseAssoc(Assoc):
    """ Parse, Fetch, and Convert data from CTD into Triples

    Attributes:
      assoc_id (str): Association Curie (Prefix:ID)
      chem_id (str): Chemical Curie
      phenotype_id (str): Phenotype Curie
      pub_list (str,list): One or more publication curies
      rel (str): Property relating assoc_id and chem_id
      evidence (str): Evidence curie
    """

    def __init__(self, assoc_id, chem_id, phenotype_id, pub_list, relationship, evidence):
        super().__init__()
        self.annot_id = assoc_id
        self.chem_id = chem_id
        self.phenotype_id = phenotype_id
        self.pub_list = pub_list
        self.pub_id = None
        self.evidence = evidence
        self.rel = relationship
        self.cu = CurieUtil(curie_map.get())

        self.setSubject(chem_id)
        self.setObject(phenotype_id)
        return

    def addAssociationNodeToGraph(self, graph):
        """
        The reified relationship between a chemical and a phenotype
        is decorated with some provenance information.

        :param g (rdflib.Graph()): Graph in which to add triples
        :return: graph (rdflib.Graph()): Graph containing reified relationship
        """

        self.addAssociationToGraph(graph)

        return graph