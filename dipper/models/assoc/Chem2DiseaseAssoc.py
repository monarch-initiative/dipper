from dipper.models.assoc.Association import Assoc


class Chem2DiseaseAssoc(Assoc):
    """
    Attributes:
    assoc_id (str): Association Curie (Prefix:ID)
    chem_id (str): Chemical Curie
    phenotype_id (str): Phenotype Curie
    pub_list (str,list): One or more publication curies
    rel (str): Property relating assoc_id and chem_id
    evidence (str): Evidence curie

    """

    def __init__(self, graph, definedby, chem_id, phenotype_id, rel_id=None):
        super().__init__(graph, definedby)
        self.chem_id = chem_id
        self.phenotype_id = phenotype_id

        self.set_subject(chem_id)
        self.set_object(phenotype_id)
        if rel_id is None:
            rel_id = self.gu.object_properties['has_phenotype']
        self.set_relationship(rel_id)

        return

    def set_association_id(self, assoc_id=None):
        """
        This will set the association ID based on the internal parts
        of the association.
        To be used in cases where an external association identifier
        should be used.

        :param assoc_id:

        :return:

        """

        if assoc_id is None:
            self.assoc_id = self.make_c2p_assoc_id()

        return

    def make_c2p_assoc_id(self):

        assoc_id = self.make_association_id(self.definedby, self.sub, self.rel,
                                            self.obj)

        return assoc_id
