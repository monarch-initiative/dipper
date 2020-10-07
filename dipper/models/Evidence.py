import logging
from dipper.models.Model import Model
from dipper.graph.Graph import Graph

LOG = logging.getLogger(__name__)
# note: currently no log issued


class Evidence:
    """
    To model evidence as the basis for an association.
    This encompasses:
        * measurements taken from the lab, and their significance.
            these can be derived from papers or other agents.
        * papers
    >1 measurement may result from an assay,
        each of which may have it's own significance

    """

    def __init__(self, graph, association):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".format(graph))
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        self.association = association

    def add_supporting_evidence(self, evidence_line, evidence_type=None, label=None):
        """
        Add supporting line of evidence node to association id

        :param evidence_line: curie or iri, evidence line
        :param evidence_type: curie or iri, evidence type if available
        :return: None
        """
        self.graph.addTriple(
            self.association,
            self.globaltt['has_supporting_evidence_line'],
            evidence_line)
        if evidence_type is not None:
            self.model.addIndividualToGraph(evidence_line, label, evidence_type)

    def add_evidence(self, evidence_line, evidence_type=None, label=None):
        """
        Add line of evidence node to association id

        :param evidence_line: curie or iri, evidence line
        :param evidence_type: curie or iri, evidence type if available
        :return: None
        """
        self.graph.addTriple(
            self.association, self.globaltt['has_evidence_line'], evidence_line)
        if evidence_type is not None:
            self.model.addIndividualToGraph(evidence_line, label, evidence_type)

    def add_data_individual(
            self, data_curie, label=None, ind_type=None, data_curie_category=None
    ):
        """
        Add data individual
        :param data_curie: str either curie formatted or long string,
                           long strings will be converted to bnodes
        :param data_curie_category: a biolink category CURIE for data_curie
        :param type: str curie
        :param label: str
        :return: None
        """
        part_length = len(data_curie.split(':'))
        if part_length == 0:
            curie = "_:{}".format(data_curie)
        elif part_length > 2:
            raise ValueError("Misformatted curie {}".format(data_curie))
        else:
            curie = data_curie

        self.model.addIndividualToGraph(
            curie, label, ind_type, ind_category=data_curie_category)

    def add_supporting_data(self, evidence_line, measurement_dict):
        """
        Add supporting data
        :param evidence_line:
        :param data_object: dict, where keys are curies or iris
        and values are measurement values for example:
            {
              "_:1234" : "1.53E07"
              "_:4567": "20.25"
            }
        Note: assumes measurements are RDF:Type 'ed elsewhere
        :return: None
        """
        for measurement in measurement_dict:
            self.graph.addTriple(
                evidence_line, self.globaltt['has_evidence_item'], measurement)

            if measurement_dict[measurement] != '':
                self.graph.addTriple(
                    measurement, self.globaltt['has_value'],  # has measurement value ??
                    measurement_dict[measurement], object_is_literal=True)

    def add_supporting_publication(
            self, evidence_line, publication, label=None, pub_type=None
    ):
        """
        <evidence> <has_supporting_reference> <source>
        <source> <rdf:type> <type>
        <source> <rdfs:label> "label"
        :param evidence_line: str curie
        :param publication: str curie
        :param label: optional, str type as curie
        :param type: optional, str type as curie
        :return:
        """
        self.graph.addTriple(
            evidence_line,
            self.globaltt['has_supporting_reference'],
            publication)

        self.model.addIndividualToGraph(publication, label, pub_type)

    def add_source(
            self,
            evidence_line,
            source,
            label=None,
            src_type=None,
            source_category=None
    ):
        """
        Applies the triples:
        <evidence> <dc:source> <source>
        <source> <rdf:type> <type>
        <source> <rdfs:label> "label"

        TODO this should belong in a higher level class
        :param evidence_line: str curie
        :param source: str source as curie
        :param label: optional, str type as curie
        :param type: optional, str type as curie
        :return: None
        """
        self.graph.addTriple(
            evidence_line,
            self.globaltt['Source'],
            source,
            subject_category=None,
            object_category=source_category
        )

        self.model.addIndividualToGraph(
            source, label, src_type, ind_category=source_category
        )
