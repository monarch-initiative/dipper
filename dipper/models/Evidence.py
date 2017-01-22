import logging
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


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
    evidence_types = {
        'measurement datum': 'IAO:0000109',
        'zscore': 'STATO:0000104',
        'pvalue': 'OBI:0000175',
        'fold_change': 'STATO:0000169',
        'assay': 'OBI:0000070',
        'statistical_hypothesis_test': 'OBI:0000673',
        'effect_size': 'STATO:0000085',
        'percent_change': 'STATO:percent_change',
        'blood test evidence': 'ECO:0001016'
    }

    object_properties = {
        'has_evidence': 'SEPIO:0000006',
        'has_supporting_evidence': 'SEPIO:0000007',
        'has_evidence': 'SEPIO:0000006',
        'has_supporting_data': 'SEPIO:0000084',
        'is_evidence_for': 'SEPIO:0000031',
        'is_refuting_evidence_for': 'SEPIO:0000033',
        'is_supporting_evidence_for': 'SEPIO:0000032',
        'is_evidence_supported_by': 'SEPIO:000010',
        'is_evidence_with_support_from': 'SEPIO:0000059',
        'has_significance': 'STATO:has_significance'
    }

    data_property = {
        'has_value': 'STATO:0000129',
        'has_measurement': 'IAO:0000004'
    }

    def __init__(self, graph):

        self.graph = graph
        self.graph_utils = GraphUtils(curie_map.get())

        return

    def add_supporting_evidence(self, assoc_id, evidence_line, type=None, label=None):
        """
        Add supporting line of evidence node to association id

        :param assoc_id: curie or iri, association id
        :param evidence_line: curie or iri, evidence line
        :return: None
        """
        self.graph_utils.addTriple(self.graph, assoc_id,
                                   self.object_properties['has_supporting_evidence'],
                                   evidence_line)
        if type is not None:
            self.graph_utils.addIndividualToGraph(self.graph, evidence_line,
                                                  label, type)
        return

    def add_evidence(self, assoc_id, evidence_line, type=None, label=None):
        """
        Add line of evidence node to association id

        :param assoc_id: curie or iri, association id
        :param evidence_line: curie or iri, evidence line
        :return: None
        """
        self.graph_utils.addTriple(self.graph, assoc_id,
                                   self.object_properties['has_evidence'],
                                   evidence_line)
        if type is not None:
            self.graph_utils.addIndividualToGraph(self.graph, evidence_line,
                                                  label, type)
        return

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
            self.graph_utils.addTriple(self.graph, evidence_line,
                                       self.object_properties['has_supporting_data'],
                                       measurement)

            self.graph_utils.addTriple(self.graph, measurement,
                                       self.data_property['has_value'],
                                       measurement_dict[measurement], True)
        return
