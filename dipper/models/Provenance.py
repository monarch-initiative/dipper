import logging
from dipper.models.Model import Model
from dipper.graph.Graph import Graph

logger = logging.getLogger(__name__)


class Provenance:
    """
    To model provenance as the basis for an association.
    This encompasses:
        * Process history leading to a claim being made,
          including processes through which evidence is evaluated
        * Processes through which information used as evidence is created.

    Provenance metadata includes accounts of who conducted these processes,
     what entities participated in them, and when/where they occurred.

    """
    provenance_types = {
        'assertion': 'SEPIO:0000001',
        'assay': 'OBI:0000070',
        'organization': 'foaf:organization',
        'person': 'foaf:person',
        'statistical_hypothesis_test': 'OBI:0000673',
        'mixed_model': 'STATO:0000189',
        'project': 'VIVO:Project',
        'study': 'OBI:0000471',
        'variant_classification_guideline': 'SEPIO:0000037',
        'assertion_process': 'SEPIO:0000003',
        'xref': 'OIO:hasdbxref'
    }

    object_properties = {
        'has_provenance': 'SEPIO:0000011',
        'has_participant': 'RO:0000057',
        'has_input': 'RO:0002233',
        'has_agent': 'SEPIO:0000017',
        'created_by': 'SEPIO:0000018',
        'date_created': 'SEPIO:0000021',
        'is_assertion_supported_by': 'SEPIO:0000111',
        'is_asserted_in': 'SEPIO:0000015',
        'output_of': 'RO:0002353',
        'specified_by': 'SEPIO:0000041',
        'created_at_location': 'SEPIO:0000019',
        'created_with_resource': 'SEPIO:0000022',
        'measures': 'SEPIO:0000114',
        'has_supporting_study': 'SEPIO:0000085',
        'asserted_by': 'SEPIO:0000130',
        'created_on': 'pav:createdOn'
    }

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)

        return

    def add_date_created(self, prov_type, date):
        self.graph.addTriple(object_is_literal=True,
                             subject_id=prov_type,
                             predicate_id=Provenance.object_properties['created_on'],
                             obj=date
                             )
        return

    def add_study_parts(self, study, study_parts):
        for part in study_parts:
            self.graph.addTriple(
                study, self.model.object_properties['has_part'], part)
        return

    def add_study_to_measurements(self, study, measurements):
        for measurement in measurements:
            self.graph.addTriple(
                measurement, self.object_properties['output_of'], study)
        return

    def add_study_measure(self, study, measure):
        self.graph.addTriple(
            study, self.object_properties['measures'], measure)
        return

    def add_assertion(self, assertion, agent, agent_label, date=None):
        """
        Add assertion to graph
        :param assertion:
        :param agent:
        :param evidence_line:
        :param date:
        :return: None
        """
        self.model.addIndividualToGraph(
            assertion, None, self.provenance_types['assertion'])

        self.add_agent_to_graph(agent, agent_label,
                                self.provenance_types['organization'])

        self.graph.addTriple(
            assertion, self.object_properties['created_by'], agent)

        if date is not None:
            self.graph.addTriple(
                self.graph, assertion,
                self.object_properties['date_created'], date)

        return

    def add_agent_to_graph(self, agent_id, agent_label, agent_type=None,
                           agent_description=None):

        if agent_type is None:
            agent_type = self.provenance_types['organization']
        self.model.addIndividualToGraph(
            agent_id, agent_label, agent_type, agent_description)

        return

    def add_assay_to_graph(self, assay_id, assay_label, assay_type=None,
                           assay_description=None):
        if assay_type is None:
            assay_type = self.provenance_types['assay']
        self.model.addIndividualToGraph(
            assay_id, assay_label, assay_type, assay_description)

        return
