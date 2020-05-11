import logging
from dipper.models.Model import Model
from dipper.graph.Graph import Graph

LOG = logging.getLogger(__name__)
# note currently no log issued


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

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".format(graph))
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        return

    def add_date_created(self, prov_type, date):
        self.graph.addTriple(
            object_is_literal=True, subject_id=prov_type,
            predicate_id=self.globaltt['created_on'], obj=date)
        return

    def add_study_parts(self, study, study_parts):
        for part in study_parts:
            self.graph.addTriple(
                study, self.globaltt['has_part'], part)
        return

    def add_study_to_measurements(self, study, measurements):
        for measurement in measurements:
            self.graph.addTriple(
                measurement, self.globaltt['output_of'], study)
        return

    def add_study_measure(self, study, measure, object_is_literal=None):
        self.graph.addTriple(
            study, self.globaltt['measures_parameter'], measure, object_is_literal)
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
        self.model.addIndividualToGraph(assertion, None, self.globaltt['assertion'])

        self.add_agent_to_graph(agent, agent_label, self.globaltt['organization'])

        self.graph.addTriple(
            assertion, self.globaltt['created_by'], agent)

        if date is not None:
            self.graph.addTriple(
                self.graph, assertion, self.globaltt['Date Created'], date)

        return

    def add_agent_to_graph(
            self, agent_id, agent_label, agent_type=None, agent_description=None):

        if agent_type is None:
            agent_type = self.globaltt['organization']
        self.model.addIndividualToGraph(
            agent_id, agent_label, agent_type, agent_description)

        return

    def add_assay_to_graph(
            self, assay_id, assay_label, assay_type=None, assay_description=None):
        if assay_type is None:
            assay_type = self.globaltt['assay']
        self.model.addIndividualToGraph(
            assay_id, assay_label, assay_type, assay_description)

        return
