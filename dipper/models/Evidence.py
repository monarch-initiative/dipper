import logging
import re
from datetime import datetime
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

__author__ = 'nlw'

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

    def _add_measurement_data(self, assay_id, measurement_unit,
                             significance=None):
        """

        This is a legacy function to handle MPD and was never
        implemented.  Leaving here for now but will fix/delete
        as needed

        Adds an object to measurement_datums like:
        assay_id : {
            # MB says each assay should be an instance of OBI Assay
            type: OBI:0000070
            label: appmeth
            unit: unit id or literal?
            significance: {
                  type:  significance class id
                  value: ####
                  unit:  unit id (optional)
              } (optional)
        }

        :param assay_id:
        :param measurement_value:
        :param measurement_unit:
        :param significance:
        :return:
        """

        if assay_id not in self.measurement_datums:
            # @N @M This needs to be a list, not a set,
            # because multiple values may be possible for a given assay?
            self.measurement_datums[assay_id] = list()

        ms = self.get_score(self.prov_types['measurement datum'],
                            measurement_unit)
        if significance is not None:
            ms['significance'] = significance

        self.measurement_datums[assay_id].append(ms)

        # as a convenience, we return the measurement data
        return ms


    def _add_measurement_data_to_graph(self, graph, measurement_data):
        """
        This is a legacy function to handle MPD and was never
        implemented.  Leaving here for now but will fix/delete
        as needed
        :param graph:
        :param measurement_data:
        :return:
        """

        # assay_id : {
        #     type: "measurement datum" class
        #     value: #####,
        #     unit: unit id or literal?
        #     significance: {
        #           type:  significance class id
        #           value: ####
        #           unit:  unit id (optional)
        #       } (optional)
        # }

        # if no prov_id, then make it a random blank node
        if self.prov_id is None:
            t = datetime.now()
            t_string = t.strftime("%Y-%m-%d-%H-%M")
            self.prov_id = '_'+str(self.agent)+t_string

        # # TODO deal with units
        # for m in measurement_data:
        #     s = measurement_data[m]
        #     # we assume that the assay has already been added
        #     # as an Individual
        #       with properties elsewhere
        #
        #     for i in s:
        #
        #         logger.debug("\tTYPE: "+str(i.get('type'))+"\tVALUE: "+
        #                      str(i.get('value'))+"\tUNIT: "+
        #                      str(i.get('unit'))+" TYPE: None\tDESCRIPTION:"+
        #                      str(i.get('description')))
        #         mid = self.get_score_id(i.get('type'), i.get('value'),
        #                                 i.get('unit'))
        #
        #         self.graph_utils.addIndividualToGraph(graph, mid, None,
        #                                       i.get('description'))
        #         self.graph_utils.addTriple(
        #           graph, mid, self.data_property['has_value'],
        #                           i.get('value'))
        #         sig = i.get('significance')
        #         if sig is not None:
        #             sid = self.get_score_id(i.get('type'), i.get('value'),
        #                                     i.get('unit'))
        #             self.graph_utils.addIndividualToGraph(graph, sid, None,
        #                                          s.get('type'))
        #             self.graph_utils.addTriple(graph, mid,
        #                               self.rel_properties['has_significance'],
        #                               sid)
        #             self.graph_utils.addTriple(graph, sid,
        #                               self.data_property['has_value'],
        #                               sig.get('value'))
        #         # add the measurement as part of this provenance
        #         self.graph_utils.addTriple(graph, self.prov_id,
        #                           # TODO should this be "has_measurement"?
        #                           self.graph_utils.object_properties['has_part'], mid)
        #         self.graph_utils.addTriple(graph, self.prov_id,
        #                           self.graph_utils.object_properties['has_agent'],
        #                           self.agent)

        return
