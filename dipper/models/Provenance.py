__author__ = 'nlw'

from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
import logging
import re
from datetime import datetime
import hashlib



logger = logging.getLogger(__name__)


class Provenance:
    """
    To model provenance as the basis for an association.
    This encompases:
        * measurements taken from the lab, and their significance.  these can be derived from papers or other agents.
        * papers
    >1 measurement may result from an assay, each of which may have it's own significance


    Status:  IN PROGRESS  (as observed from the incomplete identifiers

    TODO: add in

    """
    prov_types = {
        'measurement datum': 'IAO:0000109',
        'zscore': 'STATO:0000104',
        'pvalue': 'OBI:0000175',
        'assay': 'OBI:0000070',
        'agent': 'IAO:xxxxxxx',
        'organization': 'foaf:organization',
        'person': 'foaf:person',
        'statistical_hypothesis_test': 'OBI:0000673'
    }

    rel_properties = {
        'has_significance': 'STATO:has_significance',  # FIXME using has_specified_output
        'has_agent': 'RO:has_agent'  # FIXME

    }

    data_property = {
        'has_value': 'STATO:0000129',  # FIXME should this be in RO?
        'has_measurement': 'IAO:0000004'
    }

    def __init__(self, prov_type=None):

        if prov_type is None:
            self.prov_type = 'OBAN:provenance'
        self.prov_id = None
        self.measurement_datums = {}
        self.agent = None
        self.reference = None  # TODO this will be papers in the future
        self.gu = GraphUtils(curie_map.get())

        return

    def set_prov_id(self, prov_id):
        self.prov_id = prov_id
        return

    def get_prov_id(self):

        return self.prov_id

    def add_measurement_data(self, assay_id, measurement_unit, significance=None):
        """
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
            # @N @M This needs to be a list, not a set, because multiple values may be possible for a given assay?
            self.measurement_datums[assay_id] = list()

        ms = self.get_score(self.prov_types['measurement datum'], measurement_unit)
        if significance is not None:
            ms['significance'] = significance

        self.measurement_datums[assay_id].append(ms)

        # as a convenience, we return the measurement data
        return ms

    def get_zscore(self, zscore_value):
        """
        This will get you a significance scoring object
        :param zscore_value:
        :return:
        """
        s = self.get_score(self.prov_types['zscore'], zscore_value)
        return s

    def get_score_id(self, score_type, score_value, score_unit=None):
        s = '-'.join((re.sub(':','',score_type), str(score_value), str(score_unit)))

        """
        a method to create unique identifiers based on very long strings
        currently implemented with md5
        :param long_string:
        :return:
        @N, not sure what this should actually be?
        """
        # FIXME for now, this will do md5.  probably not the best long-term solution
        # note others available: md5(), sha1(), sha224(), sha256(), sha384(), and sha512()

        byte_string = s.encode("utf-8")

        return ':'.join(('MONARCH', hashlib.md5(byte_string).hexdigest()))

    def get_score(self, score_type, score_value, score_unit=None):

        # TODO if score_type or score_value is None, then throw error

        score = {
            'value':score_value,
            'type': score_type
        }

        if score_unit is not None:
            score['unit'] = score_unit

        return score

    def add_provenance_to_graph(self, graph):

        # we make the assumption that the agent has already been added to the graph, ok?
        self.add_measurement_data_to_graph(graph, self.measurement_datums)

        return

    def add_agent_to_graph(self, graph, agent_id, agent_label, agent_type=None, agent_description=None):

        if agent_type is None:
            agent_type = self.prov_types['agent']
        self.gu.addIndividualToGraph(graph, agent_id, agent_label, agent_type, agent_description)
        self.agent = agent_id

        return

    def add_measurement_data_to_graph(self, graph, measurement_data):

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
        #     # we assume that the assay has already been added as an Individual with properties elsewhere
        #
        #     for i in s:
        #
        #         logger.debug("\tTYPE: "+str(i.get('type'))+"\tVALUE: "+str(i.get('value'))+"\tUNIT: "+str(i.get('unit'))+" TYPE: None\tDESCRIPTION:"+str(i.get('description')))
        #         mid = self.get_score_id(i.get('type'), i.get('value'), i.get('unit'))
        #
        #         self.gu.addIndividualToGraph(graph, mid, None, i.get('description'))
        #         self.gu.addTriple(graph, mid, self.data_property['has_value'], i.get('value'))
        #         sig = i.get('significance')
        #         if sig is not None:
        #             sid = self.get_score_id(i.get('type'), i.get('value'), i.get('unit'))
        #             self.gu.addIndividualToGraph(graph, sid, None, s.get('type'))
        #             self.gu.addTriple(graph, mid, self.rel_properties['has_significance'], sid)
        #             self.gu.addTriple(graph, sid, self.data_property['has_value'], sig.get('value'))
        #         # add the measurement as part of this provenance
        #         self.gu.addTriple(graph, self.prov_id, self.gu.object_properties['has_part'], mid)  # TODO should this be "has_measurement"?
        #         self.gu.addTriple(graph, self.prov_id, self.gu.object_properties['has_agent'], self.agent)

        return

    def add_assay_to_graph(self, graph, assay_id, assay_label, assay_type=None, assay_description=None):
        if assay_type is None:
            assay_type = self.prov_types['assay']
        self.gu.addIndividualToGraph(graph, assay_id, assay_label, assay_type, assay_description)

        return