#!/usr/bin/env python3

import unittest
import logging
# import os
# from rdflib import Graph
# from tests import test_general, test_source
from tests.test_source import SourceTestCase
from dipper.sources.IMPC import IMPC
from rdflib.namespace import URIRef
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class IMPCTestCase(SourceTestCase):

    def setUp(self):
        self.source = IMPC()
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # @unittest.skip('test not yet defined')
    # def test_hpotest(self):
    #    logger.info("An IMPC-specific test")
    #
    #    return


class EvidenceProvenanceTestCase(unittest.TestCase):

    def setUp(self):
        self.assoc_curie = 'MONARCH:test_association'
        self.eco_id = 'ECO:0000059'

        self.test_set_1 = ('MGI:1920145', 'Setd5', 'WTSI', 'MEFW', 'male',
                           'heterozygote', 'MGI:4432631', 'Setd5<tm1a(EUCOMM)Wtsi>',
                           'targeted mutation 1a', 'Wellcome Trust Sanger Institute',
                           'MGI:2159965', 'C57BL/6N', 'MGP',
                           'Wellcome Trust Sanger Institute Mouse Genetics Project',
                           'MGP Select Pipeline', 'MGP_001', 'MGP_XRY_001', 'X-ray',
                           'IMPC_XRY_008_001', 'Number of ribs right', 'MP:0005390',
                           'skeleton phenotype', 'MP:0000480', 'increased rib number',
                           '1.637023E-010', '', '8.885439E-007',
                           'Wilcoxon rank sum test with continuity correction', 'IMPC')

        # Generate test curies, these are otherwise generated
        # within _add_evidence() and _add_study_provenance()
        self.study_curie = "_:study"
        self.evidence_curie = "_:evidence"

        # IRIs for testing sparql output
        curie_dict = curie_map.get()
        curie_util = CurieUtil(curie_dict)
        self.assoc_iri = URIRef(curie_util.get_uri(self.assoc_curie))

        return

    def test_evidence_model(self):
        """
        Functional test for _add_evidence()
        """
        impc = IMPC()
        impc.load_bindings()
        impc_map = impc._get_impc_mappings()

        (p_value, percentage_change, effect_size) = self.test_set_1[24:27]
        phenotyping_center = self.test_set_1[2]

        impc._add_evidence(self.assoc_curie, self.eco_id, impc_map, p_value,
                           percentage_change, effect_size, self.study_curie,
                           phenotyping_center)

        sparql_query = """
                      SELECT ?assoc
                      WHERE {
                            ?assoc OBO:SEPIO_0000007 ?evidenceline .
                            ?evidenceline a OBO:ECO_0000059 ;
                                OBO:BFO_0000051 ?measure1 ;
                                OBO:BFO_0000051 ?measure2 ;
                                OBO:SEPIO_0000011 _:study  .

                            ?measure1 a OBO:OBI_0000175 ;
                                OBO:RO_0002353 _:study ;
                                OBO:STATO_0000129 1.637023e-10 .

                            ?measure2 a OBO:STATO_0000085 ;
                                OBO:RO_0002353 _:study ;
                                OBO:STATO_0000129 "8.885439E-007" .

                      }
                      """
        sparql_output = impc.graph.query(sparql_query)
        expected_results = [(self.assoc_iri,)]

        self.assertEqual(list(sparql_output), expected_results)

    def test_provenance_model(self):
        """
        Functional test for _add_study_provenance()
        """
        impc = IMPC()
        impc.load_bindings()
        impc_map = impc._get_impc_mappings()
        parameter_map = impc._get_parameter_mappings()

        (phenotyping_center, colony) = self.test_set_1[2:4]
        (project_fullname, pipeline_name, pipeline_stable_id,
         procedure_stable_id, procedure_name, parameter_stable_id,
         parameter_name) = self.test_set_1[13:20]
        (statistical_method, resource_name) = self.test_set_1[27:29]

        impc._add_study_provenance(
            impc_map, parameter_map, phenotyping_center, colony,
            project_fullname, pipeline_name, pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name)

        sparql_query = """
                      SELECT *
                      WHERE {
                          <https://www.mousephenotype.org/impress/procedures/15> a owl:NamedIndividual ;
                              rdfs:label "MGP Select Pipeline" .

                          <https://www.mousephenotype.org/impress/protocol/175/15> a owl:NamedIndividual ;
                              rdfs:label "X-ray" .

                          <http://www.sanger.ac.uk/> a foaf:organization ;
                              rdfs:label "WTSI" .

                          <http://www.sanger.ac.uk/science/data/mouse-genomes-project> a VIVO:Project ;
                              rdfs:label "Wellcome Trust Sanger Institute Mouse Genetics Project" .

                          <https://www.mousephenotype.org/impress/parameterontologies/1867/175> a owl:NamedIndividual ;
                              rdfs:label "Number of ribs right" .

                          ?study a OBO:OBI_0000471 ;
                              OBO:BFO_0000051 OBO:STATO_0000076 ;
                              OBO:BFO_0000051 <https://www.mousephenotype.org/impress/procedures/15> ;
                              OBO:BFO_0000051 <https://www.mousephenotype.org/impress/protocol/175/15> ;
                              OBO:SEPIO_0000114 <https://www.mousephenotype.org/impress/parameterontologies/1867/175> ;
                              OBO:BFO_0000050 <http://www.sanger.ac.uk/science/data/mouse-genomes-project> ;
                              OBO:RO_0002233 ?colony ;
                              OBO:SEPIO_0000017 <http://www.sanger.ac.uk/> .

                          ?colony a owl:NamedIndividual ;
                              rdfs:label "MEFW" .
                      }
                      """

        sparql_output = impc.graph.query(sparql_query)
        # Should output a single row for ?study and ?colony
        # print(sparql_output)
        # >> [[rdflib.term.BNode('f3015fe2476ce825a8c6978af6222d87'),
        #      rdflib.term.BNode('9328ff6b6455b01254a5548c3cfcc8c4')]]

        self.assertEqual(len(list(sparql_output)[0]), 2)

    def tearDown(self):
        return


if __name__ == '__main__':
    unittest.main()
