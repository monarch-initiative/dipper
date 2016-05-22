#!/usr/bin/env python3

import unittest
import logging
import csv
import gzip
from tests.test_source import SourceTestCase
from dipper.sources.IMPC import IMPC
from rdflib import URIRef, BNode
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
                           'targeted mutation 1a, Wellcome Trust Sanger Institute',
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

        (p_value, percentage_change, effect_size) = self.test_set_1[23:26]
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
         parameter_name) = self.test_set_1[12:19]
        (statistical_method, resource_name) = self.test_set_1[26:28]

        impc._add_study_provenance(
            impc_map, parameter_map, phenotyping_center, colony,
            project_fullname, pipeline_name, pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name)

        sparql_query = """
                      SELECT ?study ?colony
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

        # This may change if we change our approach for
        # making blank node iris, it might be better
        # to check the length of the output (see test_provenance_mode)
        study = BNode('9328ff6b6455b01254a5548c3cfcc8c4')
        colony = BNode('f3015fe2476ce825a8c6978af6222d87')
        expected_output = [(study, colony)]

        self.assertEqual(list(sparql_output), expected_output)

    def test_assertion_model(self):
        """
        Functional test for _add_study_provenance()
        """
        impc = IMPC()
        impc.load_bindings()
        impc_map = impc._get_impc_mappings()

        impc._add_assertion_provenance(self.assoc_curie,
                                       self.evidence_curie, impc_map)

        sparql_query = """
                      SELECT *
                      WHERE {
                          :MONARCH_test_association OBO:SEPIO_0000015 ?assertion.
                          ?assertion a OBO:SEPIO_0000001 ;
                              OBO:SEPIO_0000018 <http://www.mousephenotype.org/> ;
                              OBO:SEPIO_0000111 _:evidence  .

                          <http://www.mousephenotype.org/> a foaf:organization ;
                              rdfs:label "International Mouse Phenotyping Consortium" .
                      }
                      """

        sparql_output = impc.graph.query(sparql_query)
        # Test that query passes and returns one row
        self.assertEqual(len(list(sparql_output)), 1)

    def test_random_data_set(self):
        """
        Download dataset using fetch(), then take a row of data and
        run through evidence and provenance functions to test the output

        Line of data is hardcoded, but theoretically should work on any line
        """
        line_to_test = 1129
        count = 0
        # init impc (make this a function?)
        impc = IMPC()
        impc.load_bindings()
        impc_map = impc._get_impc_mappings()
        parameter_map = impc._get_parameter_mappings()

        # fetch file
        impc.fetch(True)
        file_path = '/'.join((impc.rawdir, impc.files['all']['file']))
        with gzip.open(file_path, 'rt') as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for row in filereader:
                count += 1
                if count == line_to_test:
                    test_set = row
                    self.test_set_1 = row
                    break

        # Some DRY violation with the above tests
        (phenotyping_center, colony) = row[2:4]
        (project_fullname, pipeline_name, pipeline_stable_id,
         procedure_stable_id, procedure_name, parameter_stable_id,
         parameter_name) = row[12:19]
        (statistical_method, resource_name) = row[26:28]

        (p_value, percentage_change, effect_size) = self.test_set_1[23:26]

        impc._add_evidence(self.assoc_curie, self.eco_id, impc_map, p_value,
                           percentage_change, effect_size, self.study_curie,
                           phenotyping_center)

        impc._add_study_provenance(
            impc_map, parameter_map, phenotyping_center, colony,
            project_fullname, pipeline_name, pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name)

        sparql_query = """
                      SELECT *
                      WHERE {
                          ?assoc OBO:SEPIO_0000007 ?evidenceline .
                          ?evidenceline a OBO:ECO_0000059 ;
                              OBO:SEPIO_0000011 _:study  .

                          ?study a OBO:OBI_0000471 ;
                              OBO:SEPIO_0000114 ?param ;
                              OBO:BFO_0000050 ?project ;
                              OBO:RO_0002233 ?colony ;
                              OBO:SEPIO_0000017 ?agent .
                      }
                      """

        sparql_output = impc.graph.query(sparql_query)
        # Test that query passes and returns one row
        self.assertEqual(len(list(sparql_output)), 1)

    def tearDown(self):
        return


if __name__ == '__main__':
    unittest.main()
