#!/usr/bin/env python3

import unittest
import logging
import csv
import gzip
import json
from tests.test_source import SourceTestCase
from dipper.sources.IMPC import IMPC
from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.TestUtils import TestUtils
from dipper import curie_map
from dipper.graph.RDFGraph import RDFGraph, URIRef


logging.basicConfig()
logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


class IMPCTestCase(SourceTestCase):

    def setUp(self):
        self.source = IMPC('rdf_graph', True)  # Skolem Yes
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


class EvidenceProvenanceTestCase(unittest.TestCase):

    def setUp(self):
        self.test_util = TestUtils()
        self.assoc_curie = 'MONARCH:test_association'
        self.eco_id = 'ECO:0000015'

        self.test_set_1 = (
            'MGI:1920145', 'Setd5', 'WTSI', 'MEFW', 'male', 'heterozygote',
            'MGI:4432631', 'Setd5<tm1a(EUCOMM)Wtsi>',
            'targeted mutation 1a, Wellcome Trust Sanger Institute',
            'MGI:2159965', 'C57BL/6N', 'MGP',
            'Wellcome Trust Sanger Institute Mouse Genetics Project',
            'MGP Select Pipeline', 'MGP_001', 'MGP_XRY_001', 'X-ray',
            'IMPC_XRY_008_001', 'Number of ribs right', 'MP:0005390',
            'skeleton phenotype', 'MP:0000480', 'increased rib number',
            '1.637023E-010', '', '8.885439E-007',
            'Wilcoxon rank sum test with continuity correction',
            'International Mouse Phenotyping Consortium')

        # Generate test curies, these are otherwise generated
        # within _add_evidence() and _add_study_provenance()
        # these blank nodes are hardcoded as NOT Skolemized  ...
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
        impc = IMPC('rdf_graph', True)
        impc.graph = RDFGraph(True)  # Reset graph
        # Test graph is empty
        self.assertTrue(len(list(impc.graph)) == 0)

        impc_map = impc.open_and_parse_yaml(impc.map_files['impc_map'])

        (p_value, percentage_change, effect_size) = self.test_set_1[23:26]

        impc._add_evidence(self.assoc_curie, self.eco_id, impc_map, p_value,
                           percentage_change, effect_size, self.study_curie)

        triples = """
    :MONARCH_test_association SEPIO:0000007 <https://monarchinitiative.org/.well-known/genid/b097a98087df7a99> .

    <https://monarchinitiative.org/.well-known/genid/b097a98087df7a99> a ECO:0000015 ;
        SEPIO:0000084 <https://monarchinitiative.org/.well-known/genid/b89ee584330837c9>,
            <https://monarchinitiative.org/.well-known/genid/bc0eeccdea27a1d8> ;
        SEPIO:0000085 <https://monarchinitiative.org/.well-known/genid/study> .

    <https://monarchinitiative.org/.well-known/genid/bc0eeccdea27a1d8> a OBI:0000175 ;
        RO:0002353 <https://monarchinitiative.org/.well-known/genid/study> ;
        STATO:0000129 1.637023e-10 .

    <https://monarchinitiative.org/.well-known/genid/b89ee584330837c9> a STATO:0000085 ;
        RO:0002353 <https://monarchinitiative.org/.well-known/genid/study> ;
        STATO:0000129 "8.885439E-007" .
        """

        self.assertTrue(self.test_util.test_graph_equality(
            triples, impc.graph))

    def test_provenance_model(self):
        """
        Functional test for _add_study_provenance()
        """
        impc = IMPC('rdf_graph', True)
        impc.graph = RDFGraph(True)
        self.assertTrue(len(list(impc.graph)) == 0)

        impc_map = impc.open_and_parse_yaml(impc.map_files['impc_map'])
        impress_map = json.loads(
            impc.fetch_from_url(
                impc.map_files['impress_map']).read().decode('utf-8'))

        (phenotyping_center, colony) = self.test_set_1[2:4]
        (project_fullname, pipeline_name, pipeline_stable_id,
         procedure_stable_id, procedure_name, parameter_stable_id,
         parameter_name) = self.test_set_1[12:19]
        (statistical_method, resource_name) = self.test_set_1[26:28]

        impc._add_study_provenance(
            impc_map, impress_map, phenotyping_center, colony,
            project_fullname, pipeline_name, pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name)

        triples = """
    <https://monarchinitiative.org/.well-known/genid/bbdd05a8ca155dda> a OBI:0000471 ;
      BFO:0000051 OBO:STATO_0000076,
          <https://www.mousephenotype.org/impress/protocol/175/15> ;
      BFO:0000050  IMPRESS-procedure:15 ,
          <http://www.sanger.ac.uk/science/data/mouse-genomes-project> ;
      SEPIO:0000114 <https://www.mousephenotype.org/impress/parameterontologies/1867/91> ;
      SEPIO:0000017 <http://www.sanger.ac.uk/>  .

    <https://monarchinitiative.org/.well-known/genid/bc0b26361b8687b5> a owl:NamedIndividual ;
        rdfs:label "MEFW" .

    <http://www.sanger.ac.uk/> a foaf:organization ;
        rdfs:label "WTSI" .

    <http://www.sanger.ac.uk/science/data/mouse-genomes-project> a VIVO:Project ;
        rdfs:label "Wellcome Trust Sanger Institute Mouse Genetics Project" .

    <https://www.mousephenotype.org/impress/parameterontologies/1867/91> a owl:NamedIndividual ;
        rdfs:label "Number of ribs right (X-ray)" .

    IMPRESS-procedure:15 a owl:NamedIndividual ;
        rdfs:label "MGP Select Pipeline" .

    <https://www.mousephenotype.org/impress/protocol/175/15> a owl:NamedIndividual ;
        rdfs:label "X-ray" .
"""
        # dbg
        logger.debug(
            "Reference graph: %s", impc.graph.serialize(format="turtle").decode("utf-8")
        )
        # bitrot test
        #self.assertTrue(
        #    self.test_util.test_graph_equality(triples, impc.graph))

    def test_assertion_model(self):
        """
        Functional test for _add_study_provenance()
        """

        impc = IMPC('rdf_graph', True)
        impc.graph = RDFGraph(True)
        self.assertTrue(len(list(impc.graph)) == 0)

        impc_map = impc.open_and_parse_yaml(impc.map_files['impc_map'])

        impc._add_assertion_provenance(
            self.assoc_curie, self.evidence_curie, impc_map)

        triples = """
    MONARCH:test_association SEPIO:0000015 <https://monarchinitiative.org/.well-known/genid/bcb2c00a5c2f9c43> .
    <https://monarchinitiative.org/.well-known/genid/bcb2c00a5c2f9c43> a SEPIO:0000001 ;
        SEPIO:0000018 <http://www.mousephenotype.org/> ;
        SEPIO:0000111 <https://monarchinitiative.org/.well-known/genid/evidence>  .

    <http://www.mousephenotype.org/> a foaf:organization ;
        rdfs:label "International Mouse Phenotyping Consortium" .

        """
        # dbg
        logger.debug(
            "Reference graph: %s", impc.graph.serialize(format="turtle").decode("utf-8")
        )

        self.assertTrue(self.test_util.test_graph_equality(
            triples, impc.graph))

    def test_random_data_set(self):
        """
        Download dataset using fetch(), then take a row of data and
        run through evidence and provenance functions to test the output

        Line of data is hardcoded, but theoretically should work on any line
        """
        line_to_test = 1129
        count = 0
        impc = IMPC('rdf_graph', False)   # Not Skolem
        impress_map = json.loads(
            impc.fetch_from_url(
                impc.map_files['impress_map']).read().decode('utf-8'))
        impc_map = impc.open_and_parse_yaml(impc.map_files['impc_map'])

        # fetch file
        impc.fetch(True)
        file_path = '/'.join((impc.rawdir, impc.files['all']['file']))
        with gzip.open(file_path, 'rt') as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for row in filereader:
                count += 1
                if count == line_to_test:
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
                           percentage_change, effect_size, self.study_curie)

        impc._add_study_provenance(
            impc_map, impress_map, phenotyping_center, colony,
            project_fullname, pipeline_name, pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name)

        # Note that this doesn't test much since we're dealing with
        # multiple part_of  and has_part links to individuals
        # which results in ambiguity = hard to test
        sparql_query = """
SELECT *
WHERE {
    ?assoc SEPIO:0000007 ?evidenceline .
    ?evidenceline a ECO:0000015 ;
        SEPIO:0000085 _:study .

    ?study a OBI:0000471 ;
        SEPIO:0000114 ?param ;
        SEPIO:0000017 ?agent .
}
"""
        sparql_output = impc.graph.query(sparql_query)
        # Test that query passes and returns one row
        self.assertEqual(len(list(sparql_output)), 1)

    def tearDown(self):
        return


if __name__ == '__main__':
    unittest.main()
