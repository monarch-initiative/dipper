#!/usr/bin/env python3

import unittest
import logging
import csv
import gzip
from dipper.sources.IMPC import IMPC
from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.TestUtils import TestUtils
from dipper import curie_map
from dipper.graph.RDFGraph import RDFGraph, URIRef


logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)


class EvidenceProvenanceTestCase(unittest.TestCase):

    def setUp(self):
        self.test_util = TestUtils()
        self.assoc_curie = 'MONARCH:test_association'
        self.eco_id = 'ECO:0000015'

        # Headers:
        # 01 marker_accession_id,
        # 02 marker_symbol,
        # 03 phenotyping_center,
        # 04 colony_raw,
        # 05 sex,
        # 06 zygosity,
        # 07 allele_accession_id,
        # 08 allele_symbol,
        # 09 allele_name,
        # 10 strain_accession_id,
        # 11 strain_name,
        # 12 project_name,
        # 13 project_fullname,
        # 14 pipeline_name,
        # 15 pipeline_stable_id,
        # 16 procedure_stable_id,
        # 17 procedure_name,
        # 18 parameter_stable_id,
        # 19 parameter_name,
        # 20 top_level_mp_term_id,
        # 21 top_level_mp_term_name,
        # 22 mp_term_id,
        # 23 mp_term_name,
        # 24 p_value,
        # 25 percentage_change,
        # 26 effect_size,
        # 27 statistical_method,
        # 28 resource_name

        self.test_set_1 = (
            'MGI:1920145',              # 01
            'Setd5',                    # 02
            'WTSI',                     # 03
            'MEFW',                     # 04
            'male',                     # 05
            'heterozygote',             # 06
            'MGI:4432631',              # 07
            'Setd5<tm1a(EUCOMM)Wtsi>',  # 08
            'targeted mutation 1a, Wellcome Trust Sanger Institute',    # 09
            'MGI:2159965',              # 10
            'C57BL/6N',                 # 11
            'MGP',                      # 12
            'Wellcome Trust Sanger Institute Mouse Genetics Project',   # 13
            'MGP Select Pipeline',      # 14
            'MGP_001',                  # 15
            'MGP_XRY_001',              # 16
            'X-ray',                    # 17
            'IMPC_XRY_008_001',         # 18
            'Number of ribs right',     # 19
            'MP:0005390',               # 20
            'skeleton phenotype',       # 21
            'MP:0000480',               # 22
            'increased rib number',     # 23
            '1.637023E-010',            # 24
            '',                         # 25
            '8.885439E-007',            # 26
            'Wilcoxon rank sum test with continuity correction',    # 27
            'IMPC'            # 28
        )

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

        (p_value, percentage_change, effect_size) = self.test_set_1[23:26]

        impc._add_evidence(
            self.assoc_curie, self.eco_id, p_value, percentage_change, effect_size,
            self.study_curie)

        triples = """
:MONARCH_test_association SEPIO:0000007 <https://monarchinitiative.org/.well-known/genid/b97a98087df7a99d8a38> .

<https://monarchinitiative.org/.well-known/genid/b97a98087df7a99d8a38> a ECO:0000015 ;
    SEPIO:0000084 <https://monarchinitiative.org/.well-known/genid/b41ad2bfd375c9de8888>,
        <https://monarchinitiative.org/.well-known/genid/b216606de82749b03956> ;
    SEPIO:0000085 <https://monarchinitiative.org/.well-known/genid/study> .

<https://monarchinitiative.org/.well-known/genid/b216606de82749b03956> a OBI:0000175 ;
    RO:0002353 <https://monarchinitiative.org/.well-known/genid/study> ;
    STATO:0000129 1.637023e-10 .

<https://monarchinitiative.org/.well-known/genid/b41ad2bfd375c9de8888> a STATO:0000085 ;
    RO:0002353 <https://monarchinitiative.org/.well-known/genid/study> ;
    STATO:0000129 "8.885439E-007" .
    
ECO:0000015 biolink:category biolink:EvidenceType .
OBI:0000175 biolink:category biolink:EvidenceType .
<https://monarchinitiative.org/.well-known/genid/b216606de82749b03956> biolink:category biolink:EvidenceType .
<https://monarchinitiative.org/.well-known/genid/b41ad2bfd375c9de8888> biolink:category biolink:EvidenceType .
<https://monarchinitiative.org/.well-known/genid/b97a98087df7a99d8a38> biolink:category biolink:EvidenceType .
<https://monarchinitiative.org/.well-known/genid/study> biolink:category biolink:EvidenceType .
<https://monarchinitiative.org/.well-known/genid/study> biolink:category biolink:InformationContentEntity .
<https://monarchinitiative.org/MONARCH_test_association> biolink:category biolink:Association .
    
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

        (phenotyping_center, colony) = self.test_set_1[2:4]
        (project_fullname, pipeline_name, pipeline_stable_id,
         procedure_stable_id, procedure_name, parameter_stable_id,
         parameter_name) = self.test_set_1[12:19]
        (statistical_method, resource_name) = self.test_set_1[26:28]

        impc._add_study_provenance(
            phenotyping_center, colony,
            project_fullname,
            pipeline_name, pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name)

        # dbg
        LOG.info(
            "Provenance graph as turtle:\n%s\n",
            impc.graph.serialize(format="turtle").decode("utf-8")
        )

        triples = """
<https://monarchinitiative.org/.well-known/genid/b0b26361b8687b5ad9ef> a owl:NamedIndividual ;
    rdfs:label "MEFW" ;
    biolink:category biolink:PopulationOfIndividualOrganisms .

<https://monarchinitiative.org/.well-known/genid/bdd05a8ca155ddaf415e> a OBI:0000471 ;
    BFO:0000050 <http://www.sanger.ac.uk/science/data/mouse-genomes-project>,
        IMPC-pipe:MGP_001 ;
    BFO:0000051 STATO:0000076,
        IMPC-proc:MGP_XRY_001 ;
    SEPIO:0000017 <http://www.sanger.ac.uk/> ;
    SEPIO:0000114 <https://www.mousephenotype.org/impress/OntologyInfo?action=list&procID=MGP_XRY_001#IMPC_XRY_008_001> ;
    biolink:category biolink:EvidenceType,
        biolink:InformationContentEntity .

STATO:0000076 biolink:category biolink:EvidenceType .

<http://www.sanger.ac.uk/> a foaf:organization ;
    rdfs:label "WTSI" ;
    biolink:category biolink:Provider .

<http://www.sanger.ac.uk/science/data/mouse-genomes-project> a VIVO:Project ;
    rdfs:label "Wellcome Trust Sanger Institute Mouse Genetics Project" ;
    biolink:category biolink:Provider .

<https://www.mousephenotype.org/impress/OntologyInfo?action=list&procID=MGP_XRY_001#IMPC_XRY_008_001> a owl:NamedIndividual ;
    rdfs:label "Number of ribs right (X-ray)" ;
    biolink:category biolink:EvidenceType,
        biolink:InformationContentEntity .

IMPC-pipe:MGP_001 a owl:NamedIndividual ;
    rdfs:label "MGP Select Pipeline" ;
    biolink:category biolink:EvidenceType .

IMPC-proc:MGP_XRY_001 a owl:NamedIndividual ;
    rdfs:label "X-ray" ;
    biolink:category biolink:EvidenceType .
"""

        # dbg
        LOG.info(
            "Reference graph: %s",
            impc.graph.serialize(format="turtle").decode("utf-8")
        )
        self.assertTrue(
            self.test_util.test_graph_equality(triples, impc.graph))

    def test_assertion_model(self):
        """
        Functional test for _add_study_provenance()
        """

        impc = IMPC('rdf_graph', True)
        impc.graph = RDFGraph(True)
        self.assertTrue(len(list(impc.graph)) == 0)

        impc._add_assertion_provenance(self.assoc_curie, self.evidence_curie)

        triples = """
    MONARCH:test_association SEPIO:0000015 <https://monarchinitiative.org/.well-known/genid/bf92df374a884963e805> .
    <https://monarchinitiative.org/.well-known/genid/bf92df374a884963e805> a SEPIO:0000001 ;
        SEPIO:0000018 <https://www.mousephenotype.org/> ;
        SEPIO:0000111 <https://monarchinitiative.org/.well-known/genid/evidence>  .

    <https://www.mousephenotype.org/> a foaf:organization ;
        rdfs:label "International Mouse Phenotyping Consortium" .

SEPIO:0000001 biolink:category biolink:InformationContentEntity .
<https://monarchinitiative.org/.well-known/genid/bf92df374a884963e805> biolink:category biolink:InformationContentEntity .
<https://monarchinitiative.org/.well-known/genid/evidence> biolink:category biolink:EvidenceType .
MONARCH:test_association biolink:category biolink:Association .
<https://www.mousephenotype.org/> biolink:category biolink:Provider .


        """
        # dbg
        LOG.info(
            "Assertion graph:\n %s\n", impc.graph.serialize(
                format="turtle").decode("utf-8")
        )

        self.assertTrue(self.test_util.test_graph_equality(triples, impc.graph))

    @unittest.skip("Timeouts on travis")
    def test_random_data_set(self):
        """
        Download dataset using fetch(), then take a row of data and
        run through evidence and provenance functions to test the output

        Line of data is hardcoded, but theoretically should work on any line
        """
        line_to_test = 1129
        count = 0
        impc = IMPC('rdf_graph', False)   # Not Skolem
        self.test_set_N = []
        # fetch file
        # impc.fetch(True)
        file_path = '/'.join((impc.rawdir, impc.files['all']['file']))
        with gzip.open(file_path, 'rt') as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for row in filereader:
                count += 1
                if count < line_to_test:
                    continue
                elif count == line_to_test:
                    self.test_set_N = row
                elif count > line_to_test:
                    LOG.info("stopped at line:\t%s\n", count)
                    break

        # Some DRY violation with the above tests
        (phenotyping_center, colony) = self.test_set_N[2:4]
        (project_fullname, pipeline_name, pipeline_stable_id,
         procedure_stable_id, procedure_name, parameter_stable_id,
         parameter_name) = self.test_set_N[12:19]
        (statistical_method, resource_name) = self.test_set_N[26:28]

        (p_value, percentage_change, effect_size) = self.test_set_N[23:26]

        # adding evidence
        impc._add_evidence(
            self.assoc_curie, self.eco_id, p_value, percentage_change, effect_size,
            self.study_curie)

        # adding  study
        impc._add_study_provenance(
            phenotyping_center, colony, project_fullname,
            pipeline_name,
            pipeline_stable_id,
            procedure_stable_id, procedure_name,
            parameter_stable_id, parameter_name,
            statistical_method, resource_name, line_to_test)

        # Note that this doesn't test much since we're dealing with
        # multiple part_of  and has_part links to individuals
        # which results in ambiguity = hard to test

        # dbg
        LOG.info(
            "Row %i graph as ntriples:\n%s\n",
            line_to_test, impc.graph.serialize(format="ntriples").decode("utf-8")
        )

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
        LOG.info(
            "Test that query for row %i passes and returns one row", int(line_to_test))

        # print("Sparql Output: %s\n", list(sparql_output) )
        # it is an array with one list with five vars in it

        self.assertEqual(len(list(sparql_output)), 1)

    def tearDown(self):
        return


if __name__ == '__main__':
    unittest.main()
