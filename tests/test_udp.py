#!/usr/bin/env python3

import unittest
from unittest.mock import mock_open
from unittest.mock import MagicMock
import logging
from dipper.sources.UDP import UDP
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph


logging.basicConfig()
logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


class UDPTestCase(unittest.TestCase):
    """
    Test UDP parser
    """

    def setUp(self):
        self.test_util = TestUtils()
        return

    def tearDown(self):
        return

    def test_dbsnp_indel_resolution(self):
        """
        unit test for _get_rs_id()
        Test that we can resolve indels that
        have different insertion sequence(s)
        for one rsid
        15	51766637	374313651	in-del	-/A/AA/AAA/AAAA/CAAA/TAAA
        """
        udp = UDP('rdf_graph', True)
        rs_map = udp._parse_rs_map_file(udp.map_files['dbsnp_map'])
        variant_type = 'indel'
        variant = {
            'build': 'hg19',
            'chromosome': 'chr15',
            'reference_allele': '-',
            'variant_allele': 'AAAA',
            'position': '51766637'
        }
        rsid = udp._get_rs_id(variant, rs_map, variant_type)

        self.assertEqual(rsid, '374313651')

    def test_dbsnp_snp_mapping(self):
        """
        unit test for _get_rs_id()
        Test that we can resolve snps in dbsnp
        to rsids
        """
        udp = UDP('rdf_graph', True)
        rs_map = udp._parse_rs_map_file(udp.map_files['dbsnp_map'])
        variant_type = 'snp'
        variant = {
            'build': 'hg19',
            'chromosome': 'chr15',
            'reference_allele': 'A',
            'variant_allele': 'C',
            'position': '54624219'
        }
        rsid = udp._get_rs_id(variant, rs_map, variant_type)

        self.assertEqual(rsid, '755532609')

    def test_patient_phenotype_model(self):
        """
        functional test for _parse_patient_phenotypes()
        """
        udp = UDP('rdf_graph', True)
        udp.graph = RDFGraph(True)

        # test that graph is empty
        self.assertTrue(len(list(udp.graph)) == 0)

        mock_lines = [
            'patient_1\tHP:000001\tyes',
            'patient_1\tHP:000002\tno'
        ]
        mock_data = MagicMock()
        mock_data.__iter__.return_value = iter(mock_lines)

        mock_file = mock_open(mock=mock_data)
        udp._parse_patient_phenotypes(mock_file)
        triples = """
        <https://monarchinitiative.org/MONARCH_patient_1> a foaf:Person ;
            rdfs:label "patient_1" ;
            RO:0002200 DOID:4,
              HP:000001 .
        """

        self.assertTrue(self.test_util.test_graph_equality(
            triples, udp.graph))

    def test_variant_model(self):
        """
        functional test for _parse_patient_variants()
        """
        udp = UDP('rdf_graph', True)
        udp.graph = RDFGraph(True)
        # test that graph is empty
        self.assertTrue(len(list(udp.graph)) == 0)

        data = ['patient_1',
                'family_1',
                '1',
                'HG19',
                '155230432',
                'G',
                'A',
                'Maternal',
                'Biallelic',
                'Non-synonymous;DOWNSTREAM',
                'CLK2',
                '',
                '',
                '',
                '',
                '',
                '',
                '',
                'Compound heterozygous',
                'Heterozygous',
                '',
                '0.002747253',
                '']
        test_data = "\t".join(data)
        mock_lines = [test_data]
        mock_data = MagicMock()
        mock_data.__iter__.return_value = iter(mock_lines)

        mock_file = mock_open(mock=mock_data)

        udp._parse_patient_variants(mock_file)

        triples = """
        <https://monarchinitiative.org/MONARCH_patient_1> GENO:0000222 <https://monarchinitiative.org/.well-known/genid/ba5f377fc8c95d4a6d7a> .

        <https://monarchinitiative.org/.well-known/genid/b41e8da0787b45e24c4f> a SO:0001059 ;
            rdfs:label "hg19chr1(CLK2):g.155230432G>A" ;
            GENO:0000418 HGNC:2069 ;
            RO:0002162 NCBITaxon:9606 ;
            owl:sameAs dbSNP:rs11557757 .

        <https://monarchinitiative.org/.well-known/genid/ba5f377fc8c95d4a6d7a> a GENO:0000000 ;
            rdfs:label "patient_1 genotype" ;
            GENO:0000382 <https://monarchinitiative.org/.well-known/genid/b41e8da0787b45e24c4f> .
        """

        self.assertTrue(self.test_util.test_graph_equality(triples, udp.graph))
