#!/usr/bin/env python3

import unittest
from unittest.mock import mock_open
from unittest.mock import MagicMock
import logging
from dipper.sources.UDP import UDP

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class UDPTestCase(unittest.TestCase):
    """
    Test UDP parser
    """

    def setUp(self):
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
        udp = UDP()
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
        udp = UDP()
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
        mock_lines = [
            'patient_1\tHP:000001\tyes',
            'patient_1\tHP:000002\tno'
        ]
        mock_data = MagicMock()
        mock_data.__iter__.return_value = iter(mock_lines)

        mock_file = mock_open(mock=mock_data)
        udp = UDP()
        udp.load_bindings()
        udp._parse_patient_phenotypes(mock_file)
        sparql_query = """
            SELECT *
            WHERE {
                :patient_1 a foaf:Person ;
                    rdfs:label "patient_1" ;
                    OBO:RO_0002200 OBO:DOID_4,
                         OBO:HP_000001 .
            }
        """
        sparql_output = udp.graph.query(sparql_query)
        # Test that query passes and returns one row
        self.assertEqual(len(list(sparql_output)), 1)
