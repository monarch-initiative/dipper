#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.UDP import UDP
from rdflib import URIRef, BNode
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class UDPTestCase(unittest.TestCase):

    def setUp(self):

        self.test_set_1 = ('')
        return

    def tearDown(self):
        return

    def test_dbsnp_indel_resolution(self):
        """
        Functional test for _get_rs_id()
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
        Functional test for _get_rs_id()
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