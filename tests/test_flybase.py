#!/usr/bin/env python3
import unittest
import logging
import os
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.rdf2dot import rdf2dot
from dipper.utils.TestUtils import TestUtils
from dipper.sources.FlyBase import FlyBase

logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)

TEST_PATH = os.path.join(os.path.dirname(__file__), 'resources/flybase')
NT_PATH  = TEST_PATH + "/nt/"
DOT_PATH = TEST_PATH + "/dot/"
RAW_PATH = TEST_PATH + "/input/"
TTL_PATH = TEST_PATH + "/expected/"


# Alleles
alleles = [
    'FBal0195705',
    'FBal0256668'
]


class FlyBaseTestCase(unittest.TestCase):

    def test_parse(self):
        """
        Runs FlyBase.parse() and outputs dot file for each allele
        """
        for allele in alleles:
            with self.subTest(allele_id=allele):
                flybase = FlyBase('rdf_graph', True)
                flybase.graph = RDFGraph(True)
                flybase.rawdir = RAW_PATH + '/' + allele
                flybase.parse()
                dot_file_path = DOT_PATH + allele + ".dot"
                with open(dot_file_path, 'w') as dot_file:
                    rdf2dot(flybase.graph, dot_file)


if __name__ == '__main__':
    unittest.main()
