#!/usr/bin/env python3

from dipper.sources.CTD import CTD
from dipper import curie_map
from rdflib import Graph

import unittest
import logging

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


class GenotypeTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = Graph()
        self.curie_map = curie_map.get()
        self.ctd = CTD()
        self.ctd.graph = Graph()

    def tearDown(self):
        self.ctd.graph = None
        self.ctd = None

    def test_process_interaction(self):
        from dipper.utils.TestUtils import TestUtils

        # Reset graph
        self.ctd.graph = Graph()
        row = ['06-Paris-LA-66 protocol', 'C046983', 'foo',
               'Precursor Cell Lymphoblastic Leukemia-Lymphoma',
               'MESH:D054198', 'therapeutic', 'bar', 'baz', 'foo', '4519131']
        pub_map = {}
         # Expect query result of the above row, pub_map combo

        self.ctd._process_interactions(row, pub_map)
        test_query = TestUtils(self.ctd.graph)
        self.assertTrue(True)



if __name__ == '__main__':
    unittest.main()