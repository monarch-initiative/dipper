#!/usr/bin/env python3

from dipper.sources.CTD import CTD
from dipper import curie_map
from rdflib import Graph

import unittest
import logging

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


class InteractionsTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = Graph()
        self.curie_map = curie_map.get()
        self.ctd = CTD()
        self.ctd.graph = Graph()

        row1 = ['06-Paris-LA-66 protocol', 'C046983', 'foo',
               'Precursor Cell Lymphoblastic Leukemia-Lymphoma',
               'MESH:D054198', 'therapeutic', 'bar', 'baz', 'foo', '4519131']
        row2 = ['10,10-bis(4-pyridinylmethyl)-9(10H)-anthracenone',
                'C112297', 'foo', 'Hyperkinesis', 'MESH:D006948',
                'marker/mechanism', 'bar', 'baz', 'foo', '19098162']

        pub_map = {}
        self.ctd._process_interactions(row1, pub_map)
        self.ctd._process_interactions(row2, pub_map)

    def tearDown(self):
        self.ctd.graph = None
        self.ctd = None

    def test_therapeutic_relationship(self):
        from dipper.utils.TestUtils import TestUtils
        from dipper.utils.CurieUtil import CurieUtil
        from rdflib.namespace import URIRef

        cu = CurieUtil(self.curie_map)
        chem_id = 'C046983'
        curie = 'MESH:'+chem_id
        chem_uri = URIRef(cu.get_uri(curie))
        disease_id = 'MESH:D054198'
        disease_uri = URIRef(cu.get_uri(disease_id))
        relationship = 'MONARCH:treats'
        relationship_uri = URIRef(cu.get_uri(relationship))

        self.assertTrue((chem_uri, relationship_uri, disease_uri)
                        in self.ctd.graph)

if __name__ == '__main__':
    unittest.main()