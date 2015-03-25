#!/usr/bin/env python3

from dipper.sources.CTD import CTD
from dipper import curie_map
from rdflib import Graph

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
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

        # Make testutils object and load bindings
        test_query = TestUtils(self.ctd.graph)
        self.ctd.load_bindings()

        # Expected structure
        sparql_query = """
                       SELECT ?assoc ?pubmed ?disease ?chemical
                       WHERE {
                       ?assoc a Annotation: ;
                           dc:evidence OBO:ECO_0000033 ;
                           dc:source ?pubmed ;
                           :hasObject ?disease ;
                           :hasPredicate :MONARCH_treats ;
                           :hasSubject ?chemical .}
                       """

        # SPARQL variables to check
        cu = CurieUtil(self.curie_map)
        direct_evidence = 'therapeutic'
        chem_id = 'MESH:C046983'
        chem_uri = URIRef(cu.get_uri(chem_id))
        disease_id = 'MESH:D054198'
        disease_uri = URIRef(cu.get_uri(disease_id))
        assoc_id = self.ctd.make_id('ctd' + chem_id.replace('MESH:', '')
                                    + disease_id + direct_evidence)
        assoc_uri = URIRef(cu.get_uri(assoc_id))
        pubmed_id = 'PMID:4519131'
        pubmed_uri = URIRef(cu.get_uri(pubmed_id))

        # Expected output from query
        expected_output = [[assoc_uri, pubmed_uri, disease_uri, chem_uri]]

        # Query graph
        sparql_output = test_query.query_graph(sparql_query)

        self.assertEqual(expected_output, sparql_output)

if __name__ == '__main__':
    unittest.main()