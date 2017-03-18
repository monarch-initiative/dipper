#!/usr/bin/env python3

import unittest
from dipper.sources.StringDB import StringDB
from dipper.sources.Ensembl import Ensembl
import pandas as pd
from rdflib import URIRef


class StringTestFakeData(unittest.TestCase):

    def setUp(self):
        # Test set with two proteins from same species
        self.test_set_1 = \
            [['9606.ENSP00000000233', '9606.ENSP00000003084',
             0, 0, 0, 0, 50, 0, 150, 150]]

        # Test set with deprecated protein id
        self.test_set_2 = \
            [['9606.ENSP00000000233', '9606.ENSP00000006101',
             0, 0, 0, 0, 50, 0, 150, 150]]

        self.columns = [
            'protein1', 'protein2', 'neighborhood', 'fusion',
            'cooccurence', 'coexpression', 'experimental',
            'database', 'textmining', 'combined_score'
        ]

        ensembl = Ensembl('rdf_graph', True)
        self.protein_list = ensembl.fetch_protein_list(9606)

        return

    def tearDown(self):
        return

    def testFakeDataSet1(self):
        string_db = StringDB('rdf_graph', True)
        string_db.graph.bind_all_namespaces()
        ensembl = Ensembl('rdf_graph', True)
        protein_list = ensembl.fetch_protein_list(9606)
        dataframe = pd.DataFrame(data=self.test_set_1, columns=self.columns)

        string_db._process_protein_links(dataframe, protein_list, 9606)

        sparql_query = """
                      SELECT ?prot
                      WHERE {
                          ?prot RO:0002434 ENSEMBL:ENSP00000003084 .
                          ENSEMBL:ENSP00000003084 RO:0002434 ?prot .
                      }
                      """
        sparql_output = string_db.graph.query(sparql_query)
        results = list(sparql_output)
        expected = [(URIRef(string_db.graph._getNode("ENSEMBL:ENSP00000000233")),)]
        self.assertEqual(results, expected)

    def testFakeDataSet2(self):
        string_db = StringDB('rdf_graph', True)
        dataframe = pd.DataFrame(data=self.test_set_2, columns=self.columns)
        string_db._process_protein_links(dataframe, self.protein_list, 9606)
        self.assertEqual(len(string_db.graph), 3)

