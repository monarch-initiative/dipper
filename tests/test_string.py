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
             0, 0, 0, 0, 300, 0, 150, 800]]

        # Test set with deprecated protein id
        self.test_set_2 = \
            [['9606.ENSP00000000233', '9606.ENSP00000006101',
              0, 0, 0, 0, 300, 0, 150, 800]]

        self.columns = [
            'protein1', 'protein2', 'neighborhood', 'fusion',
            'cooccurence', 'coexpression', 'experimental',
            'database', 'textmining', 'combined_score'
        ]

        ensembl = Ensembl('rdf_graph', True)
        self.protein_list = ensembl.fetch_protein_gene_map(9606)

        return

    def tearDown(self):
        return

    def testFakeDataSet1(self):
        string_db = StringDB('rdf_graph', True)
        string_db.graph.bind_all_namespaces()
        ensembl = Ensembl('rdf_graph', True)
        prot_map = ensembl.fetch_protein_gene_map(9606)
        for key in prot_map.keys():
            prot_map[key] = "ENSEMBL:{}".format(prot_map[key])

        print("Finished fetching ENSP IDs, "
              "fetched {} proteins".format(len(prot_map.keys())))
        dataframe = pd.DataFrame(data=self.test_set_1, columns=self.columns)

        string_db._process_protein_links(dataframe, prot_map, 9606)

        sparql_query = """
                      SELECT ?prot
                      WHERE {
                          ?prot RO:0002434 ENSEMBL:ENSG00000004059 .
                      }
                      """
        sparql_output = string_db.graph.query(sparql_query)
        results = list(sparql_output)
        expected = [(URIRef(string_db.graph._getNode("ENSEMBL:ENSG00000001626")),)]
        self.assertEqual(results, expected)

    def testFakeDataSet2(self):
        """
        Dataset contains a deprecated protein ID
        that we expect if filtered out by ensembl biomart
        We test that this returns a graph with 3 triples:
        MonarchData:string.ttl a owl:Ontology ;
        owl:versionIRI <https://archive.monarchinitiative.org/.../string.ttl> ;
        owl:versionInfo "some version"
        :return:
        """
        string_db = StringDB('rdf_graph', True)
        dataframe = pd.DataFrame(data=self.test_set_2, columns=self.columns)
        string_db._process_protein_links(dataframe, self.protein_list, 9606)
        self.assertEqual(len(string_db.graph), 3)

