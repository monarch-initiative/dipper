#!/usr/bin/env python3

import unittest
from dipper.sources.StringDB import StringDB
from dipper.sources.Ensembl import Ensembl
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.TestUtils import TestUtils
import pandas as pd


class StringTestFakeData(unittest.TestCase):

    def setUp(self):
        self.test_util = TestUtils()
        # Test set with two proteins from same species
        self.test_set_1 = [[
            '9606.ENSP00000000233', '9606.ENSP00000003084',
            0, 0, 0, 0, 300, 0, 150, 800]]

        # Test set with deprecated protein id
        self.test_set_2 = [[
            '9606.ENSP00000000233', '9606.ENSP00000006101',
            0, 0, 0, 0, 300, 0, 150, 800]]

        self.columns = [
            'protein1', 'protein2', 'neighborhood', 'fusion', 'cooccurence',
            'coexpression', 'experimental', 'database', 'textmining', 'combined_score']

        ensembl = Ensembl('rdf_graph', True)
        self.protein_list = ensembl.fetch_protein_gene_map('9606')

        return

    def tearDown(self):
        return

    def testFakeDataSet1(self):
        string_db = StringDB('rdf_graph', True)
        string_db.graph = RDFGraph(True)
        self.assertEqual(len(string_db.graph), 0)

        ensembl = Ensembl('rdf_graph', True)
        prot_map = ensembl.fetch_protein_gene_map('9606')

        [prot_map.update({k: ['ENSEMBL:' + prot_map[k]]}) for k in prot_map.keys()]

        print("Finished fetching ENSP IDs, fetched {} proteins".format(len(prot_map)))

        # just looking
        # for key in prot_map:
        #    if string_db.graph.curie_regexp.match(prot_map[key]) is None:
        #        print("INVALID curie for %s from %s", prot_map[key], key)

        dataframe = pd.DataFrame(data=self.test_set_1, columns=self.columns)

        string_db._process_protein_links(dataframe, prot_map, '9606')

        # g1 <interacts with> g2
        triples = """
ENSEMBL:ENSG00000001626 RO:0002434 ENSEMBL:ENSG00000004059 .
ENSEMBL:ENSG00000001626  biolink:category biolink:Gene .
ENSEMBL:ENSG00000004059 biolink:category biolink:Gene .
        """

        self.assertTrue(self.test_util.test_graph_equality(triples, string_db.graph))

    def testFakeDataSet2(self):
        """
        Dataset contains a deprecated protein ID
        that we expect if filtered out by ensembl biomart
        We test that this returns an empty graph
        :return:
        """
        string_db = StringDB('rdf_graph', True)
        string_db.graph = RDFGraph()
        self.assertEqual(len(string_db.graph), 0)

        dataframe = pd.DataFrame(data=self.test_set_2, columns=self.columns)
        string_db._process_protein_links(dataframe, self.protein_list, '9606')
        self.assertEqual(len(string_db.graph), 0)
