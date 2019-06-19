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
    """
    Test FlyBase methods

    Note setUp and tearDown do not run before/after TestCase.subTest !
    need something similar to @pytest.mark.parametrize
    """

    def setUp(self):
        self.flybase = FlyBase('rdf_graph', True)
        self.flybase.graph = RDFGraph(True)

    def tearDown(self):
        self.flybase = None

    def test_parse(self):
        """
        Runs FlyBase.parse() and outputs dot file for each allele
        This is less of a unit test and more for viewing the
        output of an entire run on a single allele,
        dot files can be converted to images using
        scripts/dot-to-svg.sh
        """
        graph_opts = {
            'rankdir': 'LR'
        }
        for allele in alleles:
            with self.subTest(allele_id=allele):
                self.tearDown()
                self.setUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase.parse()
                dot_file_path = DOT_PATH + allele + ".dot"
                with open(dot_file_path, 'w') as dot_file:
                    rdf2dot(self.flybase.graph, dot_file, graph_opts)

    def test_allele_gene(self):
        """
        test FlyBase._process_allele_gene()
        """
        for allele in alleles:
            with self.subTest(allele_id=allele):
                self.tearDown()
                self.setUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_allele_gene(limit=None)

                reference_ttl = TTL_PATH + allele + '/' + 'allele_gene.ttl'
                with open(reference_ttl, 'r') as ref_fh:
                    ref_graph = "\n".join(ref_fh.readlines())

                # debug
                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                self.assertTrue(TestUtils.test_graph_equality(
                    ref_graph, self.flybase.graph))

    def test_disease_model(self):
        """
        test FlyBase._process_disease_model()
        """
        for allele in alleles:
            with self.subTest(allele_id=allele):
                self.tearDown()
                self.setUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_disease_model(limit=None)

                reference_ttl = TTL_PATH + allele + '/' + 'disease_model.ttl'
                with open(reference_ttl, 'r') as ref_fh:
                    ref_graph = "\n".join(ref_fh.readlines())

                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                self.assertTrue(TestUtils.test_graph_equality(
                    ref_graph, self.flybase.graph))

    def test_allele_phenotype(self):
        """
        test FlyBase._process_allele_phenotype()
        """
        for allele in alleles:
            with self.subTest(allele_id=allele):
                self.tearDown()
                self.setUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_allele_phenotype(limit=None)

                reference_ttl = TTL_PATH + allele + '/' + 'allele_phenotype.ttl'
                with open(reference_ttl, 'r') as ref_fh:
                    ref_graph = "\n".join(ref_fh.readlines())

                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                self.assertTrue(TestUtils.test_graph_equality(
                    ref_graph, self.flybase.graph))

    def test_gene_xref(self):
        """
        test FlyBase._process_gene_xref()
        """
        for allele in alleles:
            with self.subTest(allele_id=allele):
                self.tearDown()
                self.setUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_gene_xref(limit=None)

                reference_ttl = TTL_PATH + allele + '/' + 'gene_xref.ttl'
                with open(reference_ttl, 'r') as ref_fh:
                    ref_graph = "\n".join(ref_fh.readlines())

                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                self.assertTrue(TestUtils.test_graph_equality(
                    ref_graph, self.flybase.graph))


if __name__ == '__main__':
    unittest.main()
