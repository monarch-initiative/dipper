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
NT_PATH = TEST_PATH + "/nt/"
DOT_PATH = TEST_PATH + "/dot/"
RAW_PATH = TEST_PATH + "/input/"
TTL_PATH = TEST_PATH + "/expected/"

# Alleles
ALLELES = [
    'FBal0195705',
    'FBal0256668'
]


class FlyBaseTestCase(unittest.TestCase):
    """
    Test FlyBase methods

    Input source data is generated with scripts/test-sets/flybase.sh
    All tests here are functional, and will fail for a variety of reasons,
    including new or updated data connected to an allele

    Note setUp and tearDown do not run before/after TestCase.subTest !
    need something similar to @pytest.mark.parametrize

    Most of these tests are repeated code, should determine
    a way to dry it out
    """

    def setUp(self):
        self.flybase = FlyBase('rdf_graph', True)
        self.flybase.graph = RDFGraph(True)

    def tearDown(self):
        self.flybase = None

    def tearDownAndSetUp(self):
        """
        In lieu of the setUp and tearDown working for subTests
        """
        self.tearDown()
        self.setUp()

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
        for allele in ALLELES:
            with self.subTest(allele_id=allele):
                self.tearDownAndSetUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase.parse()
                dot_file_path = DOT_PATH + allele + ".dot"
                with open(dot_file_path, 'w') as dot_file:
                    rdf2dot(self.flybase.graph, dot_file, graph_opts)

    def test_allele_gene(self):
        """
        test FlyBase._process_allele_gene()
        """
        for allele in ALLELES:
            with self.subTest(allele_id=allele):
                self.tearDownAndSetUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_allele_gene(limit=None)

                # debug
                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                reference_ttl = TTL_PATH + allele + '/' + 'allele_gene.ttl'

                self.assertTrue(TestUtils.test_graph_equality(
                    reference_ttl, self.flybase.graph))

    def test_disease_model(self):
        """
        test FlyBase._process_disease_model()
        """
        for allele in ALLELES:
            with self.subTest(allele_id=allele):
                self.tearDownAndSetUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_disease_model(limit=None)

                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                reference_ttl = TTL_PATH + allele + '/' + 'disease_model.ttl'

                self.assertTrue(TestUtils.test_graph_equality(
                    reference_ttl, self.flybase.graph))

    def test_allele_phenotype(self):
        """
        test FlyBase._process_allele_phenotype()
        """
        for allele in ALLELES:
            with self.subTest(allele_id=allele):
                self.tearDownAndSetUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_allele_phenotype(limit=None)

                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                reference_ttl = TTL_PATH + allele + '/' + 'allele_phenotype.ttl'

                self.assertTrue(TestUtils.test_graph_equality(
                    reference_ttl, self.flybase.graph))

    def test_gene_xref(self):
        """
        test FlyBase._process_gene_xref()
        """
        for allele in ALLELES:
            with self.subTest(allele_id=allele):
                self.tearDownAndSetUp()
                self.flybase.rawdir = RAW_PATH + '/' + allele
                self.flybase._process_gene_xref(limit=None)
                LOG.debug(
                    "Reference graph: %s",
                    self.flybase.graph.serialize(format="turtle").decode("utf-8"))

                reference_ttl = TTL_PATH + allele + '/' + 'gene_xref.ttl'
                self.assertTrue(TestUtils.test_graph_equality(
                    reference_ttl, self.flybase.graph))


if __name__ == '__main__':
    unittest.main()
