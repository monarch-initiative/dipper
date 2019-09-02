#!/usr/bin/env python3
import unittest
import logging
import os
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.rdf2dot import rdf2dot
from dipper.sources.GWASCatalog import GWASCatalog
from dipper.utils.TestUtils import TestUtils

logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)

TEST_PATH = os.path.join(os.path.dirname(__file__), 'resources/gwascatalog')
DOT_PATH = TEST_PATH + "/dot/"
RAW_PATH = TEST_PATH + "/input/"
TTL_PATH = TEST_PATH + "/expected/"


# Variants
VARIANTS = [
    'rs1329573',  # rs1329573-?; rs7020413-?; rs3824344-?; rs3758171-?
    'kgp8851185',
    'rs1491921',  # rs1491921-C
]


class GWASCatalogTestCase(unittest.TestCase):
    """
    Test gwas catalog methods

    Input source data is generated with scripts/test-sets/gwas.sh
    All tests here are functional, and will fail for a variety of reasons,
    including new or updated data connected to an allele

    Note setUp and tearDown do not run before/after TestCase.subTest !
    need something similar to @pytest.mark.parametrize

    """

    def setUp(self):
        self.gwascatalog = GWASCatalog('rdf_graph', True)
        self.gwascatalog.graph = RDFGraph(True)

    def tearDown(self):
        self.gwascatalog = None

    def tearDownAndSetUp(self):
        """
        In lieu of the setUp and tearDown working for subTests
        """
        self.tearDown()
        self.setUp()

    def test_parse(self):
        """
        Runs WormBase.parse() and outputs dot file for each allele
        This is less of a unit test and more for viewing the
        output of an entire run on a single allele,
        dot files can be converted to images using
        scripts/dot-to-svg.sh
        """
        for variant in VARIANTS:
            with self.subTest(variant_id=variant):
                self.tearDownAndSetUp()
                self.gwascatalog.rawdir = RAW_PATH + '/' + variant
                self.gwascatalog.parse()
                dot_file_path = DOT_PATH + variant + ".dot"
                with open(dot_file_path, 'w') as dot_file:
                    rdf2dot(self.gwascatalog.graph, dot_file)

                # debug
                LOG.debug(
                    "Reference graph: %s",
                    self.gwascatalog.graph.serialize(format="turtle").decode("utf-8"))

                reference_ttl = TTL_PATH + variant + '.ttl'

                self.assertTrue(TestUtils.test_graph_equality(
                    reference_ttl, self.gwascatalog.graph))


if __name__ == '__main__':
    unittest.main()
