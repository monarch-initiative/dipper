#!/usr/bin/env python3
import unittest
import logging
import os
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.rdf2dot import rdf2dot
from dipper.sources.WormBase import WormBase

logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)

TEST_PATH = os.path.join(os.path.dirname(__file__), 'resources/wormbase')
DOT_PATH = TEST_PATH + "/dot/"
RAW_PATH = TEST_PATH + "/input/"

# Genes
GENES = [
    'WBGene00001414',
    'WBGene00004967',
    'WBGene00003916',
]


class WormBaseTestCase(unittest.TestCase):
    """
    Test WormBase methods

    Input source data is generated with scripts/test-sets/flybase.sh
    All tests here are functional, and will fail for a variety of reasons,
    including new or updated data connected to an allele

    Note setUp and tearDown do not run before/after TestCase.subTest !
    need something similar to @pytest.mark.parametrize

    """

    def setUp(self):
        self.wormbase = WormBase('rdf_graph', True)
        self.wormbase.graph = RDFGraph(True)

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
        Runs WormBase.parse() and outputs dot file for each allele
        This is less of a unit test and more for viewing the
        output of an entire run on a single allele,
        dot files can be converted to images using
        scripts/dot-to-svg.sh
        """
        graph_opts = {
            'rankdir': 'LR'
        }
        for gene in GENES:
            with self.subTest(gene_id=gene):
                self.tearDownAndSetUp()
                self.wormbase.rawdir = RAW_PATH + '/' + gene
                self.wormbase.version_num = 'test_version'
                self.wormbase.parse()
                dot_file_path = DOT_PATH + gene + ".dot"
                with open(dot_file_path, 'w') as dot_file:
                    rdf2dot(self.wormbase.graph, dot_file, graph_opts)


if __name__ == '__main__':
    unittest.main()
