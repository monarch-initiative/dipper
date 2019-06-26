#!/usr/bin/env python3

import unittest
import logging
from dipper.graph.RDFGraph import RDFGraph
from dipper import curie_map

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class GeneralGraphTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = RDFGraph()
        self.curie_map = curie_map.get()

    def tearDown(self):
        self.graph = None

    def test_curieprefixes(self):
        """
        This will ensure that we can create identifiers for all of the
        defined curie prefixes using the GraphUtils.getNode() method
        :return:

        """
        # add one id per curie as classes to the graph
        for p in self.curie_map.keys():
            testid = p+':testme'
            n = self.graph._getnode(testid)
            m = "prefix \""+p+"\" has an error...can't create graph node"
            self.assertTrue(n is not None, m)

        return


if __name__ == '__main__':
    unittest.main()
