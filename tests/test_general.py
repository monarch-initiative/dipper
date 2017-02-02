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
            n = self.graph._getNode(testid)
            m = "prefix \""+p+"\" has an error...can't create graph node"
            self.assertTrue(n is not None, m)

        return

    def readGraphFromTurtleFile(self, f):
        """
        This will read the specified file into a graph.  A simple parsing test.
        :param f:
        :return:

        """
        import os
        vg = RDFGraph()
        p = os.path.abspath(f)
        logger.info("Testing reading turtle file from %s", p)
        vg.parse(f, format="turtle")
        logger.info('Found %s graph nodes in %s', len(vg), p)
        self.assertTrue(len(vg) > 0, "No nodes found in "+p)

        return

    def readGraphIntoOWL(self, f):
        """
        test if the ttl can be parsed by owlparser
        this expects owltools to be accessible from commandline
        :param f: file of ttl
        :return:
        """

        import subprocess
        from subprocess import check_call

        status = check_call(["owltools", f], stderr=subprocess.STDOUT)
        # returns zero is success!
        if status != 0:
            logger.error(
                'finished verifying with owltools with status %s', status)
        self.assertTrue(status == 0)

        return

if __name__ == '__main__':
    unittest.main()
