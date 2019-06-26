#!/usr/bin/env python3

import unittest
import logging
from rdflib import URIRef
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

    def test_make_category_triple(self):
        """
        test that method adds category triple to graph correctly
        """
        subj = "http://www.google.com"
        pred = "http://w3id.org/biolink/vocab/category"
        obj = "http://w3id.org/biolink/vocab/namedThing"

        self.graph._make_category_triple(subj, obj)
        self.assertEqual(len(self.graph), 1, "method didn't make a triple")

        for this_subj, this_pred, this_obj in self.graph.triples((None, None, None)):
            self.assertEqual(URIRef(subj), this_subj)
            self.assertEqual(URIRef(pred), this_pred)
            self.assertEqual(URIRef(this_obj), this_obj)
            continue

        self.graph._make_category_triple(subj, obj, predicate="rdf:type")
        self.assertEqual(len(self.graph), 2, "method didn't make a triple")
        for this_subj, this_pred, this_obj in self.graph.triples(
                (None, URIRef("rdf:type"), None)):
            self.assertEqual(URIRef(pred), this_pred)


    def test_is_literal(self):
        """
        test that method infers type (either literal or CURIE) correctly
        """
        self.assertEqual(self.graph._is_literal("1"), True)
        self.assertEqual(self.graph._is_literal("foo:bar"), False)
        self.assertEqual(self.graph._is_literal("http://www.zombo.com/"), False)
        self.assertEqual(self.graph._is_literal("https://www.zombo.com/"), False)
        self.assertEqual(self.graph._is_literal("ftp://ftp.1000genomes.ebi.ac.uk/"),
                                                False)


if __name__ == '__main__':
    unittest.main()
