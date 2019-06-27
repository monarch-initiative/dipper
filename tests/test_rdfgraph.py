#!/usr/bin/env python3

import os
import unittest
import logging
from rdflib import URIRef
from dipper.graph.RDFGraph import RDFGraph

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class RDFGraphTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = RDFGraph()

        self.test_cat_subj = "http://www.google.com"
        self.test_cat_default_pred = "http://w3id.org/biolink/vocab/category"
        self.test_cat_default_category = "http://w3id.org/biolink/vocab/NamedThing"

        self.test_cat_nondefault_category = "http://w3id.org/biolink/vocab/Gene"

    def tearDown(self):
        self.graph = None

    def test_add_triple_makes_triple(self):
        """
        test that addTriple() makes at least one triple
        """
        self.graph.addTriple(subject_id=self.test_cat_subj,
                             predicate_id="rdf:type",
                             obj="rdf:class")
        self.assertTrue(len(self.graph) > 0, "addTriples() didn't make >=1 triple")

    def test_add_triple_subject_category_assignment(self):
        """
        test that addTriple() correctly assigns subject category
        """
        self.graph.addTriple(subject_id=self.test_cat_subj,
                             predicate_id="rdf:comment",
                             obj="website",
                             subject_category=self.test_cat_nondefault_category)
        self.assertEqual(len(self.graph), 2,
                         "addTriples() didn't make exactly two triples " +
                         "(should be one for the triple itself" +
                         "and one for the subject category)")
        for this_subj, this_pred, this_obj in self.graph.triples(
                (URIRef(self.test_cat_subj), URIRef(self.test_cat_default_pred), None)):
            self.assertEqual(URIRef(self.test_cat_nondefault_category), this_obj)
            break

    def test_add_triple_object_category_assignment(self):
        """
        test that addTriple() correctly assigns obj category
        """
        self.graph.addTriple(subject_id=self.test_cat_subj,
                             predicate_id="rdf:type",
                             obj="rdf:class",
                             object_category=self.test_cat_nondefault_category)
        self.assertEqual(2, len(self.graph),
                         "addTriples() didn't make exactly two triples " +
                         "(should be one for the triple itself" +
                         "and one for the object category)")

        for this_subj, this_pred, this_obj in self.graph.triples(
                (URIRef(self.test_cat_subj), URIRef(self.test_cat_default_pred), None)):
            self.assertEqual(URIRef(self.test_cat_nondefault_category), this_obj)
            break

    def read_graph_from_turtle_file(self, f):
        """
        This will read the specified file into a graph.  A simple parsing test.
        :param f:
        :return:

        """
        vg = RDFGraph()
        p = os.path.abspath(f)
        logger.info("Testing reading turtle file from %s", p)
        vg.parse(f, format="turtle")
        logger.info('Found %s graph nodes in %s', len(vg), p)
        self.assertTrue(len(vg) > 0, "No nodes found in "+p)

        return

    def read_graph_into_owl(self, f):
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


    def test_make_category_triple_default(self):
        """
        test that method adds category triple to graph correctly (default pred and obj)
        """
        self.graph._make_category_triple(self.test_cat_subj)

        self.assertEqual(len(self.graph), 1, "method didn't make a triple")
        for this_subj, this_pred, this_obj in self.graph.triples((None, None, None)):
            self.assertEqual(URIRef(self.test_cat_subj), this_subj)
            self.assertEqual(URIRef(self.test_cat_default_pred), this_pred)
            self.assertEqual(URIRef(self.test_cat_default_category), this_obj)
            break

    def test_make_category_triple_nondefault_category(self):
        """
        test that method adds category triple to graph correctly
        """
        category = "http://w3id.org/biolink/vocab/gene"
        self.graph._make_category_triple(self.test_cat_subj, category)
        self.assertEqual(len(self.graph), 1, "method didn't make a triple")
        for this_subj, this_pred, this_obj in self.graph.triples((None, None, None)):
            self.assertEqual(URIRef(category), this_obj)
            break

    def test_make_category_triple_nondefault_pred(self):
        """
        test that method adds category triple to graph correctly (non default pred)
        """
        nondefault_pred = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
        self.graph._make_category_triple(self.test_cat_subj, self.test_cat_default_category,
                                         predicate=nondefault_pred)
        self.assertEqual(len(self.graph), 1, "method didn't make a triple")
        for this_subj, this_pred, this_obj in self.graph.triples((None, None, None)):
            self.assertEqual(URIRef(nondefault_pred), this_pred)
            break

    def test_make_category_triple_category_none_should_emit_named_thing(self):
        """
        test that method adds category triple to graph correctly (default pred and obj)
        """
        self.graph._make_category_triple(self.test_cat_subj, category=None)

        self.assertEqual(len(self.graph), 1, "method didn't make a triple")
        for this_subj, this_pred, this_obj in self.graph.triples((None, None, None)):
            self.assertEqual(URIRef(self.test_cat_default_category), this_obj)
            break

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