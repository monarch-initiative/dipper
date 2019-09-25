#!/usr/bin/env python3

import unittest
import logging
import rdflib
from dipper.utils import GraphUtils

logging.basicConfig(level=logging.WARNING)
LOG = logging.getLogger(__name__)


class GraphUtilsTestCase(unittest.TestCase):

    def setUp(self):
        self.graph_object = rdflib.Graph()
        self.graph_object2 = rdflib.Graph()
        self.graph_util = GraphUtils.GraphUtils
        self.test_graph = self.graph_object.parse(
            "tests/resources/graphutils/gu_test_graph.ttl", format="ttl")
        self.test_graph_id = self.test_graph.identifier
        self.test_graph2 = self.graph_object2.parse(
            "tests/resources/graphutils/gu_test_graph2.ttl", format="ttl")
        self.test_graph2_id = self.test_graph2.identifier
        self.enemyOf = rdflib.term.URIRef(
            'http://www.perceive.net/schemas/relationship/enemyOf')
        self.name = rdflib.term.URIRef('http://xmlns.com/foaf/0.1/name')

    def tearDown(self):
        pass

    def test_count_predicates(self):
        self.assertTrue(hasattr(self.graph_util.count_predicates, '__call__'),
                        "count_predicates isn't a method or doesn't exist")
        count = self.graph_util.count_predicates(self.test_graph)
        self.assertEqual(
            count.get(self.name), 1,
            "didn't get correct count for 'name'")
        self.assertEqual(
            count.get(self.enemyOf), 2,
            "didn't get correct count for 'enemy of'")

    def test_compare_graph_predicates(self):
        self.assertTrue(hasattr(self.graph_util.compare_graph_predicates,
                                '__call__'),
                        "compare_graph_predicates isn't a method " +
                        "or doesn't exist")

        compare = self.graph_util.compare_graph_predicates(
            self.test_graph,
            self.test_graph2)
        self.assertTrue(
            compare.get(self.enemyOf).get(str(self.test_graph_id)) == 2,
            "testing > 1 count, didn't get correct count for 'enemyOf'")
        self.assertTrue(
            compare.get(self.name).get(str(self.test_graph_id)) == 1,
            "testing hit on both graphs, " +
            "didn't get correct count for 'name' (graph 1)")
        self.assertTrue(
            compare.get(self.name).get(str(self.test_graph2_id)) == 1,
            "testing hit on both graphs, " +
            "didn't get correct count for 'name' (graph 2)")


if __name__ == '__main__':
    unittest.main()
