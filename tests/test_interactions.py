#!/usr/bin/env python3

import unittest
import logging
from rdflib import Graph
from dipper.sources.CTD import CTD
from dipper import curie_map

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class InteractionsTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = Graph()

        self.curie_map = curie_map.get()
        self.ctd = CTD('rdf_graph', True)
        self.ctd.g = self.ctd.graph
        row1 = [
            '06-Paris-LA-66 protocol', 'C046983', 'foo',
            'Precursor Cell Lymphoblastic Leukemia-Lymphoma', 'MESH:D054198',
            'therapeutic', 'bar', 'baz', 'foo', '4519131']
        row2 = [
            '10,10-bis(4-pyridinylmethyl)-9(10H)-anthracenone', 'C112297',
            'foo', 'Hyperkinesis', 'MESH:D006948', 'marker/mechanism',
            'bar', 'baz', 'foo', '19098162']

        self.ctd._process_interactions(row1)
        self.ctd._process_interactions(row2)

    def tearDown(self):
        self.ctd.graph = None
        self.ctd = None

    def test_therapeutic_relationship(self):
        from dipper.utils.TestUtils import TestUtils
        from dipper.utils.GraphUtils import GraphUtils

        # Make testutils object and load bindings
        test_query = TestUtils(self.ctd.graph)

        # Expected structure
        sparql_query = """
                       SELECT ?assoc ?pubmed ?disease ?chemical
                       WHERE {
                       ?assoc a Annotation: ;
                           dc:evidence OBO:ECO_0000033 ;
                           dc:source ?pubmed ;
                           :hasObject ?disease ;
                           :hasPredicate OBO:RO_0002606 ;
                           :hasSubject ?chemical .}
                       """

        # SPARQL variables to check
        chem_id = 'MESH:D009538'
        chem_uri = self.graph._getNode(chem_id)
        disease_id = 'OMIM:188890'
        disease_uri = self.graph._getNode(disease_id)
        pubmed_id = 'PMID:16785264'
        pubmed_uri = self.graph._getNode(pubmed_id)
        rel_id = self.model.object_properties['substance_that_treats']
        eco = 'ECO:0000033'
        # TODO PYLINT  make_association_id() does not exist in CTD
        # there is "_make_association()" with a different sig

        assoc_id = self.ctd.make_association_id(
            'ctd', chem_id, rel_id, disease_id, eco, pubmed_id)
        assoc_uri = self.graph._getNode(assoc_id)

        # Expected output from query
        expected_output = [assoc_uri, pubmed_uri, disease_uri, chem_uri]

        # Query graph
        sparql_output = test_query.query_graph(sparql_query)

        self.assertTrue(expected_output in sparql_output)

        logger.info("Test finished.")


if __name__ == '__main__':
    unittest.main()
