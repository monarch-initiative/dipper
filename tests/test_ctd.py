#!/usr/bin/env python3

from dipper.sources.CTD import CTD
from tests.test_source import SourceTestCase

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class CTDTestCase(SourceTestCase):
    def setUp(self):
        self.source = CTD()
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


    def test_therapeutic_relationship(self):
        from dipper.utils.TestUtils import TestUtils
        from dipper.utils.GraphUtils import GraphUtils
        from dipper import curie_map

        # Make testutils object and load ttl
        test_query = TestUtils(self.source.graph)
        test_query.load_testgraph_from_turtle(self.source)

        # Expected structure
        # TODO can this be unified OBAN and the Annot models to be automatically generated?
        sparql_query = """
                       SELECT ?assoc ?pubmed ?disease ?chemical
                       WHERE {
                       ?assoc a Annotation: ;
                           dc:evidence OBO:ECO_0000033 ;
                           dc:source ?pubmed ;
                           :hasObject ?disease ;
                           :hasPredicate :MONARCH_treats ;
                           :hasSubject ?chemical .}
                       """

        # SPARQL variables to check
        gu = GraphUtils(curie_map.get())
        direct_evidence = 'therapeutic'
        chem_id = 'MESH:D009538'
        chem_uri = gu.getNode(chem_id)
        disease_id = 'OMIM:188890'
        disease_uri = gu.getNode(disease_id)
        # consider replacing with make_ctd_chem_disease_assoc_id()
        assoc_id = self.source.make_id('ctd' + chem_id.replace('MESH:', '') #MONARCH_ed73d862405982b83a9a016c2b02ba71
                                    + disease_id + direct_evidence)
        assoc_uri = gu.getNode(assoc_id)
        pubmed_id = 'PMID:16785264'
        pubmed_uri = gu.getNode(pubmed_id)

        # One of the expected outputs from query
        expected_output = [assoc_uri, pubmed_uri, disease_uri, chem_uri]

        # Query graph
        sparql_output = test_query.query_graph(sparql_query)

        self.assertTrue(expected_output in sparql_output, "did not find expected association: "+assoc_id)

        logger.info("Test query data finished.")

    # OTHER CTD-SPECIFIC TESTS GO HERE
    # @unittest.skip('CTD-specific tests not yet defined')
    # def test_ctd(self):
    #    logger.info("A CTD-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()