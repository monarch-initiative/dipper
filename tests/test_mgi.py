#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.MGI import MGI
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph
import os


logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


class MGITestCase(SourceTestCase):

    def setUp(self):
        self.source = MGI('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


class EvidenceTestCase(unittest.TestCase):

    def setUp(self):
        """
        Because _process_evidence_view uses
        self.rawdir to find the evidence file,
        the defaults are overriden here to
        point to our test file
        Note the file name must match what is in
        that method - evidence_view
        """
        self.test_util = TestUtils()
        self.mgi = MGI('rdf_graph', True)
        self.mgi.rawdir = os.path.join(os.path.dirname(__file__), 'resources/mgi')
        self.mgi.idhash['annot']['6901981'] = ':association'

    def tearDown(self):
        self.mgi = None
        return

    def test_sex_specificity_model(self):
        self.mgi.graph = RDFGraph(True)  # Reset graph
        self.mgi._process_evidence_view(limit=None)
        logger.debug("Reference graph: %s",
                     self.mgi.graph.serialize(format="turtle")
                                   .decode("utf-8")
        )
        expected_triples = """
        :association RO:0002558 ECO:0000006 ;
            dc:source J:74619 ;
            :has_sex_specificity PATO:0000384 .

        J:74619 a IAO:0000310 .
        """
        self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, self.mgi.graph))


if __name__ == '__main__':
    unittest.main()
