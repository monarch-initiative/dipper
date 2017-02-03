#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.GeneOntology import GeneOntology
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class GeneOntologyTestCase(SourceTestCase):

    def setUp(self):
        self.source = GeneOntology('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()
