#!/usr/bin/env python3

import unittest
import logging
# import os
# from rdflib import Graph
# from tests import test_general, test_source
from tests.test_source import SourceTestCase
from dipper.sources.Orphanet import Orphanet
# from dipper import curie_map

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class OrphanetTestCase(SourceTestCase):

    def setUp(self):
        self.source = Orphanet('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()
