#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.WormBase import WormBase


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class WormBaseTestCase(SourceTestCase):
    def setUp(self):
        self.source = WormBase('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()
