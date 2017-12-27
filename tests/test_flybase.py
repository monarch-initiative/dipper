#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.FlyBase import FlyBase
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class FlyBaseTestCase(SourceTestCase):
    def setUp(self):
        self.source = FlyBase('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()
