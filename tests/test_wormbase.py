#!/usr/bin/env python3

from dipper.sources.WormBase import WormBase
from tests.test_source import SourceTestCase

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class WormBaseTestCase(SourceTestCase):
    def setUp(self):
        self.source = WormBase()
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()