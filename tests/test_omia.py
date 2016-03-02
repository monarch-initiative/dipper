#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.OMIA import OMIA


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class OMIATestCase(SourceTestCase):

    def setUp(self):
        self.source = OMIA()
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()
