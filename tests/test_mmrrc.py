#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.MMRRC import MMRRC


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class MMRRCTestCase(SourceTestCase):
    def setUp(self):
        self.source = MMRRC()
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


if __name__ == '__main__':
    unittest.main()
