#!/usr/bin/env python3

from dipper.sources.HGNC import HGNC
from tests.test_source import SourceTestCase

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class HGNCTestCase(SourceTestCase):

    def setUp(self):
        self.source = HGNC()
        self.source.test_ids = self._get_conf()['test_ids']['gene']
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()