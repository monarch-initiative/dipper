#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.OMA import OMA

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class OMATestCase(SourceTestCase):

    def setUp(self):
        self.source = OMA('rdf_graph', True)
        self.source.test_ids = self._get_testids['protein']
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


if __name__ == '__main__':
    unittest.main()
