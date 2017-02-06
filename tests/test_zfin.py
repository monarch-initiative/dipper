#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.ZFIN import ZFIN

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class ZFINTestCase(SourceTestCase):

    def setUp(self):
        self.source = ZFIN('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    @unittest.skip(
        'Will eventually write test to check if phenotype sextuples' +
        'are mapped to ZP ids')
    def test_allZPAvailable(self):
        """
        This test will identify if there are
        any missing ZP terms in the mapping file
        :return:

        """
        # TODO add this test to check if all phenotype sextuples
        # are mapped to ZP ids

        return

if __name__ == '__main__':
    unittest.main()
