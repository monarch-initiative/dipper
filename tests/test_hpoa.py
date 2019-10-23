#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.HPOAnnotations import HPOAnnotations

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class HPOATestCase(SourceTestCase):
    def setUp(self):
        self.source = HPOAnnotations('rdf_graph', True)
        self.source.test_ids = self.all_test_ids['disease']
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    @unittest.skip('test not yet defined')
    def test_hpotest(self):
        logger.info("A HPO-specific test")

        return

        # TODO test if we know all of the aspect codes


if __name__ == '__main__':
    unittest.main()
