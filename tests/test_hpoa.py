#!/usr/bin/env python3

from dipper.sources.HPOAnnotations import HPOAnnotations
from dipper import curie_map
from rdflib import Graph
from tests import test_general, test_source
from tests.test_source import SourceTestCase

import unittest
import logging
import os

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class HPOATestCase(SourceTestCase):
    def setUp(self):
        self.source = HPOAnnotations()
        self.source.test_ids = self._get_conf()['test_ids']['disease']
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