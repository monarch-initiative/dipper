#!/usr/bin/env python3

from dipper.sources.MGI import MGI
from dipper import curie_map
from rdflib import Graph
from tests import test_general, test_source
from tests.test_source import SourceTestCase

import unittest
import logging
import os

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class MGITestCase(SourceTestCase):

    def setUp(self):
        self.source = MGI()
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are hitting all parts of the code
    #@unittest.skip('test not yet defined')
    #def test_mgitest(self):
    #    logger.info("An MGI-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()