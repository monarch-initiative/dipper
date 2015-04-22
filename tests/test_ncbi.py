#!/usr/bin/env python3

from dipper.sources.NCBIGene import NCBIGene
from dipper import curie_map
from rdflib import Graph
from tests import test_general, test_source
from tests.test_source import SourceTestCase

import unittest
import logging
import os

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class NCBITestCase(SourceTestCase):

    def setUp(self):
        self.source = NCBIGene()
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are hitting all parts of the code
    #@unittest.skip('test not yet defined')
    #def test_ncbitest(self):
    #    logger.info("An NCBI-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()