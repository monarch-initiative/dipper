#!/usr/bin/env python3

import unittest
import logging
# import os
# from rdflib import Graph
# from tests import test_general, test_source
from tests.test_source import SourceTestCase
from dipper.sources.IMPC import IMPC
# from dipper import curie_map

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class IMPCTestCase(SourceTestCase):

    def setUp(self):
        self.source = IMPC()
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # @unittest.skip('test not yet defined')
    # def test_hpotest(self):
    #    logger.info("An IMPC-specific test")
    #
    #    return


class EvidenceProvenanceTestCase():

    def setUp(self):
        return

    def tearDown(self):
        return

    # @unittest.skip('test not yet defined')
    # def test_hpotest(self):
    #    logger.info("An IMPC-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()
