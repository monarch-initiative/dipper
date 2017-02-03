#!/usr/bin/env python3

import unittest
import logging
# import os
# from rdflib import Graph
# from tests import test_general, test_source
from tests.test_source import SourceTestCase
from dipper.sources.MPD import MPD
# from dipper import curie_map

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class MPDTestCase(SourceTestCase):

    def setUp(self):
        self.source = MPD('rdf_graph', True)
        self.source.settestonly(True)
        return

    def tearDown(self):
        self.source = None
        return

    # @unittest.skip('test not yet defined')
    # def test_hpotest(self):
    #    logger.info("An MPD-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()
