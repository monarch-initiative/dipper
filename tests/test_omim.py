#!/usr/bin/env python3

import unittest
import logging
# import os
# from rdflib import Graph
# from tests import test_general, test_source
from tests.test_source import SourceTestCase
from dipper.sources.OMIM import OMIM
# from dipper import curie_map

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class OMIMTestCase(SourceTestCase):

    def setUp(self):
        self.source = OMIM('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are
    # hitting all parts of the code
    # @unittest.skip('test not yet defined')
    # def test_omimtest(self):
    #    logger.info("An OMIM-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()
