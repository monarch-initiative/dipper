#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.Panther import Panther

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class PantherTestCase(SourceTestCase):

    def setUp(self):
        self.source = Panther('rdf_graph', True)
        self.source.test_ids = self._get_conf()['test_ids']['protein']
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are
    # hitting all parts of the code
    # @unittest.skip('test not yet defined')
    # def test_panthertest(self):
    #    logger.info("A Panther-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()
