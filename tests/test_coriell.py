#!/usr/bin/env python3
import unittest
import logging
from dipper.sources.Coriell import Coriell
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class CoriellTestCase(SourceTestCase):

    def setUp(self):
        self.source = Coriell('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # @unittest.skip('Coriell-specific tests not yet defined')
    # def test_coriell(self):
    #    logger.info("A Coriell-specific test")
    #
    #    return

if __name__ == '__main__':
    unittest.main()
