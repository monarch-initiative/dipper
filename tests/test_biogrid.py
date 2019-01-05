#!/usr/bin/env python3
import unittest
import logging
from dipper.sources.BioGrid import BioGrid
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
LOG = logging.getLogger(__name__)


class BioGridTestCase(SourceTestCase):

    def setUp(self):
        self.source = BioGrid('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # @unittest.skip('Biogrid-specific tests not yet defined')
    # def test_biogrid(self):
    #    LOG.info("A BioGrid-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()
