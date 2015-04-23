#!/usr/bin/env python3

from dipper.sources.UCSCBands import UCSCBands
from tests.test_source import SourceTestCase

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class UCSCBandsTestCase(SourceTestCase):

    def setUp(self):
        self.source = UCSCBands()
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are hitting all parts of the code
    #@unittest.skip('test not yet defined')
    #def test_ucscbandstest(self):
    #    logger.info("A UCSCBands-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()