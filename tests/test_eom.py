#!/usr/bin/env python3

from dipper.sources.EOM import EOM
from tests.test_source import SourceTestCase

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class EOMTestCase(SourceTestCase):

    def setUp(self):
        self.source = EOM()
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are hitting all parts of the code
    #@unittest.skip('test not yet defined')
    #def test_eomtest(self):
    #    logger.info("A EOM-specific test")
    #
    #    return

if __name__ == '__main__':
    unittest.main()