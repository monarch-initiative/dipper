#!/usr/bin/env python3

from dipper.sources.ClinVar import ClinVar
from tests.test_source import SourceTestCase

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

class ClinVarTestCase(SourceTestCase):

    def setUp(self):
        self.source = ClinVar()
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # @unittest.skip('Clinvar-specific tests not yet defined')
    # def test_clinvar(self):
    #    logger.info("A ClinVar-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()