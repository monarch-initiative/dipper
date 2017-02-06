#!/usr/bin/env python3
import unittest
import logging
from dipper.sources.ClinVar import ClinVar
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class ClinVarTestCase(SourceTestCase):

    def setUp(self):
        self.source = ClinVar('rdf_graph', True)
        self.source.gene_ids = self._get_conf()['test_ids']['gene']
        self.source.disease_ids = self._get_conf()['test_ids']['disease']
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
