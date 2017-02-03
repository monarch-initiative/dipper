#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.AnimalQTLdb import AnimalQTLdb
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class AnimalQTLdbTestCase(SourceTestCase):

    def setUp(self):
        self.source = AnimalQTLdb('rdf_graph', True)
        self.source.gene_ids = self._get_conf()['test_ids']['gene']
        self.source.disease_ids = self._get_conf()['test_ids']['disease']
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

if __name__ == '__main__':
    unittest.main()
