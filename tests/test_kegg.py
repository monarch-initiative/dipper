#!/usr/bin/env python3

from dipper.sources.KEGG import KEGG
from dipper import curie_map
from rdflib import Graph
from tests import test_general, test_source
from tests.test_source import SourceTestCase

import unittest
import logging
import os

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class KEGGTestCase(SourceTestCase):

    def setUp(self):
        self.source = KEGG()
        self.source.disease_ids = self._get_conf()['test_ids']['disease']
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests to make sure we are hitting all parts of the code


if __name__ == '__main__':
    unittest.main()