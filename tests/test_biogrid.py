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

    def test_interactor_to_gene_curie(self):
        self.assertEqual('NCBIGene:3645446', self.source._interactor_to_gene_curie(
                                             'entrez gene/locuslink:3645446'))
        self.assertEqual('BIOGRID:4383875', self.source._interactor_to_gene_curie(
                                             'biogrid:4383875'))
        self.assertEqual(None,  self.source._interactor_to_gene_curie('NOTAGENEID'))


if __name__ == '__main__':
    unittest.main()
