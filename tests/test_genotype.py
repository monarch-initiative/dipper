#!/usr/bin/env python3

from dipper.models.Genotype import Genotype
from dipper import curie_map
from rdflib import Graph

import unittest
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class GenotypeTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = Graph()
        self.curie_map = curie_map.get()
        self.genotype = Genotype(self.graph)

    def tearDown(self):
        self.genotype = None

    def test_addGenotype(self):
        from rdflib.namespace import RDFS,URIRef
        from rdflib import Literal
        from dipper.utils.CurieUtil import CurieUtil
        cu = CurieUtil(self.curie_map)
        id = 'MGI:5515892'
        label = \
            'Pmp22<Tr-2J>/Pmp22<+> [C57BL/6J-Pmp22<Tr-2J>/GrsrJ]'
        self.genotype.addGenotype(id, label)
        self.assertTrue((URIRef(cu.get_uri(id)), RDFS['label'],
                         Literal(label)) in self.genotype.graph)


if __name__ == '__main__':
    unittest.main()
