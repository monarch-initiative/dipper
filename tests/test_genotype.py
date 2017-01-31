#!/usr/bin/env python3

import unittest
import logging
from dipper.models.Genotype import Genotype
from dipper.graph.RDFGraph import RDFGraph
from dipper import curie_map

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class GenotypeTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = RDFGraph()
        self.curie_map = curie_map.get()
        self.genotype = Genotype(self.graph)

    def tearDown(self):
        self.genotype = None

    def test_addGenotype(self):
        from rdflib.namespace import RDFS, URIRef
        from rdflib import Literal
        from dipper.utils.CurieUtil import CurieUtil
        cutil = CurieUtil(self.curie_map)
        gid = 'MGI:5515892'
        label = \
            'Pmp22<Tr-2J>/Pmp22<+> [C57BL/6J-Pmp22<Tr-2J>/GrsrJ]'
        self.genotype.addGenotype(gid, label)
        self.assertTrue(
            (URIRef(cutil.get_uri(gid)), RDFS['label'],
             Literal(label)) in self.genotype.graph)


if __name__ == '__main__':
    unittest.main()
