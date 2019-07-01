#!/usr/bin/env python3

import unittest
import logging
from rdflib.namespace import RDFS, URIRef
from rdflib import Literal
from dipper import curie_map
from dipper.models.Genotype import Genotype
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.CurieUtil import CurieUtil
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class GenotypeTestCase(unittest.TestCase):

    def setUp(self):
        self.graph = RDFGraph()
        self.curie_map = curie_map.get()
        self.genotype = Genotype(self.graph)
        self.cutil = CurieUtil(self.curie_map)
        self.test_cat_pred = self.cutil.get_uri(blv.category.value)
        self.test_cat_genotype_category = self.cutil.get_uri(blv.genotype.value)
        self.test_cat_background_category = self.cutil.get_uri(
            blv.population_of_individual_organisms.value)

    def tearDown(self):
        self.genotype = None

    def test_addGenotype(self):
        cutil = CurieUtil(self.curie_map)
        gid = 'MGI:5515892'
        label = \
            'Pmp22<Tr-2J>/Pmp22<+> [C57BL/6J-Pmp22<Tr-2J>/GrsrJ]'
        self.genotype.addGenotype(gid, label)
        self.assertTrue(
            (URIRef(cutil.get_uri(gid)), RDFS['label'],
             Literal(label)) in self.genotype.graph)

    def test_addGenomicBackgroundToGenotype(self):
        """
         test that addGenomicBackgroundToGenotype() correctly assigns
         subject/object category
         """
        genotype_id = "GENO:0000002"
        background_id = "GENO:0000002" # no idea what a good example background ID is
        self.genotype.addGenomicBackgroundToGenotype(
            background_id=background_id, genotype_id=genotype_id)

        geno_triples = list(self.graph.triples((
                            URIRef(self.cutil.get_uri(genotype_id)),
                            URIRef(self.test_cat_pred),
                            URIRef(self.test_cat_genotype_category))))
        self.assertEqual(len(geno_triples), 1,
                         "addTriples() didn't make exactly 1 genotype category triple")
        self.assertEqual(geno_triples[0][2], URIRef(self.test_cat_genotype_category),
                         "addTriples() didn't assign the right genotype category")

        background_triples = list(self.graph.triples((
                                            URIRef(self.cutil.get_uri(background_id)),
                                            URIRef(self.test_cat_pred),
                                            URIRef(self.test_cat_background_category))))
        self.assertEqual(len(background_triples), 1,
                         "addTriples() didn't make exactly 1 genotype category triple")
        self.assertEqual(background_triples[0][2],
                         URIRef(self.test_cat_background_category),
                         "addTriples() didn't assign the right background category")


if __name__ == '__main__':
    unittest.main()
