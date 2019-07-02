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
        self.test_cat_genotype_category = self.cutil.get_uri(blv.Genotype.value)
        self.test_cat_background_category = self.cutil.get_uri(
            blv.PopulationOfIndividualOrganisms.value)

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

    def test_addGenomicBackgroundToGenotype_adds_genotype(self):
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


    def test_addGenomicBackgroundToGenotype_adds_categories(self):
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

    def test_addParts(self):
        """
        """
        if part_relationship is None:
            part_relationship = self.globaltt['has_part']
        # Fail loudly if parent or child identifiers are None
        if parent_id is None:
            raise TypeError('Attempt to pass None as parent')
        elif part_id is None:
            raise TypeError('Attempt to pass None as child')
        elif part_relationship is None:
            part_relationship = self.globaltt['has_part']

        self.graph.addTriple(parent_id, part_relationship, part_id,
                             subject_category=subject_category,
                             object_category=object_category)

        return

if __name__ == '__main__':
    unittest.main()
