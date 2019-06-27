#!/usr/bin/env python3

import unittest
import logging
from rdflib import URIRef
from dipper.graph.RDFGraph import RDFGraph
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

logger = logging.getLogger(__name__)


class ModelTestCase(unittest.TestCase):

    def setUp(self):
        g = RDFGraph()
        self.model = Model(g)

        # stuff to make test triples
        self.test_cat_subj_curie = "MGI:1234"
        self.test_cat_subj = "http://www.informatics.jax.org/accession/MGI:1234"
        self.test_cat_default_pred = "http://w3id.org/biolink/vocab/category"
        self.test_named_indiv = "http://www.w3.org/2002/07/owl#NamedIndividual"
        self.test_label_pred = "http://www.w3.org/2000/01/rdf-schema#label"
        self.test_label = "some label"

    def tearDown(self):
        self.graph = None

    def test_addIndividualToGraph_assign_label(self):
        self.model.addIndividualToGraph(self.test_cat_subj_curie,
                                        "some label")

        label_triple = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj), URIRef(self.test_label_pred), None)))
        self.assertEqual(len(label_triple), 1, "method didn't assign label")
        self.assertEqual(str(label_triple[0][2]), self.test_label,
                         "method didn't assign correct label")

    def test_addIndividualToGraph_assign_type_named_individual(self):
        self.model.addIndividualToGraph(self.test_cat_subj_curie,
                                        "some label")
        triples = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj), None, URIRef(self.test_named_indiv))))
        self.assertEqual(len(triples), 1,
                         "method didn't assign type as named individual")

    def test_addIndividualToGraph_assign_category(self):
        self.model.addIndividualToGraph(self.test_cat_subj_curie,
                                        "some label",
                                        ind_category=blv.genotype.value)

        triples = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj), URIRef(self.test_cat_default_pred), None)))
        self.assertEqual(len(triples), 1, "method didn't assign category")

