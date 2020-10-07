#!/usr/bin/env python3

import unittest
import logging
from rdflib import URIRef, Literal
from dipper import curie_map
from dipper.graph.RDFGraph import RDFGraph
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv
from dipper.utils.CurieUtil import CurieUtil

logger = logging.getLogger(__name__)


class ModelTestCase(unittest.TestCase):

    def setUp(self):
        g = RDFGraph()
        self.model = Model(g)

        this_curie_map = curie_map.get()
        self.cutil = CurieUtil(this_curie_map)

        # stuff to make test triples
        self.test_cat_subj_curie = "MGI:1234"
        self.test_cat_subj = self.cutil.get_uri("MGI:1234")
        self.test_cat_default_pred = self.cutil.get_uri("biolink:category")
        self.test_named_indiv = self.cutil.get_uri("owl:NamedIndividual")
        self.test_label_pred = self.cutil.get_uri("rdfs:label")
        self.test_label = "some label"

        self.test_comment_IRI = self.cutil.get_uri("rdfs:comment")
        self.test_comment = 'bonus eruptus'

    def tearDown(self):
        self.graph = None

    def test_addIndividualToGraph_assign_label(self):
        self.model.addIndividualToGraph(self.test_cat_subj_curie,
                                        "some label")

        label_triple = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj),
             URIRef(self.test_label_pred),
             None)))

        self.assertEqual(len(label_triple), 1, "method didn't assign label")
        self.assertEqual(str(label_triple[0][2]), self.test_label,
                         "method didn't assign correct label")

    def test_addIndividualToGraph_assign_type_named_individual(self):
        self.model.addIndividualToGraph(self.test_cat_subj_curie,
                                        "some label")

        triples = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj),
             None,
             URIRef(self.test_named_indiv))))

        self.assertEqual(len(triples), 1,
                         "method didn't assign type as named individual")

    def test_addIndividualToGraph_assign_category(self):
        self.model.addIndividualToGraph(self.test_cat_subj_curie,
                                        "some label",
                                        ind_category=blv.terms['Genotype'])

        triples = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj),
             URIRef(self.test_cat_default_pred),
             None)))

        self.assertEqual(len(triples), 1, "method didn't assign category")

    def test_add_comment(self):
        self.model.addComment(self.test_cat_subj, self.test_comment)

        triples = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj),
             URIRef(self.test_comment_IRI),
             Literal(self.test_comment))))

        self.assertEqual(len(triples), 1, "method didn't assign comment")

    def test_add_comment_assign_subject_category(self):
        self.model.addComment(self.test_cat_subj, self.test_comment,
                              subject_category=blv.terms['Genotype'])

        triples = list(self.model.graph.triples(
            (URIRef(self.test_cat_subj),
             URIRef(self.test_cat_default_pred),
             None)))
        self.assertEqual(len(triples), 1, "method didn't assign category")
