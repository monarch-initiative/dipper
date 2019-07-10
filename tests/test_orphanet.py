#!/usr/bin/env python3

import unittest
import logging
import os
# from tests.test_source import SourceTestCase
from dipper.sources.Orphanet import Orphanet
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph

logging.getLogger().setLevel(logging.DEBUG)
LOG = logging.getLogger(__name__)


class GeneVariantDiseaseTest(unittest.TestCase):

    def setUp(self):
        """
        """
        self.test_util = TestUtils()
        self.orphanet = Orphanet('rdf_graph', True)
        self.orphanet.rawdir = os.path.join(
            os.path.dirname(__file__), 'resources/orphanet')

    def tearDown(self):
        self.orphanet = None
        return

    def test_germline_variant_to_disease(self):
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-germline.xml'

        self.orphanet._process_diseasegene(limit=None)
        LOG.debug(
            "Reference graph: %s",
            self.orphanet.graph.serialize(format="turtle").decode("utf-8")
        )
        expected_triples = """
MONARCH:ba2ac5d2153c70e2bb98 a OBAN:association ;
    RO:0002558 ECO:0000322 ;
    OBAN:association_has_object ORPHA:938475 ;
    OBAN:association_has_predicate RO:0004013 ;
    OBAN:association_has_subject HGNC:30497 .

ENSEMBL:ENSG00000166813 a owl:Class .

HGNC:30497 a owl:Class ;
    RO:0004013 ORPHA:938475 ;
    oboInOwl:hasExactSynonym "KAS1" ;
    owl:equivalentClass ENSEMBL:ENSG00000166813,
       ORPHA:268061 .

ORPHA:268061 a owl:Class .

ORPHA:938475 a owl:Class ;
    rdfs:label "too much unit testing disorder" .
    
ECO:0000322 biolink:category biolink:EvidenceType .
HGNC:30497 biolink:category biolink:Genotype .
ORPHA:938475 biolink:category biolink:PhenotypicFeature .

        """
        self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, self.orphanet.graph))
        return

    def test_germline_lof_variant_to_disease(self):
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-germline-lof.xml'

        self.orphanet._process_diseasegene(limit=None)
        LOG.debug(
            "Reference graph: %s",
            self.orphanet.graph.serialize(format="turtle").decode("utf-8")
        )
        expected_triples = """
MONARCH:b9ad1b0c562ad4db3f1e a OBAN:association ;
    RO:0002558 ECO:0000322 ;
    OBAN:association_has_object ORPHA:938475 ;
    OBAN:association_has_predicate RO:0004012 ;
    OBAN:association_has_subject ORPHA:268061 .

ORPHA:268061 RO:0004012 ORPHA:938475 ;
    oboInOwl:hasExactSynonym "KAS1" .

ORPHA:938475 a owl:Class ;
    rdfs:label "too much unit testing disorder" .
    
ECO:0000322 biolink:category biolink:EvidenceType .
ORPHA:268061 biolink:category biolink:Genotype .
ORPHA:938475 biolink:category biolink:PhenotypicFeature .
    
        """
        self.assertTrue(
            self.test_util.test_graph_equality(expected_triples, self.orphanet.graph))
        return

    def test_gene_to_disease(self):
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-no-variant.xml'

        self.orphanet._process_diseasegene(limit=None)
        LOG.debug(
            "Reference graph: %s",
            self.orphanet.graph.serialize(format="turtle") .decode("utf-8")
        )
        expected_triples = """
MONARCH:bdbeb077e365ddedda20 a OBAN:association ;
    RO:0002558 ECO:0000322 ;
    OBAN:association_has_object ORPHA:938475 ;
    OBAN:association_has_predicate RO:0004015 ;
    OBAN:association_has_subject ORPHA:268061 .

ORPHA:268061 RO:0004015 ORPHA:938475 ;
    oboInOwl:hasExactSynonym "KAS1" .

ORPHA:938475 a owl:Class ;
    rdfs:label "too much unit testing disorder" .
    
ECO:0000322 biolink:category biolink:EvidenceType .
ORPHA:268061 biolink:category biolink:Genotype .
ORPHA:938475 biolink:category biolink:PhenotypicFeature .

        """
        self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, self.orphanet.graph))
        return

    def test_unmapped_disease_assoc_type(self):
        """
        Test that a gene disease type that we have
        not mapped in translationtable/orphanet.yaml
        raises a ValueError
        """
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-no-mapping.xml'
        self.assertRaises(
            KeyError, lambda: self.orphanet._process_diseasegene(limit=None))
        return


if __name__ == '__main__':
    unittest.main()
