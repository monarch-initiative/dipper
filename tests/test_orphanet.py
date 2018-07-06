#!/usr/bin/env python3

import unittest
import logging
import os
from tests.test_source import SourceTestCase
from dipper.sources.Orphanet import Orphanet
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph


logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


class OrphanetTestCase(SourceTestCase):

    def setUp(self):
        self.source = Orphanet('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


class GeneVariantDiseaseTest(unittest.TestCase):

    def setUp(self):
        """
        """
        self.test_util = TestUtils()
        self.orphanet = Orphanet('rdf_graph', True)
        # Override so tests don't break when we update terms
        self.global_terms = self.orphanet.open_and_parse_yaml(
            os.path.join(os.path.dirname(__file__), './resources/global_terms.yaml'))
        self.orphanet.rawdir = os.path.join(os.path.dirname(__file__), 'resources/orphanet')

    def tearDown(self):
        self.orphanet = None
        return

    def test_germline_variant_to_disease(self):
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-germline.xml'

        self.orphanet._process_diseasegene(limit=None)
        logger.debug("Reference graph: %s",
                     self.orphanet.graph.serialize(format="turtle")
                                   .decode("utf-8")
        )
        expected_triples = """
:MONARCH_a37b628d8347ddb0 a OBAN:association ;
    OBO:RO_0002558 OBO:ECO_0000322 ;
    OBAN:association_has_object <http://www.orpha.net/ORDO/Orphanet_938475> ;
    OBAN:association_has_predicate OBO:RO_0002200 ;
    OBAN:association_has_subject <https://monarchinitiative.org/.well-known/genid/b56f798350412a34> .

ENSEMBL:ENSG00000166813 a owl:Class .

<http://identifiers.org/hgnc/HGNC:30497> a owl:Class .

<http://www.orpha.net/ORDO/Orphanet_268061> a owl:Class ;
    rdfs:label "KS1" ;
    dc:description "kinesin family member 7" ;
    OIO:hasExactSynonym "KAS1" ;
    rdfs:subClassOf OBO:SO_0001217 ;
    owl:equivalentClass ENSEMBL:ENSG00000166813,
        <http://identifiers.org/hgnc/HGNC:30497> .

<https://monarchinitiative.org/.well-known/genid/b56f798350412a34> a OBO:GENO_0000002 ;
    rdfs:label "germline variant of KS1" ;
    OBO:GENO_0000418 <http://www.orpha.net/ORDO/Orphanet_268061> ;
    OBO:RO_0002200 <http://www.orpha.net/ORDO/Orphanet_938475> ;
    :MONARCH_anonymous true ;
    :has_cell_origin OBO:GENO_0000900 .

<http://www.orpha.net/ORDO/Orphanet_938475> a owl:Class ;
    rdfs:label "too much unit testing disorder" .
        """
        self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, self.orphanet.graph))

    def test_germline_lof_variant_to_disease(self):
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-germline-lof.xml'

        self.orphanet._process_diseasegene(limit=None)
        logger.debug("Reference graph: %s",
                     self.orphanet.graph.serialize(format="turtle")
                                  .decode("utf-8")
        )
        expected_triples = """
:MONARCH_f3a74695c7999465 a OBAN:association ;
    OBO:RO_0002558 OBO:ECO_0000322 ;
    OBAN:association_has_object <http://www.orpha.net/ORDO/Orphanet_938475> ;
    OBAN:association_has_predicate OBO:RO_0002200 ;
    OBAN:association_has_subject <https://monarchinitiative.org/.well-known/genid/ba0884fb61004110> .

ENSEMBL:ENSG00000166813 a owl:Class .

<http://identifiers.org/hgnc/HGNC:30497> a owl:Class .

<http://www.orpha.net/ORDO/Orphanet_268061> a owl:Class ;
    rdfs:label "KS1" ;
    dc:description "kinesin family member 7" ;
    OIO:hasExactSynonym "KAS1" ;
    rdfs:subClassOf OBO:SO_0001217 ;
    owl:equivalentClass ENSEMBL:ENSG00000166813,
        <http://identifiers.org/hgnc/HGNC:30497> .

<https://monarchinitiative.org/.well-known/genid/ba0884fb61004110> a OBO:GENO_0000002 ;
    rdfs:label "germline loss of function variant of KS1" ;
    OBO:GENO_0000418 <http://www.orpha.net/ORDO/Orphanet_268061> ;
    OBO:RO_0002200 <http://www.orpha.net/ORDO/Orphanet_938475> ;
    :MONARCH_anonymous true ;
    :has_cell_origin OBO:GENO_0000900 ;
    :has_functional_consequence OBO:SO_0002054 .

<http://www.orpha.net/ORDO/Orphanet_938475> a owl:Class ;
    rdfs:label "too much unit testing disorder" .
        """
        self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, self.orphanet.graph))

    def test_gene_to_disease(self):
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-no-variant.xml'

        self.orphanet._process_diseasegene(limit=None)
        logger.debug("Reference graph: %s",
                     self.orphanet.graph.serialize(format="turtle")
                                  .decode("utf-8")
        )
        expected_triples = """
<https://monarchinitiative.org/MONARCH_2ee6feba4c6269aa> a OBAN:association ;
    OBO:RO_0002558 OBO:ECO_0000322 ;
    OBAN:association_has_object <http://www.orpha.net/ORDO/Orphanet_938475> ;
    OBAN:association_has_predicate OBO:RO_0002326 ;
    OBAN:association_has_subject <http://www.orpha.net/ORDO/Orphanet_268061> .

ENSEMBL:ENSG00000166813 a owl:Class .

<http://identifiers.org/hgnc/HGNC:30497> a owl:Class .

<http://www.orpha.net/ORDO/Orphanet_268061> a owl:Class ;
    rdfs:label "KS1" ;
    OBO:RO_0002326 <http://www.orpha.net/ORDO/Orphanet_938475> ;
    dc:description "kinesin family member 7" ;
    OIO:hasExactSynonym "KAS1" ;
    rdfs:subClassOf OBO:SO_0001217 ;
    owl:equivalentClass ENSEMBL:ENSG00000166813,
        <http://identifiers.org/hgnc/HGNC:30497> .

<http://www.orpha.net/ORDO/Orphanet_938475> a owl:Class ;
    rdfs:label "too much unit testing disorder" .
        """
        self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, self.orphanet.graph))

    def test_unmapped_disease_assoc_type(self):
        """
        Test that a gene disease type that we have
        not mapped in translationtable/orphanet.yaml
        raises a ValueError
        """
        self.orphanet.graph = RDFGraph()  # Reset graph
        self.orphanet.files['disease-gene']['file'] = 'orph-no-mapping.xml'
        self.assertRaises(ValueError, lambda: self.orphanet._process_diseasegene(limit=None))


if __name__ == '__main__':
    unittest.main()
