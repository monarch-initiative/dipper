#!/usr/bin/env python3

import unittest
from dipper.sources.Reactome import Reactome
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph


class ReactomeTestCase(unittest.TestCase):

    def setUp(self):
        self.test_util = TestUtils()
        self.test_set_1 = \
            ('ENSBTAP00000013354', 'R-BTA-3000480',
             'http://www.reactome.org/PathwayBrowser/#/R-BTA-3000480',
             'Scavenging by Class A Receptors',	'IEA', 'Bos taurus')
        return

    def tearDown(self):
        return

    def testEnsemblReactomeParser(self):
        reactome = Reactome('rdf_graph', True)
        reactome.graph = RDFGraph(True)
        self.assertTrue(len(list(reactome.graph)) == 0)

        eco_map = Reactome.get_eco_map(Reactome.map_files['eco_map'])
        (gene, pathway_id, pathway_iri, pathway_label,
         go_ecode, species_name) = self.test_set_1
        reactome._add_component_pathway_association(
            eco_map, gene, 'ENSEMBL', pathway_id,
            'REACT', pathway_label, go_ecode)

        triples = """
        ENSEMBL:ENSBTAP00000013354 RO:0002331 REACT:R-BTA-3000480 .
        
        :MONARCH_b582c188b7ec20016206 a OBAN:association ;
            OBO:RO_0002558 ECO:0000501 ;
            OBAN:association_has_object REACT:R-BTA-3000480 ;
            OBAN:association_has_predicate RO:0002331 ;
            OBAN:association_has_subject ENSEMBL:ENSBTAP00000013354 .

        REACT:R-BTA-3000480 a owl:Class ;
            rdfs:label "Scavenging by Class A Receptors" ;
            rdfs:subClassOf GO:0009987,
                PW:0000001 .
                
        ECO:0000501 biolink:category biolink:EvidenceType .
        REACT:R-BTA-3000480 biolink:category biolink:Pathway .
        ENSEMBL:ENSBTAP00000013354 biolink:category biolink:Gene .
        """
        self.assertTrue(self.test_util.test_graph_equality(
            triples, reactome.graph))
