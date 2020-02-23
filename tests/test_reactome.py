#!/usr/bin/env python3

import unittest
from dipper.sources.Reactome import Reactome
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph


class ReactomeTestCase(unittest.TestCase):

    def setUp(self):
        self.test_util = TestUtils()
        self.test_set_1 = \
            ('ENSBTAP00000013354',
             'R-BTA-3000480',
             'http://www.reactome.org/PathwayBrowser/#/R-BTA-3000480',
             'Scavenging by Class A Receptors',
             'IEA',
             'Bos taurus')
        self.gaf_eco = {"IEA": "ECO:0000501"}
        return

    def tearDown(self):
        return

    def testEnsemblReactomeParser(self):
        '''

        '''
        reactome = Reactome('rdf_graph', True)
        reactome.graph = RDFGraph(True)
        self.assertTrue(len(list(reactome.graph)) == 0)
        # reactome.parse_gaf_eco('gaf-eco-mapping')

        (gene,
         pathway_id,
         pathway_iri,
         pathway_label,
         go_ecode, species_name) = self.test_set_1
        reactome._add_component_pathway_association(
            'ENSEMBL:' + gene, 'REACT:' + pathway_id, pathway_label,
            self.gaf_eco[go_ecode])

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
        """
        self.assertTrue(self.test_util.test_graph_equality(triples, reactome.graph))
