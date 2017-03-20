#!/usr/bin/env python3

import unittest
from dipper.sources.Reactome import Reactome


class ReactomeTestCase(unittest.TestCase):

    def setUp(self):
        self.test_set_1 = \
            ('ENSBTAP00000013354', 'R-BTA-3000480',
             'http://www.reactome.org/PathwayBrowser/#/R-BTA-3000480',
             'Scavenging by Class A Receptors',	'IEA', 'Bos taurus')
        return

    def tearDown(self):
        return

    def testEnsemblReactomeParser(self):
        reactome = Reactome('rdf_graph', True)
        eco_map = Reactome.get_eco_map(Reactome.map_files['eco_map'])
        (gene, pathway_id, pathway_iri, pathway_label,
         go_ecode, species_name) = self.test_set_1
        reactome._add_gene_pathway_association(eco_map, gene, pathway_id,
                                               pathway_label, go_ecode)

        sparql_query = """
                      SELECT ?assoc
                      WHERE {
                           ENSEMBL:ENSBTAP00000013354 OBO:RO_0002331 REACT:R-BTA-3000480 .

                           ?assoc a OBAN:association ;
                               OBO:RO_0002558 OBO:ECO_0000501 ;
                               OBAN:association_has_object REACT:R-BTA-3000480 ;
                               OBAN:association_has_predicate OBO:RO_0002331 ;
                               OBAN:association_has_subject ENSEMBL:ENSBTAP00000013354 .

                           REACT:R-BTA-3000480 a owl:Class ;
                               rdfs:label "Scavenging by Class A Receptors" ;
                               rdfs:subClassOf OBO:GO_0009987,
                                   OBO:PW_0000001 .

                      }
                      """
        sparql_output = reactome.graph.query(sparql_query)

        self.assertEqual(len(list(sparql_output)), 1)
