#!/usr/bin/env python3
import unittest
import logging
from dipper.sources.CTD import CTD
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.TestUtils import TestUtils

logging.basicConfig()
logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


class CTDTestCase(unittest.TestCase):
    def setUp(self):
        self.test_util = TestUtils()
        self.source = CTD('rdf_graph', True)
        self.source.graph = RDFGraph(True)
        self.test_row = [
            'Nicotine',
            'D009538',
            '',
            'TOBACCO ADDICTION, SUSCEPTIBILITY TO',
            'OMIM:188890',
            'therapeutic',
            '',
            '',
            '',
            '12345|56789'
        ]
        return

    def tearDown(self):
        self.source = None
        return

    def test_therapeutic_relationship(self):
        # test that graph is empty
        self.assertTrue(len(list(self.source.graph)) == 0)

        self.source._process_interactions(self.test_row)

        triples = """
            :MONARCH_b6c289df47cb72653f79 a OBAN:association ;
                RO:0002558 ECO:0000033 ;
                dc:source PMID:12345, PMID:56789 ;
                OBAN:association_has_object OMIM:188890 ;
                OBAN:association_has_predicate RO:0002606 ;
                OBAN:association_has_subject MESH:D009538 .
            
            MESH:D009538 a owl:Class ;
                rdfs:label "Nicotine" ;
                biolink:category biolink:ChemicalSubstance ;
                RO:0002606 OMIM:188890 .
                
            PMID:12345 a IAO:0000013 .

            PMID:56789 a IAO:0000013 .

            OMIM:188890 a owl:Class ;
                biolink:category biolink:DiseaseOrPhenotypicFeature .
        """
        # test exact contents of graph
        self.assertTrue(self.test_util.test_graph_equality(
            triples, self.source.graph))


if __name__ == '__main__':
    unittest.main()
