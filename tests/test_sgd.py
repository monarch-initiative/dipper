#!/usr/bin/env python3

import unittest
from dipper.sources.SGD import SGD
from dipper.utils.TestUtils import TestUtils
from dipper.graph.RDFGraph import RDFGraph


class SGDTestCase(unittest.TestCase):
    def setUp(self):
        self.test_util = TestUtils()
        self.test_set_1 = {'Allele': 'atp6-L183R (L183R)',
                           'Chemical': 'glycerol',
                           'Condition': 'elevated temperature (35 deg C)|nonfermentable carbon source',
                           'Details': 'similar results obtained with atp6-L247R, and atp6-W136R, all '
                                      'corresponding to human NARP syndrome mutants',
                           'Experiment Type': 'classical genetics',
                           'Feature Name': 'Q0085',
                           'Feature Type': 'ORF',
                           'Gene Name': 'ATP6',
                           'Mutant Type': 'reduction of function',
                           'Phenotype': 'respiratory growth: decreased rate',
                           'Reference': 'PMID: 21715656|SGD_REF: S000145858',
                           'Reporter': ' ',
                           'SGDID': 'S000007268',
                           'Strain Background': 'Other'}

        return

    def tearDown(self):
        return

    def testSGDParser(self):
        sgd = SGD('rdf_graph', True)
        sgd.graph = RDFGraph(True)
        record = self.test_set_1
        sgd.make_association(record)

        description = [
            'genomic_background: {}'.format(record['Strain Background']),
            'allele: {}'.format(record['Allele']),
            'chemical: {}'.format(record['Chemical']),
            'condition: {}'.format(record['Condition']),
            'details: {}'.format(record['Details']),
            'feature_name: {}'.format(record['Feature Name']),
            'gene_name: {}'.format(record['Gene Name']),
            'mutant_type: {}'.format(record['Mutant Type']),
            'reporter: {}'.format(record['Reporter']),
        ]
        description = " | ".join(description).strip()

        triples = """
        :MONARCH_95158d413dd73476 a OBAN:association ;
            OBO:RO_0002558 OBO:APO_0000020 ;
            dc:description "{0}";
            dc:source PMID:21715656 ;
            OBAN:association_has_object MONARCH:OBO_APO_0000309OBO_APO_0000245 ;
            OBAN:association_has_predicate OBO:RO_0002200 ;
            OBAN:association_has_subject SGD:S000007268 .
            
        SGD:S000007268 rdfs:label "ATP6" ;
        RO:0002200 MONARCH:OBO_APO_0000309OBO_APO_0000245 .

        APO:0000020 rdfs:label "classical genetics" .

        PMID:21715656 a OBO:IAO_0000311 ;
        owl:sameAs SGD_REF:S000145858 .

        MONARCH:OBO_APO_0000309OBO_APO_0000245 rdfs:label "respiratory growth:decreased rate" ;
        rdfs:subClassOf UPHENO:0001001 .

        """.format(description)
        # test exact contents of graph
        self.assertTrue(self.test_util.test_graph_equality(
            triples, sgd.graph))

