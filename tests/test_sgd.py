#!/usr/bin/env python3

import unittest
from dipper.sources.SGD import SGD

test_set = [
    {'Allele': 'atp6-L183R (L183R)',
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
     'Strain Background': 'Other'},
    {'Allele': 'atp6-S250P (mutation analogous to a human MT-ATP6 mutation found '
               'in multiple patients with LS, NARP, CMT or spinocerebellar ataxia)',
     'Chemical': 'oligomycin (0.25 ug/ml)',
     'Condition': ' ',
     'Details': ' ',
     'Experiment Type': 'classical genetics',
     'Feature Name': 'Q0085',
     'Feature Type': 'ORF',
     'Gene Name': 'ATP6',
     'Mutant Type': 'reduction of function',
     'Phenotype': 'resistance to chemicals: decreased',
     'Reference': 'PMID: 24316278|SGD_REF: S000156155',
     'Reporter': ' ',
     'SGDID': 'S000007268',
     'Strain Background': 'Other'},
    {'Allele': 'atp6-S250P (mutation analogous to a human MT-ATP6 mutation found '
               'in multiple patients with LS, NARP, CMT or spinocerebellar ataxia)',
     'Chemical': ' ',
     'Condition': ' ',
     'Details': ' ',
     'Experiment Type': 'classical genetics',
     'Feature Name': 'Q0085',
     'Feature Type': 'ORF',
     'Gene Name': 'ATP6',
     'Mutant Type': 'reduction of function',
     'Phenotype': 'respiratory growth: normal',
     'Reference': 'PMID: 24316278|SGD_REF: S000156155',
     'Reporter': ' ',
     'SGDID': 'S000007268',
     'Strain Background': 'Other'},
    {'Allele': 'atp6-L252P (mutation analogous to a human MT-ATP6 mutation found '
               'in a patient with severe Leigh syndrome)',
     'Chemical': ' ',
     'Condition': ' ',
     'Details': ' ',
     'Experiment Type': 'classical genetics',
     'Feature Name': 'Q0085',
     'Feature Type': 'ORF',
     'Gene Name': 'ATP6',
     'Mutant Type': 'reduction of function',
     'Phenotype': 'respiratory growth: absent',
     'Reference': 'PMID: 24316278|SGD_REF: S000156155',
     'Reporter': ' ',
     'SGDID': 'S000007268',
     'Strain Background': 'Other'},
    {'Allele': 'atp6-L183R (L183R)',
     'Chemical': ' ',
     'Condition': ' ',
     'Details': ' ',
     'Experiment Type': 'classical genetics',
     'Feature Name': 'Q0085',
     'Feature Type': 'ORF',
     'Gene Name': 'ATP6',
     'Mutant Type': 'unspecified',
     'Phenotype': 'respiratory growth: decreased rate',
     'Reference': 'PMID: 17855363|SGD_REF: S000124079',
     'Reporter': ' ',
     'SGDID': 'S000007268',
     'Strain Background': 'Other (MR6)'},
    {'Allele': 'atp6-L247R (L247R; equivalent to pathogenic mutation in human '
               'ATP6)',
     'Chemical': 'oligomycin (0.25 ug/ml)',
     'Condition': 'nonfermentable carbon source (2% glycerol)',
     'Details': ' ',
     'Experiment Type': 'classical genetics',
     'Feature Name': 'Q0085',
     'Feature Type': 'ORF',
     'Gene Name': 'ATP6',
     'Mutant Type': 'unspecified',
     'Phenotype': 'resistance to chemicals: decreased',
     'Reference': 'PMID: 20056103|SGD_REF: S000132686',
     'Reporter': ' ',
     'SGDID': 'S000007268',
     'Strain Background': 'Other'}
]

class SGDTestCase(unittest.TestCase):
    def setUp(self):
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
        sgd.make_association(record=self.test_set_1)
        sparql_query = """
        SELECT ?assoc WHERE {
        ?assoc a OBAN:association ;
            OBO:RO_0002558 OBO:APO_0000020 ;
            dc:description ?description;
            OBAN:association_has_object <https://monarchinitiative.org/MONARCH_APO_0000309APO_0000245> ;
            OBAN:association_has_predicate OBO:RO_0002200 ;
            OBAN:association_has_subject SGD:S000007268 .
        }
        """

        sparql_output = sgd.graph.query(sparql_query)
        self.assertEqual(len(list(sparql_output)), 1)

