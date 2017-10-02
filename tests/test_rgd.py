#!/usr/bin/env python3

import unittest
from dipper.sources.RGD import RGD


class RGDTestCase(unittest.TestCase):
    def setUp(self):
        self.test_set_1 = {'aspect': 'N',
                           'date': '2006-10-26',
                           'evidence': {'has_supporting_reference': ['RGD:1581841', 'PMID:12799311'],
                                        'type': 'IED',
                                        'with_support_from': []},
                           'negated': False,
                           'object': {'id': 'MP:0003340', 'taxon': 'NCBITaxon:10116'},
                           'provided_by': 'RGD',
                           'qualifiers': [],
                           'relation': {'id': None},
                           'source_line': 'RGD\t2535\tEdnra\t\tMP:0003340\tRGD:1581841|PMID:12799311\t'
                                          'IED\t\tN\tendothelin receptor type A\t\tgene\ttaxon:10116\t'
                                          '20061026\tRGD\t\t\n',
                           'subject': {'fullname': 'endothelin receptor type A',
                                       'id': 'RGD:2535',
                                       'label': 'Ednra',
                                       'synonyms': [],
                                       'taxon': {'id': 'NCBITaxon:10116'},
                                       'type': 'gene'},
                           'subject_extensions': [{'filler': '\n', 'property': 'isoform'}]}

        return

    def tearDown(self):
        return

    def testRGDParser(self):
        rgd = RGD('rdf_graph', True)
        rgd.make_association(record=self.test_set_1)
        sparql_query = """
        SELECT ?assoc WHERE {
        ?assoc a OBAN:association ;
        OBO:RO_0002558 OBO:ECO_0005611 ;
        dc:source <http://rgd.mcw.edu/rgdweb/report/reference/main.html?id=1581841> ;
        OBAN:association_has_object OBO:MP_0003340 ;
        OBAN:association_has_predicate OBO:RO_0002200 ;
        OBAN:association_has_subject <http://rgd.mcw.edu/rgdweb/report/gene/main.html?id=2535> ;
        pav:createdOn "2006-10-26" .
        <http://rgd.mcw.edu/rgdweb/report/gene/main.html?id=2535> OBO:RO_0002200 OBO:MP_0003340 .
        <http://rgd.mcw.edu/rgdweb/report/reference/main.html?id=1581841> a OBO:IAO_0000311 ;
                owl:sameAs <http://www.ncbi.nlm.nih.gov/pubmed/12799311> .
        }
        """
        #
        sparql_output = rgd.graph.query(sparql_query)
        self.assertEqual(len(list(sparql_output)), 1)





