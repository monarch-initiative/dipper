#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.RGD import RGD
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.TestUtils import TestUtils

logging.basicConfig()
logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


class RGDTestCase(unittest.TestCase):
    def setUp(self):
        self.test_util = TestUtils()
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
        rgd.graph = RDFGraph(True)

        self.assertTrue(len(list(rgd.graph)) == 0)

        rgd.make_association(record=self.test_set_1)
        triples = """
    :MONARCH_b4650e8c3d865f11a1a5 a OBAN:association ;
        RO:0002558 ECO:0005611 ;
        dc:source RGDRef:1581841 ;
        OBAN:association_has_object OBO:MP_0003340 ;
        OBAN:association_has_predicate OBO:RO_0002200 ;
        OBAN:association_has_subject RGD:2535 ;
        pav:createdOn "2006-10-26" .
    
    RGD:2535 OBO:RO_0002200 MP:0003340 .
        RGDRef:1581841 a IAO:0000311 ;
        owl:sameAs PMID:12799311 .
        """
        # dbg
        logger.debug("Reference graph: %s",
                     rgd.graph.serialize(format="turtle")
                              .decode("utf-8")
        )
        self.assertTrue(self.test_util.test_graph_equality(
            triples, rgd.graph))
