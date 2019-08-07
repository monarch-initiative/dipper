#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.Source import Source
from dipper.graph.RDFGraph import RDFGraph
from dipper import curie_map as curiemap
from rdflib import URIRef, Literal
from inspect import getdoc

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class DatasetTestCase(unittest.TestCase):
    """
    For testing requirements of the Dataset description

    Dataset creates a graph describing the metadata associated with the dataset in
    question, which should follow the HCLS specification for dataset descriptions
    https://www.w3.org/TR/2015/NOTE-hcls-dataset-20150514/

    Triples in the Dataset graph that describe the metadata are created as side effects
    in Source.py (mostly in the constructor and in get_files()). To test, we make an
    instance of a test ingest class (below) that creates the Dataset graph by
    interacting with Dataset and Source objects in the same way as other source classes
    (CTD, RGD, OMIM, etc.), then test the resulting Dataset graph contained in the test
    ingest object
    """

    def setUp(self):
        self.curie_map = curiemap.get()
        self.identifier = "fakeingest"
        self.ingest_url = "http://fakeingest.com"
        self.ingest_title = "this ingest title"

        # load source and fetch files to make dataset graph containing metadata
        self.source = FakeIngestClass("rdf_graph",
                                      are_bnodes_skolemized=False,
                                      identifier=self.identifier,
                                      ingest_url=self.ingest_url,
                                      ingest_title=self.ingest_title)
        self.source.fetch()

        # expected summary level IRI
        self.summary_level_IRI = URIRef(self.curie_map.get("MONARCH") + ":" +
                                        self.identifier)

    def tearDown(self):
        pass

    def test_has_dataset_attribute(self):
        self.assertTrue(hasattr(self.source, "dataset"),
                        "source object doesn't have dataset attribute")

    def test_dataset_has_graph(self):
        self.assertIsInstance(self.source.graph, RDFGraph,
                              "dataset doesn't contain an RDF graph")

    def test_summary_level_type(self):
        triples = list(self.source.dataset.graph.triples(
            (URIRef(self.summary_level_IRI),
             URIRef(self.curie_map.get("rdf") + "type"),
             URIRef(self.curie_map.get("dctypes") + "Dataset"))))
        self.assertTrue(len(triples) == 1, "missing summary level type triple")

    def test_summary_level_title(self):
        triples = list(self.source.dataset.graph.triples(
            (URIRef(self.summary_level_IRI),
             URIRef(self.curie_map.get("dcterms") + "title"),
             Literal(self.ingest_title))))
        self.assertTrue(len(triples) == 1, "missing summary level title triple")

    def test_summary_level_description(self):
        all_triples = list(self.source.dataset.graph.triples((None, None, None)))
        triples = list(self.source.dataset.graph.triples(
            (URIRef(self.summary_level_IRI),
             URIRef(self.curie_map.get("dc") + "description"),
             Literal("Fake ingest to test metadata in Dataset graph"))))
        self.assertTrue(len(triples) == 1, "missing summary level description")


    def test_summary_level_web_page(self):
        triples = list(self.source.dataset.graph.triples(
            (URIRef(self.summary_level_IRI),
             URIRef(self.curie_map.get("foaf") + "page"),
             URIRef(self.ingest_url))))
        self.assertTrue(len(triples) == 1, "missing summary level ingest URL triple")


if __name__ == '__main__':
    unittest.main()


class FakeIngestClass(Source):
    """
    Fake ingest to test metadata in Dataset graph
    """
    BASE_URL = 'https://data.monarchinitiative.org/'
    # using robots.txt b/c it's a trivially small file/empty file that we control
    files = {
        'test_file': {
            'file': 'test_file.txt',
            'url': BASE_URL + 'robots.txt'},
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 identifier=None,
                 ingest_url=None,
                 ingest_title=None,
                 ingest_desc=None
                 ):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            name=identifier,
            ingest_url=ingest_url,
            ingest_title=ingest_title
        )

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)
        return

    def parse(self):
        pass
