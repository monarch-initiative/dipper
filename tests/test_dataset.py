#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.Source import Source
from dipper.graph.RDFGraph import RDFGraph
from dipper import curie_map as curiemap
from rdflib.namespace import *
from rdflib import URIRef

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
        self.ingest_id = "fakeingest"

        # load source and fetch files to make dataset graph containing metadata
        self.source = FakeIngestClass("rdf_graph", False, self.ingest_id)
        self.source.fetch()

        # expected summary level IRI
        self.summary_level_IRI = ":".join(URIRef(self.curie_map.get("MONARCH"),
                                                 self.ingest_id))

    def tearDown(self):
        pass

    def test_has_dataset_attribute(self):
        self.assertTrue(hasattr(self.source, "dataset"),
                        "source object doesn't have dataset attribute")

    def test_dataset_has_graph(self):
        self.assertIsInstance(self.source.graph, RDFGraph,
                              "dataset doesn't contain a graph")

    def test_summary_level_type(self):
        all_triples = list(self.source.dataset.graph.triples((None, None, None)))
        summary_level_type_triple = list(self.source.dataset.graph.triples(
            (URIRef(self.summary_level_IRI), RDF.type,
             URIRef("http://purl.org/dc/dcmitype/Dataset"))))

        self.assertTrue(len(summary_level_type_triple) == 1,
                        "missing summary level type triples")


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

    def __init__(self, graph_type, are_bnodes_skolemized, identifier):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            identifier
        )

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)
        return

    def parse(self):
        pass
