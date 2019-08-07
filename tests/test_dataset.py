#!/usr/bin/env python3

import unittest
import logging
from dipper.models.Dataset import Dataset
from dipper.sources.Source import Source

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
        self.source = FakeIngestClass("rdf_graph", False)

    def tearDown(self):
        pass

    def test_has_dataset_attribute(self):
        self.assertTrue(hasattr(self.source, "dataset"),
                        "source object doesn't have dataset attribute")

if __name__ == '__main__':
    unittest.main()

class FakeIngestClass(Source):
    """
    Fake ingest to test metadata in Dataset graph
    """
    BASE_URL = 'ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/'
    files = {
        'rat_gene2mammalian_phenotype': {
            'file': 'rattus_genes_mp',
            'url': BASE_URL + 'rattus_genes_mp'},
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'rgd',
            ingest_title='Fake ingest class',
            ingest_url='http://rgd.mcw.edu/',
            license_url=None,
            data_rights='https://rgd.mcw.edu/wg/disclaimer/',
        )
        self.dataset.set_citation('https://rgd.mcw.edu/wg/citing-rgd/')

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)
        return
