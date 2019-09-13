import unittest
from unittest.mock import MagicMock
import logging
import os
from datetime import datetime
from shutil import copyfile
from stat import ST_CTIME
from rdflib import XSD
from rdflib import URIRef, Literal

from dipper.sources.Source import Source
from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper import curie_map as curiemap

logging.basicConfig(level=logging.WARNING)
LOG = logging.getLogger(__name__)


class SourceMetadataTestCase(unittest.TestCase):
    """
    Some metadata triples to be emitted in Dataset's graph are created in
    Source.get_files() and elsewhere - these are tested below.

    This is a separate class from SourceTestCase, because the latter is called from
    many implementations of Source, and the tests below don't need to be called for
    that purpose.
    """

    @classmethod
    def setUpClass(cls):

        # feed these to FakeIngestClass, then check in tests
        cls.curie_map = curiemap.get()
        cls.identifier = "someid"
        cls.ingest_url = "http://sourceofdata.com"
        cls.ingest_title = "our transform of some source"
        cls.ingest_logo_url = cls.ingest_url + "/logo.png"
        cls.license_url = "https://choosealicense.com/licenses/apache-2.0/"
        cls.data_rights = "https://www.gnu.org/licenses/gpl-2.0.html"

        cls.theseFiles = {
            'test_file': {
                'file': 'test_file.txt',
                'url': 'http://fakeingest.com/remote_file.txt',
                'path_to_mock_download_file':
                    'tests/resources/fakeingest/test_file.txt'}
        }

        cls.thesePgFiles = {
            'test_file': {
                'file': 'test_file_pg.txt',
                'url': 'http://fakeingest.com/remote_file.txt',
                'path_to_mock_download_file':
                    'tests/resources/fakeingest/test_file_pg.txt'}
        }

        # expected things:
        cls.expected_curie_prefix = "DatasetBase"
        cls.timestamp_date = datetime.today().strftime("%Y%m%d")

        # expected summary level IRI
        cls.summary_level_IRI = URIRef(cls.curie_map.get(cls.expected_curie_prefix)
                                       + cls.identifier)
        # expected version level IRI
        cls.version_level_IRI = URIRef(cls.summary_level_IRI + "_" + cls.timestamp_date)

        # expected distribution level IRI (for ttl resource)
        cls.distribution_level_IRI_ttl = URIRef(cls.version_level_IRI + ".ttl")

        # expected IRIs
        cls.iri_version = URIRef(cls.curie_map.get("pav") + "version")
        cls.iri_source = URIRef(cls.curie_map.get("dcterms") + "source")
        cls.iri_retrieved_on = URIRef(cls.curie_map.get("pav") + "retrievedOn")
        cls.iri_triples_count = URIRef(cls.curie_map.get("void") + "triples")

    def setUp(self):
        # load source and fetch files to make dataset graph containing metadata
        self.source = FakeIngestClass("rdf_graph",
                                      are_bnodes_skolemized=False,
                                      identifier=self.identifier,
                                      ingest_url=self.ingest_url,
                                      ingest_title=self.ingest_title,
                                      ingest_logo=self.ingest_logo_url,
                                      license_url=self.license_url,
                                      data_rights=self.data_rights,
                                      files=self.theseFiles)
        self.source.fetch()
        # don't write metadata out with main graph in tests, because this screws up
        # triples counts during testing
        self.source.write(write_metadata_in_main_graph=False)

        # same as above, but using postgres class
        # load source and fetch files to make dataset graph containing metadata
        self.pg_source = FakeIngestUsingPostgres("rdf_graph",
                                                 are_bnodes_skolemized=False,
                                                 identifier=self.identifier,
                                                 ingest_url=self.ingest_url,
                                                 ingest_title=self.ingest_title,
                                                 ingest_logo=self.ingest_logo_url,
                                                 license_url=self.license_url,
                                                 data_rights=self.data_rights,
                                                 files=self.thesePgFiles)
        self.pg_source.fetch()
        # don't write metadata out with main graph in tests, because this screws up
        # triples counts during testing
        self.pg_source.write(write_metadata_in_main_graph=False)

    def test_version_level_source_file_triple(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI,
             self.iri_source,
             URIRef(self.source.files.get("test_file").get("url")))))
        self.assertTrue(len(triples) == 1, "missing version level file source triple")

    def test_version_level_source_version_download_timestamp(self):
        path_to_dl_file = '/'.join(
            (self.source.rawdir, self.source.files.get('test_file').get('file')))
        fstat = os.stat(path_to_dl_file)

        self.downloaded_file_timestamp = \
            datetime.utcfromtimestamp(fstat[ST_CTIME]).strftime("%Y%m%d")

        triples = list(self.source.dataset.graph.triples(
            (URIRef(self.theseFiles.get("test_file").get("url")),
             self.iri_retrieved_on,
             None)))
        self.assertTrue(len(triples) == 1,
                        "missing triple for ingest file retrieved on " +
                        "(download timestamp)")
        self.assertEqual(Literal(triples[0][2], datatype=XSD.date),
                         Literal(self.downloaded_file_timestamp, datatype=XSD.date),
                         "version level source version timestamp isn't " +
                         "the same as the download timestamp of the local file")

    def test_version_level_source_version_download_timestamp_setter(self):
        this_date = "1970-01-01"
        file_iri = self.source.files.get("test_file").get("url")
        self.source.dataset.set_ingest_source_file_version_retrieved_on(
            file_iri, this_date, datatype=XSD.date)
        triples = list(self.source.dataset.graph.triples(
            (URIRef(file_iri),
             self.iri_retrieved_on,
             Literal(this_date, datatype=XSD.date))))
        self.assertTrue(len(triples) == 1,
                        "ingest source file retrievedOn date not set correctly")

    def test_postgres_version_level_source_file_triple(self):
        triples = list(self.pg_source.dataset.graph.triples(
            (self.version_level_IRI,
             self.iri_source,
             URIRef(self.pg_source.files.get("test_file").get("url")))))
        self.assertTrue(len(triples) == 1, "missing version level file source triple")

    def test_postgres_version_level_source_version_download_timestamp(self):
        path_to_dl_file = '/'.join(
            (self.pg_source.rawdir, self.pg_source.files.get('test_file').get('file')))
        fstat = os.stat(path_to_dl_file)

        self.downloaded_file_timestamp = \
            datetime.utcfromtimestamp(fstat[ST_CTIME]).strftime("%Y%m%d")

        triples = list(self.pg_source.dataset.graph.triples(
            (URIRef(self.theseFiles.get("test_file").get("url")),
             self.iri_retrieved_on,
             None)))
        self.assertTrue(len(triples) == 1,
                        "missing triple for ingest file retrieved on " +
                        "(download timestamp)")
        self.assertEqual(Literal(triples[0][2], datatype=XSD.date),
                         Literal(self.downloaded_file_timestamp, datatype=XSD.date),
                         "version level source version timestamp isn't " +
                         "the same as the download timestamp of the local file")

    def test_postgres_version_level_source_version_download_timestamp_setter(self):
        this_date = "1970-01-01"
        file_iri = self.pg_source.files.get("test_file").get("url")
        self.pg_source.dataset.set_ingest_source_file_version_retrieved_on(
            file_iri, this_date, datatype=XSD.date)
        triples = list(self.pg_source.dataset.graph.triples(
            (URIRef(file_iri),
             self.iri_retrieved_on,
             Literal(this_date, datatype=XSD.date))))
        self.assertTrue(len(triples) == 1,
                        "ingest source file retrievedOn date not set correctly")


if __name__ == '__main__':
    unittest.main()


class FakeIngestClass(Source):
    """
    Fake ingest to test metadata in Dataset graph
    """

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 identifier=None,
                 ingest_url=None,
                 ingest_title=None,
                 ingest_desc=None,
                 ingest_logo=None,
                 license_url=None,
                 data_rights=None,
                 files=None
                 ):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            name=identifier,
            ingest_url=ingest_url,
            ingest_title=ingest_title,
            ingest_logo=ingest_logo,
            license_url=license_url,
            data_rights=data_rights
        )
        self.files = files

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit):
        pass

    # override fetch_from_url here so we don't need the network to "download" file
    def fetch_from_url(self,
                       remotefile,
                       localfile=None,
                       is_dl_forced=False,
                       headers=None):
        copyfile(
            self.files.get('test_file').get('path_to_mock_download_file'), localfile)
        return


class FakeIngestUsingPostgres(PostgreSQLSource):
    """
    Fake ingest to test metadata in Dataset graph, that uses PostgresSQLSource
    """

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 identifier=None,
                 ingest_url=None,
                 ingest_title=None,
                 ingest_desc=None,
                 ingest_logo=None,
                 license_url=None,
                 data_rights=None,
                 files=None
                 ):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            name=identifier,
            ingest_url=ingest_url,
            ingest_title=ingest_title,
            ingest_logo=ingest_logo,
            license_url=license_url,
            data_rights=data_rights
        )
        self.files = files

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit):
        pass

    def fetch_from_url(self,
                       remotefile,
                       localfile=None,
                       is_dl_forced=False,
                       headers=None):
        copyfile(
            self.files.get('test_file').get('path_to_mock_download_file'), localfile)
        return
