#!/usr/bin/env python3

import os
import datetime
import unittest
import logging
from datetime import datetime
from stat import ST_CTIME
from rdflib import XSD
from dipper.sources.Source import Source
from dipper.graph.RDFGraph import RDFGraph
from dipper import curie_map as curiemap
from rdflib import URIRef, Literal

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

    @classmethod
    def setUpClass(self):
        self.curie_map = curiemap.get()

        # parameters passed to code, to be returned in graph
        self.identifier = "fakeingest"
        self.ingest_url = "http://fakeingest.com"
        self.ingest_title = "this ingest title"
        self.ingest_logo_url = "http://fakeingest.com/logo.png"
        self.license_url = "https://choosealicense.com/licenses/mit/"
        self.license_url_default = "https://project-open-data.cio.gov/unknown-license/"
        self.data_rights = "https://www.gnu.org/licenses/gpl-3.0.html"

        # The following to be used in both DatasetTestCase and FakeIngestClass
        # using robots.txt b/c it's a trivially small file/empty file that we control
        self.theseFiles = {
            'test_file': {
                'file': 'test_file.txt',
                'url': 'https://data.monarchinitiative.org/robots.txt'},
        }

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

        # expected things:
        self.expected_curie_prefix = "MonarchData"
        self.timestamp_date = datetime.today().strftime("%Y%m%d")

        # expected summary level IRI
        self.summary_level_IRI = URIRef(self.curie_map.get(self.expected_curie_prefix)
                                        + self.identifier)
        # expected version level IRI
        self.version_level_IRI = URIRef(self.summary_level_IRI + self.timestamp_date)

        # expected distribution level IRI (for ttl resource)
        self.distribution_level_IRI_ttl = URIRef(self.version_level_IRI + ".ttl")

        # expected timestamp for version level "version" triple
        # downloaded file should end up here:
        self.downloaded_file_path = \
            '/'.join((self.source.rawdir, self.theseFiles.get("test_file").get("file")))
        fstat = os.stat(self.downloaded_file_path)
        self.downloaded_file_timestamp = \
            datetime.utcfromtimestamp(fstat[ST_CTIME]).strftime("%Y%m%d")

        # set expected IRIs for predicates and other things
        self.iri_rdf_type = URIRef(self.curie_map.get("rdf") + "type")
        self.iri_title = URIRef(self.curie_map.get("dcterms") + "title")
        self.iri_dataset = URIRef(self.curie_map.get("dctypes") + "Dataset")
        self.iri_description = URIRef(self.curie_map.get("dc") + "description")
        self.iri_publisher = URIRef(self.curie_map.get("dcterms") + "Publisher")
        self.iri_source = URIRef(self.curie_map.get("dcterms") + "source")
        self.iri_logo = URIRef(self.curie_map.get("schemaorg") + "logo")
        self.iri_mi_org = URIRef("https://monarchinitiative.org/")
        self.iri_created = URIRef(self.curie_map.get("dcterms") + "created")
        self.iri_version = URIRef(self.curie_map.get("pav") + "version")
        self.iri_retrieved_on = URIRef(self.curie_map.get("pav") + "retrievedOn")
        self.iri_creator = URIRef(self.curie_map.get("dcterms") + "creator")
        self.iri_is_version_of = URIRef(self.curie_map.get("dcterms") + "isVersionOf")
        self.iri_distribution = URIRef(self.curie_map.get("dcat") + "Distribution")
        self.iri_created_with = URIRef(self.curie_map.get("pav") + "createdWith")
        self.iri_format = URIRef(self.curie_map.get("dcterms") + "format")
        self.iri_download_url = URIRef(self.curie_map.get("dcterms") + "downloadURL")
        self.iri_license = URIRef(self.curie_map.get("dcterms") + "license")
        self.iri_data_rights = URIRef(self.curie_map.get("dcterms") + "rights")
        self.iri_cites_as_authority = URIRef(self.curie_map.get("cito") +
                                             "citesAsAuthority")

        self.iri_dipper = URIRef("https://github.com/monarch-initiative/dipper")
        self.iri_ttl_spec = URIRef("https://www.w3.org/TR/turtle/")

        # put all triples in a list for debugging below
        self.all_triples = list(self.source.dataset.graph.triples((None, None, None)))

    @classmethod
    def tearDownClass(self):
        pass

    def test_has_dataset_attribute(self):
        self.assertTrue(hasattr(self.source, "dataset"),
                        "source object doesn't have dataset attribute")

    def test_dataset_has_graph(self):
        self.assertIsInstance(self.source.graph, RDFGraph,
                              "dataset doesn't contain an RDF graph")

    def test_set_ingest_source_file_version_num(self):
        this_version = "version1234"
        file_iri = self.source.files.get("test_file").get("url")
        self.source.dataset.set_ingest_source_file_version_num(file_iri, this_version)
        triples = list(self.source.dataset.graph.triples(
            (URIRef(file_iri), self.iri_version, Literal(this_version))))
        self.assertTrue(len(triples) == 1, "ingest source file version not set")

    def test_set_ingest_source_file_version_date(self):
        this_version = "1970-01-01"
        file_iri = self.source.files.get("test_file").get("url")
        self.source.dataset.set_ingest_source_file_version_date(file_iri, this_version)

        triples = list(self.source.dataset.graph.triples(
            (URIRef(file_iri),
             self.iri_version,
             Literal(this_version,datatype=XSD.date))))
        self.assertTrue(len(triples) == 1,
                        "ingest source file version not set with literal type of date")

    def test_get_graph(self):
        self.assertIsInstance(
            self.source.dataset.get_graph(), RDFGraph,
            "get_graph() didn't return an RDF graph")

    def test_get_license(self):
        gpl2_iri = "https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"
        self.source.dataset.license_url=gpl2_iri
        self.assertEqual(self.source.dataset.get_license(),
                         gpl2_iri, "set_license didn't set license_url correctly")

    def test_set_citation(self):
        citation_iri =\
            "http://purl.obolibrary.org/obo/uberon/releases/2016-01-26/uberon.owl"
        self.source.dataset.set_citation(citation_iri)
        self.assertTrue(self.source.dataset.citation.issuperset([citation_iri]))
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI,
             URIRef(self.iri_cites_as_authority),
             URIRef(citation_iri))))
        self.assertTrue(len(triples) == 1, "missing citation triple")

    #
    # Test summary level triples:
    #
    def test_summary_level_type(self):
        triples = list(self.source.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_rdf_type, self.iri_dataset)))
        self.assertTrue(len(triples) == 1, "missing summary level type triple")

    def test_summary_level_title(self):
        triples = list(self.source.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_title, Literal(self.ingest_title))))
        self.assertTrue(len(triples) == 1, "missing summary level title triple")

    def test_summary_level_description(self):
        # by default, description is the class's docstring
        triples = list(self.source.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_description,
             Literal("Fake ingest to test metadata in Dataset graph"))))
        self.assertTrue(len(triples) == 1, "missing summary level description triple")

    def test_summary_level_publisher(self):
        triples = list(self.source.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_publisher, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1, "missing summary level publisher triple")

    def test_summary_level_source_web_page(self):
        triples = list(self.source.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_source, URIRef(self.ingest_url))))
        self.assertTrue(len(triples) == 1, "missing summary level source page triple")

    def test_summary_level_source_logo(self):
        triples = list(self.source.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_logo, URIRef(self.ingest_logo_url))))
        self.assertTrue(len(triples) == 1, "missing summary level source logo triple")

    #
    # Test version level resource triples:
    #
    def test_version_level_type(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_rdf_type, self.iri_dataset)))
        self.assertTrue(len(triples) == 1, "missing version level type triple")

    def test_version_level_title(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_title, Literal(self.version_level_IRI))))
        self.assertTrue(len(triples) == 1, "missing version level title triple")

    def test_version_level_description(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_description,
             Literal(self.version_level_IRI))))
        self.assertTrue(len(triples) == 1, "missing version level description triple")

    def test_version_level_created(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_created, None)))
        self.assertTrue(len(triples) == 1, "missing version level created triple")
        self.assertEqual(triples[0][2],
                         Literal(self.timestamp_date, datatype=XSD.date),
                         "version level created triple has the wrong timestamp")

    def test_version_level_version(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_version, None)))
        self.assertTrue(len(triples) == 1, "missing version level version triple")
        self.assertEqual(triples[0][2],
                         Literal(self.timestamp_date, datatype=XSD.date),
                         "version level version triple has the wrong timestamp")

    def test_version_level_creator(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_creator, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1, "missing version level creator triple")

    def test_version_level_publisher(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_publisher, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1, "missing version level publisher triple")

    def test_version_level_isVersionOf(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI, self.iri_is_version_of, self.summary_level_IRI)))
        self.assertTrue(len(triples) == 1, "missing version level isVersionOf triple")

    def test_version_level_source_file_triple(self):
        triples = list(self.source.dataset.graph.triples(
            (self.version_level_IRI,
             self.iri_source,
             URIRef(self.source.files.get("test_file").get("url"))
        )))
        self.assertTrue(len(triples) == 1, "missing version level file source triple")

    def test_version_level_source_version_download_timestamp(self):
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
        # test setter too
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

    #
    # test distribution level triples
    #
    def test_distribution_level_dataset_type(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_rdf_type, self.iri_dataset)))
        self.assertTrue(len(triples) == 1, "missing version level type dataset triple")

    def test_distribution_level_distribution_type(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_rdf_type,
             self.iri_distribution)))
        self.assertTrue(len(triples) == 1,
                        "missing version level type distribution triple")

    def test_distribution_level_title(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_title,
             Literal(self.distribution_level_IRI_ttl))))
        self.assertTrue(len(triples) == 1, "missing version level type title triple")

    def test_distribution_level_description(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_description,
             Literal(self.distribution_level_IRI_ttl))))
        self.assertTrue(len(triples) == 1,
                        "missing version level type description triple")

    def test_distribution_level_created(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_created, None)))
        self.assertTrue(len(triples) == 1,
                        "missing version level type created triple")
        self.assertEqual(triples[0][2], Literal(self.timestamp_date, datatype=XSD.date))

    def test_distribution_level_version(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_version, None)))
        self.assertTrue(len(triples) == 1,
                        "missing version level type version triple")
        self.assertEqual(triples[0][2], Literal(self.timestamp_date, datatype=XSD.date))

    def test_distribution_level_creator(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_creator, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level creator triple")

    def test_distribution_level_publisher(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_publisher, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level publisher triple")

    def test_distribution_level_created_with(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_created_with,
             self.iri_dipper)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level createdWith triple")

    def test_distribution_level_format(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_format,
             self.iri_ttl_spec)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level format triple")

    def test_distribution_level_download_url(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_download_url,
             self.distribution_level_IRI_ttl)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level downloadUrl triple")

    def test_distribution_level_license_url(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_license,
             URIRef(self.license_url))))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level license triple")

    def test_distribution_level_data_rights(self):
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_data_rights,
             URIRef(self.data_rights))))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level data rights triple")

    def test_distribution_level_no_license_url_default_value(self):
        self.source = FakeIngestClass("rdf_graph",
                                      are_bnodes_skolemized=False,
                                      identifier=self.identifier,
                                      ingest_url=self.ingest_url,
                                      ingest_title=self.ingest_title,
                                      ingest_logo=self.ingest_logo_url,
                                      license_url=None,  # not set
                                      data_rights=self.data_rights,
                                      files=self.theseFiles)
        self.source.fetch()
        triples = list(self.source.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_license,
             URIRef(self.license_url_default))))
        self.assertTrue(len(triples) == 1,
                        "distribution level default license triple not set")


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
        return

    def parse(self):
        pass
