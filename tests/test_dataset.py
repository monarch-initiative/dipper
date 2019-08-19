import datetime
import unittest
import logging
import os.path
from datetime import datetime
from rdflib import XSD

from rdflib import URIRef, Literal, Graph
from dipper.graph.RDFGraph import RDFGraph
from dipper.models.Dataset import Dataset
from dipper import curie_map as curiemap

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class DatasetTestCase(unittest.TestCase):
    """
    For testing metadata emitted by Dataset class

    Dataset creates a graph describing the metadata associated with the dataset in
    question, which should follow the HCLS specification for dataset descriptions
    https://www.w3.org/TR/2015/NOTE-hcls-dataset-20150514/
    """

    @classmethod
    def setUpClass(self):
        self.curie_map = curiemap.get()

        # parameters passed to code, to be returned in graph
        self.monarch_data_curie_prefix = "MonarchData"
        self.identifier = "fakeingest"
        self.ingest_description = "some ingest description"
        self.ingest_url = "http://fakeingest.com"
        self.ingest_title = "this ingest title"
        self.ingest_logo_url = "http://fakeingest.com/logo.png"
        self.license_url = "https://choosealicense.com/licenses/mit/"
        self.license_url_default = "https://project-open-data.cio.gov/unknown-license/"
        self.data_rights = "https://www.gnu.org/licenses/gpl-3.0.html"

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
        self.iri_triples_count = URIRef(self.curie_map.get("void") + "triples")
        self.iri_entities_count = URIRef(self.curie_map.get("void") + "entities")
        self.iri_distinct_subjects = URIRef(self.curie_map.get("void") +
                                            "distinctSubjects")
        self.iri_distinct_objects = URIRef(self.curie_map.get("void") +
                                           "distinctObjects")
        self.iri_properties_count = URIRef(self.curie_map.get("void") + "properties")

        self.iri_dipper = URIRef("https://github.com/monarch-initiative/dipper")
        self.iri_ttl_spec = URIRef("https://www.w3.org/TR/turtle/")

    @classmethod
    def tearDownClass(self):
        pass

    def setUp(self):
        self.dataset = Dataset(
            identifier=self.monarch_data_curie_prefix + ":" + self.identifier,
            ingest_title=self.ingest_title,
            ingest_url=self.ingest_url,
            ingest_logo=self.ingest_logo_url,
            ingest_description=self.ingest_description,
            license_url=self.license_url,
            data_rights=self.data_rights
        )

        # put all triples in a list for debugging below
        self.all_triples = list(self.dataset.graph.triples((None, None, None)))

    def tearDown(self):
        pass

    def test_dataset_has_graph(self):
        self.assertIsInstance(self.dataset.graph, Graph,
                              "dataset doesn't contain an RDF graph")

    def test_get_graph(self):
        self.assertIsInstance(
            self.dataset.get_graph(), RDFGraph,
            "get_graph() didn't return an RDF graph")

    def test_get_license(self):
        gpl2_iri = "https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"
        self.dataset.license_url=gpl2_iri
        self.assertEqual(self.dataset.get_license(),
                         gpl2_iri, "set_license didn't set license_url correctly")

    def test_set_citation(self):
        citation_iri =\
            "http://purl.obolibrary.org/obo/uberon/releases/2016-01-26/uberon.owl"
        self.dataset.set_citation(citation_iri)
        self.assertTrue(self.dataset.citation.issuperset([citation_iri]))
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI,
             URIRef(self.iri_cites_as_authority),
             URIRef(citation_iri))))
        self.assertTrue(len(triples) == 1, "missing citation triple")

    def test_set_ingest_source_file_version_num(self):
        this_version = "version1234"
        file_iri = "http://somefilesource.org/file.txt"
        self.dataset.set_ingest_source_file_version_num(file_iri, this_version)
        triples = list(self.dataset.graph.triples(
            (URIRef(file_iri), self.iri_version, Literal(this_version))))
        self.assertTrue(len(triples) == 1, "ingest source file version not set")

    def test_set_ingest_source_file_version_date(self):
        this_version = "1970-01-01"
        file_iri = "http://somefilesource.org/file.txt"
        self.dataset.set_ingest_source_file_version_date(file_iri, this_version)

        triples = list(self.dataset.graph.triples(
            (URIRef(file_iri),
             self.iri_version,
             Literal(this_version,datatype=XSD.date))))
        self.assertTrue(len(triples) == 1,
                        "ingest source file version not set with literal type of date")

    #
    # Test summary level triples:
    #
    def test_summary_level_type(self):
        triples = list(self.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_rdf_type, self.iri_dataset)))
        self.assertTrue(len(triples) == 1, "missing summary level type triple")

    def test_summary_level_title(self):
        triples = list(self.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_title, Literal(self.ingest_title))))
        self.assertTrue(len(triples) == 1, "missing summary level title triple")

    def test_summary_level_description(self):
        # by default, description is the class's docstring
        triples = list(self.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_description,
             Literal(self.ingest_description))))
        self.assertTrue(len(triples) == 1, "missing summary level description triple")

    def test_summary_level_publisher(self):
        triples = list(self.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_publisher, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1, "missing summary level publisher triple")

    def test_summary_level_source_web_page(self):
        triples = list(self.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_source, URIRef(self.ingest_url))))
        self.assertTrue(len(triples) == 1, "missing summary level source page triple")

    def test_summary_level_source_logo(self):
        triples = list(self.dataset.graph.triples(
            (self.summary_level_IRI, self.iri_logo, URIRef(self.ingest_logo_url))))
        self.assertTrue(len(triples) == 1, "missing summary level source logo triple")

    #
    # Test version level resource triples:
    #
    def test_version_level_type(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_rdf_type, self.iri_dataset)))
        self.assertTrue(len(triples) == 1, "missing version level type triple")

    def test_version_level_title(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_title, Literal(self.ingest_title))))
        self.assertTrue(len(triples) == 1, "missing version level title triple")

    def test_version_level_description(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_description,
             Literal(self.ingest_description))))
        self.assertTrue(len(triples) == 1, "missing version level description triple")

    def test_version_level_created(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_created, None)))
        self.assertTrue(len(triples) == 1, "didn't get exactly 1 version level " +
                                           "created triple")
        self.assertEqual(triples[0][2],
                         Literal(self.timestamp_date, datatype=XSD.date),
                         "version level created triple has the wrong timestamp")

    def test_version_level_version(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_version, None)))
        self.assertTrue(len(triples) == 1, "didn't get exactly one version level " +
                                           "version triple")
        self.assertEqual(triples[0][2],
                         Literal(self.timestamp_date, datatype=XSD.date),
                         "version level version triple has the wrong timestamp")

    def test_version_level_creator(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_creator, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1, "missing version level creator triple")

    def test_version_level_publisher(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_publisher, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1, "missing version level publisher triple")

    def test_version_level_isVersionOf(self):
        triples = list(self.dataset.graph.triples(
            (self.version_level_IRI, self.iri_is_version_of, self.summary_level_IRI)))
        self.assertTrue(len(triples) == 1, "missing version level isVersionOf triple")

    #
    # test distribution level triples
    #
    def test_distribution_level_dataset_type(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_rdf_type, self.iri_dataset)))
        self.assertTrue(len(triples) == 1, "missing version level type dataset triple")

    def test_distribution_level_distribution_type(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_rdf_type,
             self.iri_distribution)))
        self.assertTrue(len(triples) == 1,
                        "missing version level type distribution triple")

    def test_distribution_level_title(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_title,
             Literal(self.ingest_title))))
        self.assertTrue(len(triples) == 1, "missing version level type title triple")

    def test_distribution_level_description(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_description,
             Literal(self.ingest_description))))
        self.assertTrue(len(triples) == 1,
                        "missing version level type description triple")

    def test_distribution_level_created(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_created, None)))
        self.assertTrue(len(triples) == 1,
                        "didn't get exactly 1 version level type created triple")
        self.assertEqual(triples[0][2], Literal(self.timestamp_date, datatype=XSD.date))

    def test_distribution_level_version(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_version, None)))
        self.assertTrue(len(triples) == 1,
                        "didn't get exactly 1 version level type version triple")
        self.assertEqual(triples[0][2], Literal(self.timestamp_date, datatype=XSD.date))

    def test_distribution_level_creator(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_creator, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level creator triple")

    def test_distribution_level_publisher(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl, self.iri_publisher, self.iri_mi_org)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level publisher triple")

    def test_distribution_level_created_with(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_created_with,
             self.iri_dipper)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level createdWith triple")

    def test_distribution_level_format(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_format,
             self.iri_ttl_spec)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level format triple")

    def test_distribution_level_download_url(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_download_url,
             self.distribution_level_IRI_ttl)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level downloadUrl triple")

    def test_distribution_level_license_url(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_license,
             URIRef(self.license_url))))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level license triple")

    def test_distribution_level_data_rights(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_data_rights,
             URIRef(self.data_rights))))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level data rights triple")

    def test_distribution_level_no_license_url_default_value(self):
        self.dataset = Dataset(
            identifier=self.monarch_data_curie_prefix + ":" + self.identifier,
            ingest_title=self.ingest_title,
            ingest_url=self.ingest_url,
            ingest_logo=self.ingest_logo_url,
            ingest_description=self.ingest_description,
            license_url=None,
            data_rights=self.data_rights
        )
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_license,
             URIRef(self.license_url_default))))
        self.assertTrue(len(triples) == 1,
                        "distribution level default license triple not set")

    def test_distribution_level_triples_count(self):
        # feed test graph with 2 triples to self.dataset.compute_triples_statistics()
        exp_triples_count = 2
        test_ttl = "tests/resources/fakeingest/test_graph_simple.ttl"
        test_graph = RDFGraph()
        test_graph.parse(test_ttl,  format="turtle")

        self.dataset.compute_triples_statistics(test_graph)

        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_triples_count,
             None)))
        self.assertTrue(len(triples) == 1,
                        "didn't get exactly 1 distribution level triples count")
        self.assertEqual(triples[0][2], Literal(exp_triples_count),
                         "didn't get correct triples count")

    @unittest.skip("not implemented yet")
    def test_distribution_level_entities_count(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_entities_count,
             None)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level entities count")

    @unittest.skip("not implemented yet")
    def test_distribution_level_distinct_subject_count(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_distinct_subjects,
             None)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level distinct subject count")

    @unittest.skip("not implemented yet")
    def test_distribution_level_distinct_object_count(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_distinct_objects,
             None)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level distinct object count")

    @unittest.skip("not implemented yet")
    def test_distribution_level_properties_count(self):
        triples = list(self.dataset.graph.triples(
            (self.distribution_level_IRI_ttl,
             self.iri_properties_count,
             None)))
        self.assertTrue(len(triples) == 1,
                        "missing distribution level properties count")


if __name__ == '__main__':
    unittest.main()
