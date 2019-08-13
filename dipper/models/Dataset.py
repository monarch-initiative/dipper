import logging
from datetime import datetime
from rdflib import Literal, XSD
from dipper.graph.RDFGraph import RDFGraph
from dipper.graph.StreamedGraph import StreamedGraph
from dipper.models.Model import Model

__author__ = 'nlw'

LOG = logging.getLogger(__name__)


class Dataset:
    """
     this will produce the metadata about a dataset that is compliant with this
     specification:
     https://www.w3.org/TR/2015/NOTE-hcls-dataset-20150514/#s4_4

     Summary level: The summary level provides a description of a dataset that is
     independent of a specific version or format. (e.g. the Monarch ingest of CTD)
     CURIE for this is something like MonarchData:[SOURCE IDENTIFIER]

     Version level: The version level captures version-specific characteristics of a
     dataset. (e.g. the 01-02-2018 ingest of CTD)
     CURIE for this is something like MonarchData:[SOURCE IDENTIFIER_INGESTTIMESTAMP]

     Distribution level: The distribution level captures metadata about a specific form
     and version of a dataset (e.g. turtle file for 01-02-2018 ingest of CTD). There is
     a [distribution level resource] for each different downloadable file we emit,
     i.e. one for the TTL file, one for the ntriples file, etc.
     CURIE for this is like MonarchData:[SOURCE IDENTIFIER_INGESTTIMESTAMP].ttl
     or
     MonarchData:[SOURCE IDENTIFIER_INGESTTIMESTAMP].nt
     or
     MonarchData:[SOURCE IDENTIFIER_INGESTTIMESTAMP].[whatever file format]

     We write out at least the following triples:

     SUMMARY LEVEL TRIPLES:
     [summary level resource] - rdf:type -> dctypes:Dataset
     [summary level resource] - dct:title -> title (literal)
     [summary level resource] - dct:description -> description (literal)
                                                (use docstring from Source class)
     [summary level resource] - dcterms:source -> [source web page, e.g. omim.org]
     [summary level resource] - schemaorg:logo -> [source logo IRI]
     [summary level resource] - dct:publisher -> monarchinitiative.org
        n.b: about summary level resource triples:
        -- HCLS spec says we "should" link to our logo and web page, but I'm not,
        because it would confuse the issue of whether we are pointing to our logo/page
        or the logo/page of the data source for this ingest. Same below for
        [version level resource] and [distibution level resource] - I'm not linking to
        our page/logo down there either.
        - spec says we "should" include summary level triples describing Update
        frequency and SPARQL endpoint but I'm omitting this for now, because these are
        not clearly defined at the moment

     VERSION LEVEL TRIPLES:
     [version level resource] - rdf:type -> dctypes:Dataset
     [version level resource] - dct:title -> version title (literal)
     [version level resource] - dct:description -> version description (literal)
     [version level resource] - dct:created -> ingest timestamp [ISO 8601 compliant]
     [version level resource] - pav:version -> ingest timestamp (same one above)
     [version level resource] - dct:creator	-> monarchinitiative.org
     [version level resource] - dct:publisher -> monarchinitiative.org
     [version level resource] - dct:isVersionOf -> [summary level resource]
     [version level resource] - dcterms:source -> [source file 1 IRI]
     [version level resource] - dcterms:source -> [source file 2 IRI]
     ...

     [source file 1 IRI] - pav:retrievedOn -> [download date timestamp]
     [source file 2 IRI] - pav:version -> [source version (if set, optional)]
     [source file 2 IRI] - pav:retrievedOn -> [download date timestamp]
     [source file 2 IRI] - pav:version -> [source version (if set, optional)]
     ...

     [version level resource] - pav:createdWith -> [Dipper github URI]
     [version level resource] - void:dataset -> [distribution level resource]
        n.b: about version level resource triples:
        - spec says we "should" include Date of issue/dct:issued triple, but I'm not
        because it is redundant with this triple above:
        [version level resource] - dct:created -> time stamp
        and would introduce ambiguity and confusion if the two disagree. Same below
        for [distribution level resource] - dct:created -> tgiime stamp below
        Also omitting:
          - triples linking to our logo and page, see above.
          - License/dct:license triple, because we will make this triple via the
            [distribution level resource] below
          - Language/dct:language triple b/c it seems superfluous. Same below for
            [distribution level resource] - no language triple.
        - [version level resource] - pav:version triple is also a bit redundant
        with the pav:version triple below, but the spec requires both these triples
        - I'm omitting the [version level resource] -> pav:previousVersion because
        Dipper doesn't know this info for certain at run time. Same below for
        [distribution level resource] - pav:previousVersion.

     DISTRIBUTION LEVEL TRIPLES:
     [distribution level resource] - rdf:type -> dctypes:Dataset
     [distribution level resource] - rdf:type -> dcat:Distribution
     [distribution level resource] - dct:title -> distribution title (literal)
     [distribution level resource] - dct:description -> distribution description (lit.)
     [distribution level resource] - dct:created -> ingest timestamp[ISO 8601 compliant]
     [distribution level resource] - pav:version -> ingest timestamp (same as above)
     [distribution level resource] - dct:creator -> monarchinitiative.org
     [distribution level resource] - dct:publisher -> monarchinitiative.org
     [distribution level resource] - dct:license -> [license info, if available
                    otherwise indicate unknown]
     [distribution level resource] - pav:createdWith -> [Dipper github URI]
     [distribution level resource] - dct:format -> [IRI of ttl|nt|whatever spec]
     [distribution level resource] - dct:downloadURL -> [ttl|nt URI]
     [distribution level resource] - void:triples -> [triples count (literal)]
     [distribution level resource] - void:entities -> [entities count (literal)]
     [distribution level resource] - void:distinctSubjects -> [subject count (literal)]
     [distribution level resource] - void:distinctObjects -> [object count (literal)]
     [distribution level resource] - void:properties -> [properties count (literal)]
     ...

        n.b: about distribution level resource triples:
        - omitting Vocabularies used/void:vocabulary and Standards
        used/dct:conformTo triples, because they are described in the ttl file
        - also omitting Example identifier/idot:exampleIdentifier and
        Example resource/void:exampleResource, because we don't really have one
        canonical example of either - they're all very different.
        - [distribution level resource] - dct:created should have the exact same
        time stamp as this triple above:
        [version level resource] - dct:created -> time stamp
        - this [distribution level resource] - pav:version triple should have the
        same object as [version level resource] - pav:version triple above
        - Data source provenance/dct:source triples are above in the
        [version level resource]
        - omitting Byte size/dct:byteSize, RDF File URL/void:dataDump, and
        Linkset/void:subset triples because they probably aren't necessary for MI right
        now
        - these triples "should" be emitted, but we will do this in a later iteration:
        # of classes	void:classPartition	IRI
        # of literals	void:classPartition	IRI
        # of RDF graphs	void:classPartition	IRI
    """

    def __init__(
            self,
            identifier,       # name? should be Archive url via Source
            ingest_title,
            ingest_url,
            ingest_logo=None,
            ingest_description=None,
            license_url=None,
            data_rights=None,
            graph_type='rdf_graph',     # rdf_graph, streamed_graph
            file_handle=None
    ):

        if graph_type is None:
            self.graph = RDFGraph(None, identifier)
        elif graph_type == 'streamed_graph':
            self.graph = StreamedGraph(True, identifier, file_handle=file_handle)
        elif graph_type == 'rdf_graph':
            self.graph = RDFGraph(True, identifier)

        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        # TODO: move hard coded curies to translation table calls
        self.identifier = identifier
        self.citation = set()

        self.ingest_title = ingest_title
        if self.ingest_title is None:
            self.ingest_title = identifier

        self.ingest_url = ingest_url
        self.ingest_logo = ingest_logo
        self.ingest_description = ingest_description

        self.version = None
        self.date_issued = None

        # The data_accesed value is later used as an literal of properties
        # such as dcterms:issued, which needs to conform xsd:dateTime format.
        # TODO ... we need to have a talk about typed literals and SPARQL
        self.date_accessed = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

        self.license_url = license_url
        self.data_rights = data_rights

        # set HCLS resource CURIEs
        self.date_timestamp_iso_8601 = datetime.today().strftime("%Y%m%d")
        self.summary_level_curie = identifier
        self.version_level_curie = identifier + self.date_timestamp_iso_8601
        self.distribution_level_turtle_curie = self.version_level_curie + ".ttl"

        self.set_summary_level_triples()
        self.set_version_level_triples()
        self.set_distribution_level_triples()

        return

    def set_summary_level_triples(self):
        self.model.addType(self.summary_level_curie, 'dctypes:Dataset')
        self.graph.addTriple(self.summary_level_curie,
                             'dcterms:title',
                             self.ingest_title,
                             True)
        self.model.addTriple(self.summary_level_curie, 'dcterms:Publisher',
                             self.curie_map.get(""))
        self.model.addTriple(self.summary_level_curie,
                             "schemaorg:logo",
                             self.ingest_logo)
        self.graph.addTriple(self.summary_level_curie, 'dcterms:identifier',
                             self.identifier, True)
        if self.ingest_url is not None:
            self.graph.addTriple(self.summary_level_curie,
                                 "dcterms:source",
                                 self.ingest_url)
        if self.ingest_description is not None:
            self.model.addDescription(self.identifier, self.ingest_description)
        return

    def set_version_level_triples(self):
        self.model.addType(self.version_level_curie, 'dctypes:Dataset')
        # use version level curie itself for title and description
        self.graph.addTriple(self.version_level_curie, 'dcterms:title',
                             self.version_level_curie, True)
        self.model.addDescription(self.version_level_curie, self.version_level_curie)
        self.graph.addTriple(self.version_level_curie, 'dcterms:created',
                             Literal(self.date_timestamp_iso_8601, datatype=XSD.date))
        self.graph.addTriple(self.version_level_curie, 'pav:version',
                             Literal(self.date_timestamp_iso_8601, datatype=XSD.date))
        self.graph.addTriple(self.version_level_curie, 'dcterms:creator',
                             self.curie_map.get(""))  # eval's to MI.org
        self.graph.addTriple(self.version_level_curie, 'dcterms:Publisher',
                             self.curie_map.get(""))  # eval's to MI.org
        self.graph.addTriple(self.version_level_curie, 'dcterms:isVersionOf',
                             self.summary_level_curie)

    def set_distribution_level_triples(self):
        self.model.addType(self.distribution_level_turtle_curie, 'dctypes:Dataset')
        self.model.addType(self.distribution_level_turtle_curie, 'dcat:Distribution')
        self.graph.addTriple(self.distribution_level_turtle_curie, 'dcterms:title',
                           self.distribution_level_turtle_curie, True)
        self.model.addDescription(self.distribution_level_turtle_curie,
                                  self.distribution_level_turtle_curie)
        self.graph.addTriple(self.distribution_level_turtle_curie, 'pav:version',
                             Literal(self.date_timestamp_iso_8601, datatype=XSD.date))
        self.graph.addTriple(self.distribution_level_turtle_curie, 'dcterms:created',
                             Literal(self.date_timestamp_iso_8601, datatype=XSD.date))
        self.graph.addTriple(self.distribution_level_turtle_curie, 'dcterms:creator',
                             self.curie_map.get(""))  # eval's to MI.org
        self.graph.addTriple(self.distribution_level_turtle_curie, 'dcterms:Publisher',
                             self.curie_map.get(""))  # eval's to MI.org
        self.graph.addTriple(self.distribution_level_turtle_curie, 'pav:createdWith',
                             "https://github.com/monarch-initiative/dipper")
        self.graph.addTriple(self.distribution_level_turtle_curie, 'dcterms:format',
                             "https://www.w3.org/TR/turtle/")
        self.graph.addTriple(self.distribution_level_turtle_curie,
                             'dcterms:downloadURL',
                             self.distribution_level_turtle_curie)
        if self.license_url is None:
            self.graph.addTriple(self.distribution_level_turtle_curie,
                                     'dcterms:license',
                                     "https://project-open-data.cio.gov/unknown-license/")
        else:
            self.graph.addTriple(self.distribution_level_turtle_curie,
                                 'dcterms:license',
                                 self.license_url)

        if self.data_rights is not None:
            self.graph.addTriple(self.distribution_level_turtle_curie,
                                 'dcterms:rights',
                                 self.data_rights)

    def set_ingest_source_file_version_num(self, file_iri, version, literal_type=None):
        self.graph.addTriple(file_iri, 'pav:version', version, object_is_literal=True)

    def set_ingest_source_file_version_date(self, file_iri, date):
        self.graph.addTriple(file_iri, 'pav:version', date,
                             object_is_literal=True, literal_type=XSD.date)


    def setVersion(self, date_issued, version_id=None):
        """
        Legacy function...
            should use the other set_* for version and date

        as of 2016-10-20  used in:

        dipper/sources/HPOAnnotations.py 139:
        dipper/sources/CTD.py             99:
        dipper/sources/BioGrid.py        100:
        dipper/sources/MGI.py            255:
        dipper/sources/EOM.py             93:
        dipper/sources/Coriell.py        200:
        dipper/sources/MMRRC.py           77:

        # TODO set as deprecated

        :param date_issued:
        :param version_id:
        :return:

        """

        if date_issued is not None:
            self.set_date_issued(date_issued)
        elif version_id is not None:
            self.set_version_by_num(version_id)
        else:
            LOG.error("date or version not set!")
            # TODO throw error
            return

        if version_id is not None:
            self.set_version_by_num(version_id)
        else:
            LOG.info("set version to %s", self.version)
            self.set_version_by_date(date_issued)

        LOG.info("set version to %s", self.version)

        return

    def set_date_issued(self, date_issued):

        self.date_issued = date_issued
        self.graph.addTriple(
            self.identifier, 'dcterms:issued', date_issued,
            object_is_literal=True)
        LOG.info("setting date to %s", date_issued)

        return

    def set_version_by_date(self, date_issued=None):
        """
        This will set the version by the date supplied,
        the date already stored in the dataset description,
        or by the download date (today)
        :param date_issued:
        :return:
        """

        if date_issued is not None:
            dat = date_issued
        elif self.date_issued is not None:
            dat = self.date_issued
        else:
            dat = self.date_accessed
            LOG.info(
                "No date supplied, using download timestamp for date_issued")
        LOG.info("setting version by date to: %s", dat)
        self.set_version_by_num(dat)

        return

    def set_version_by_num(self, version_num):

        self.version = self.identifier+version_num
        self.graph.addTriple(self.version, 'dcterms:isVersionOf', self.identifier)
        self.graph.addTriple(
            self.version, 'pav:version', version_num, object_is_literal=True)

        LOG.info("setting version to %s", self.version)

        # set the monarch-generated-version of the resource-version
        # TODO sync this up with the ontology version
        if version_num != self.date_accessed:
            dipperized_version = ':' + str(self.date_accessed)
            self.graph.addTriple(
                dipperized_version, 'dcterms:isVersionOf',
                "MonarchData:" + self.identifier + ".ttl")  # fix suffix
            self.graph.addTriple(
                dipperized_version, 'pav:version',
                self.date_accessed, object_is_literal=True)
            self.graph.addTriple(
                dipperized_version, 'dcterms:issued', self.date_accessed,
                object_is_literal=True, literal_type="xsd:dateTime")
        return

    def set_ingest_source(self, url,
                          predictate='dcterms:source',
                          is_object_literal=False):
        """
        This method writes a triple to the dataset graph indicating that the ingest
        used a file or resource at [url] during the ingest.

        Triple emitted is version_level_curie dcterms:source [url]

        This triple is likely to be redundant, since this triple should also be emitted
        in Source.get_files() as remote files are being retrieved. This method is
        provided as a convenience method for sources that do their own downloading
        of files.

        :param url: a remote resource used as a source during ingest
        :param predicate: the predicate to use for the triple ["dcterms:source"]
                from spec (https://www.w3.org/TR/2015/NOTE-hcls-dataset-20150514/)
                "Use dct:source when the source dataset was used in whole or in part.
                Use pav:retrievedFrom when the source dataset was used in whole and was
                not modified from its original distribution. Use prov:wasDerivedFrom
                when the source dataset was in whole or in part and was modified from
                its original distribution."
        :return: None
        """
        self.graph.addTriple(
            self.version_level_curie, predictate, url,
            object_is_literal=is_object_literal)

    def get_graph(self):
        return self.graph

    def set_license(self, license_url):
        self.license_url = license_url
        return

    def get_license(self):
        return self.license_url

    def set_citation(self, citation_id):

        self.citation.add(citation_id)
        # TODO
        # model.addTriple(self.identifier, 'cito:citeAsAuthority', citation_id)

        return
