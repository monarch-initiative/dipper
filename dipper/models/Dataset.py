import logging
from datetime import datetime
from dipper.graph.RDFGraph import RDFGraph
from dipper.graph.StreamedGraph import StreamedGraph
from dipper.models.Model import Model

__author__ = 'nlw'

LOG = logging.getLogger(__name__)


class Dataset:
    """
     this will produce the metadata about a dataset
     following the example laid out here:
     https://www.w3.org/TR/2015/NOTE-hcls-dataset-20150514/
     http://htmlpreview.github.io/?
     https://github.com/joejimbo/HCLSDatasetDescriptions/blob/master/Overview.html#appendix_1
     (mind the wrap)

     Summary level: The summary level provides a description of a dataset that is
     independent of a specific version or format. (e.g. the Monarch ingest of CTD)

     Version level: The version level captures version-specific characteristics of a
     dataset. (e.g. the 01-02-2018 ingest of CTD)

     Distribution level: The distribution level captures metadata about a specific form
     and version of a dataset. (e.g. the turtle file for 01-02-2018 ingest of CTD)

     Write out at least the following triples:

     [summary level resource] --- rdf:type ---> dctypes:Dataset
     [summary level resource] --- dct:title ---> title (literal)
     [summary level resource] --- dct:description ---> description (literal)
                                                    (use docstring from Source class)
     [summary level resource] --- dc:source ----> [source web page, e.g. omim.org]
     [summary level resource] --- schema:logo --> [source logo IRI]
     [summary level resource] --- foaf:page  ---> monarchinitiative.org
     ^^^ this is not the download page for the transformed/ingested data, e.g. a TTL
     file
     [summary level resource] --- dct:publisher ---> monarchinitiative.org

     [version level resource] --- rdf:type ---> dctypes:Dataset
     [version level resource] --- dct:isVersionOf ---> [summary level resource]
     [version level resource] --- pav:version --> [ingest timestamp]


     [version level resource] --- dc:source ----> [source file 1 IRI]
     [version level resource] --- dc:source ----> [source file 2 IRI]
     ...
     [version level resource] --- void:dataset -> [distribution level resource]


     [distribution level resource] --- rdf:type ---> dctypes:Dataset
     [distribution level resource] --- rdf:type ---> dcat:Distribution
     [distribution level resource] --- dcat:accessURL --> [MI ttl URL]
     [distribution level resource] --- dcat:accessURL --> [MI nt URL]
     ...

     [distribution level resource] --- void:triples --> [triples count (literal)]
     [distribution level resource] --- void:distinctSubjects -> [subject count (literal)]
     [distribution level resource] --- void:distinctObjects -> [object count (literal)]
     ...


     [source file 1 IRI] -- pav:version ---> [download date timestamp]
     [source file 2 IRI] -- pav:version ---> [source version (if set, optional)]
     [source file 2 IRI] -- pav:version ---> [download date timestamp]
     [source file 2 IRI] -- pav:version ---> [source version (if set, optional)]
     ...


    """

    def __init__(
            self,
            identifier,       # name? should be Archive url via Source
            title,
            url,
            ingest_desc=None,
            license_url=None,
            data_rights=None,
            graph_type='rdf_graph',     # rdf_graph, streamed_graph
            file_handle=None):

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
        if title is None:
            self.title = identifier
        else:
            self.title = title
        self.version = None
        self.date_issued = None

        # The data_accesed value is later used as an literal of properties
        # such as dcterms:issued, which needs to conform xsd:dateTime format.
        # TODO ... we need to have a talk about typed literals and SPARQL
        self.date_accessed = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

        self.citation = set()
        self.license_url = license_url
        self.model.addType(self.identifier, 'dctypes:Dataset')
        self.graph.addTriple(self.identifier, 'dcterms:title', title, True)
        self.graph.addTriple(
            self.identifier, 'dcterms:identifier', identifier, True)
        if url is not None:
            self.graph.addTriple(self.identifier, 'foaf:page', url)
        # maybe in the future add the logo here:
        # schemaorg:logo  <uri>

        if license_url is not None:
            self.graph.addTriple(
                self.identifier, 'dcterms:license', license_url)
        else:
            LOG.debug('No license provided.')
        if data_rights is not None:
            self.graph.addTriple(
                self.identifier, 'dcterms:rights',
                data_rights, object_is_literal=True)
        else:
            LOG.debug('No rights provided.')

        if ingest_desc is not None:
            self.model.addDescription(self.identifier, ingest_desc)
        return

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

    def setFileAccessUrl(self, url, is_object_literal=False):
        self.graph.addTriple(
            self.identifier, 'dcat:accessURL', url, is_object_literal)

    def getGraph(self):
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
