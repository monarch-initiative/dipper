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
     http://htmlpreview.github.io/?
     https://github.com/joejimbo/HCLSDatasetDescriptions/blob/master/Overview.html#appendix_1
     (mind the wrap)

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
        # TODO add the license info
        # FIXME:Temporarily making this in IF statement,
        #  can revert after all current resources are updated.
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

    def add_subj_obj_count_data(self, graph):
        """
        This takes a graph from a source (OMIM, MGI, whatever) and 
        added the count data for subject, object, and subject/object
        pairs to the dataset graph itself. The intention is to add
        easily ingestible info in the dataset TTL to be consumed by
        monarch UI, or other downstream sources
        :param graph: graph from source, after ingestion
        :return: 
        """
        for key in self.graph.sub_counts.keys():
            self.add_count_triples(key, self.globaltt['count'],
                                   self.graph.sub_counts[key])

        for key in self.graph.obj_counts.keys():
            self.add_count_triples(key, self.globaltt['count'],
                                   self.graph.sub_counts[key])
        # TODO: add triples for each subj/obj pair. need to change subj/obj's
        #  into labels: make_id("subj-obj") and
        #  model.addLabel(association_iri, "some label")
            
            
    def add_count_triples(self, subj, pred, obj):
        if subj is None:
            self.graph.addTriple('http://purl.obolibrary.org/obo/BFO_0000030',
                                 pred, obj, object_is_literal=True)
        else:
            self.graph.addTriple(subj.value, pred, obj, object_is_literal=True)
        
