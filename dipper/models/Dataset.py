from datetime import datetime

__author__ = 'nlw'

from rdflib import Graph, Literal, URIRef, Namespace
from rdflib.namespace import RDF, DCTERMS, XSD, FOAF
import logging

logger = logging.getLogger(__name__)


class Dataset:

    # this will produce the metadata about a dataset
    # following the example laid out here:
    # http://htmlpreview.github.io/?https://github.com/joejimbo/HCLSDatasetDescriptions/blob/master/Overview.html#appendix_1

    namespaces = {
        'dctypes': 'http://purl.org/dc/dcmitype/',
        'pav': 'http://purl.org/pav/',
        'dcat': 'http://www.w3.org/ns/dcat#'
    }

    core_bindings = {'rdf': RDF, 'foaf': FOAF, 'xsd': XSD, 'dct': DCTERMS}

    def __init__(self, identifier, title, url, description=None, license_url=None, data_rights=None):
        DCTYPES = Namespace(self.namespaces['dctypes'])
        self.identifier = URIRef(':'+identifier)
        self.version = None
        self.date_issued = None
        self.date_accessed = None
        self.set_access_date()
        self.license = license_url
        self.graph = Graph()
        self.load_bindings()
        self.graph.add((self.identifier, RDF['type'], DCTYPES['Dataset']))
        self.graph.add((self.identifier, DCTERMS['title'], Literal(title)))
        self.graph.add((self.identifier, DCTERMS['identifier'], Literal(identifier)))
        self.graph.add((self.identifier, FOAF['page'], URIRef(url)))
        # maybe in the future add the logo here:
        # schemaorg:logo <http://www.ebi.ac.uk/rdf/sites/ebi.ac.uk.rdf/files/resize/images/rdf/chembl_service_logo-146x48.gif> .

        # TODO add the licence info
        # FIXME:Temporarily making this in IF statement, can revert after all current resources are updated.
        if license_url is not None:
            self.graph.add((self.identifier, DCTERMS['license'], URIRef(license_url)))
        else:
            logger.debug('No license provided.')
        if data_rights is not None:
            self.graph.add((self.identifier, DCTERMS['rights'], Literal(data_rights)))
        else:
            logger.debug('No rights provided.')

        if description is not None:
            self.graph.add((':'+identifier, DCTERMS['description'], description))
        return

    def load_bindings(self):
        for k in self.core_bindings:
            v = self.core_bindings[k]
            self.graph.bind(k, v)

        for k in self.namespaces.keys():
            v = self.namespaces[k]
            self.graph.bind(k, Namespace(v))

        return

    def setVersion(self, date_issued, version_id=None):
        """
        Legacy function...  should use the other set_* for version and date

        # TODO set as deprecated
        :param date_issued:
        :param version_id:
        :return:
        """
        if date_issued is not None:
            self.set_date_issued(date_issued)
        elif version_id is not None:
            # this shouldn't happen
            self.set_version_by_num(version_id)
        else:
            logger.error("No date or version set!")
            # TODO throw error
            return

        if version_id is not None:
            self.set_version_by_num(version_id)
        else:
            self.set_version_by_date(date_issued)

        logger.info("set version to %s", self.version)

        return

    def set_date_issued(self, date_issued):

        self.date_issued = date_issued
        self.graph.add((self.identifier, DCTERMS['issued'], Literal(date_issued)))
        logger.info("setting date to %s", date_issued)

        return

    def set_version_by_date(self, date_issued=None):
        """
        This will set the version by the date supplied, the date already stored in the dataset description,
        or by the download date (today)
        :param date_issued:
        :return:
        """

        if date_issued is not None:
            d = date_issued
        elif self.date_issued is not None:
            d = self.date_issued
        else:
            d = self.date_accessed
            logger.info("No date supplied for setting version; using download timestamp for date_issued")

        logger.info("setting version by date")
        self.set_version_by_num(d)

        return

    def set_version_by_num(self, version_num):
        PAV = Namespace(self.namespaces['pav'])

        self.version = URIRef(self.identifier+version_num)
        self.graph.add((self.version, DCTERMS['isVersionOf'], self.identifier))
        self.graph.add((self.version, PAV['version'], Literal(version_num)))

        logger.info("setting version to %s", self.version)

        # set the monarch-generated-version of the resource-version
        # TODO sync this up with the ontology version
        if version_num != self.date_accessed:
            self.dipperized_version = URIRef('monarch'+str(self.date_accessed))
            self.graph.add((self.dipperized_version, DCTERMS['isVersionOf'], self.version))
            self.graph.add((self.dipperized_version, PAV['version'], Literal(self.date_accessed)))
            self.graph.add((self.dipperized_version, DCTERMS['issued'], Literal(self.date_accessed, datatype=XSD.dateTime)))

        return

    def set_access_date(self):

        t = datetime.now()
        t_string = t.strftime("%Y-%m-%d-%H-%M")
        d = t_string
        self.date_accessed = d
        logger.info("Setting date of access to %s", self.date_accessed)

        return

    def setFileAccessUrl(self, url):
        DCAT = Namespace(self.namespaces['dcat'])
        self.graph.add((self.identifier, DCAT['accessURL'], URIRef(url)))
        return

    def getGraph(self):
        return self.graph

    def set_license(self, license):
        self.license = license
        return

    def get_license(self):

        return self.license