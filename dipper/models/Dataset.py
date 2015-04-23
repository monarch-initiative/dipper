__author__ = 'nicole'

from rdflib import Graph, Literal, URIRef, Namespace
from rdflib.namespace import RDF, DCTERMS, XSD, FOAF
import logging

logger = logging.getLogger(__name__)

class Dataset:

    #this will produce the metadata about a dataset
    #following the example laid out here: http://htmlpreview.github.io/?https://github.com/joejimbo/HCLSDatasetDescriptions/blob/master/Overview.html#appendix_1

    namespaces = {
        'dctypes' : 'http://purl.org/dc/dcmitype/',
        'pav' : 'http://purl.org/pav/',
        'dcat' : 'http://www.w3.org/ns/dcat#'
    }

    core_bindings = {'rdf': RDF, 'foaf': FOAF, 'xsd': XSD, 'dct': DCTERMS}

    def __init__(self,identifier, title, url, description=None, license_url=None, data_rights=None):
        DCTYPES=Namespace(self.namespaces['dctypes'])
        self.identifier = URIRef(':'+identifier)
        self.graph = Graph()
        self.load_bindings()
        self.graph.add((self.identifier,RDF['type'],DCTYPES['Dataset']))
        self.graph.add((self.identifier,DCTERMS['title'], Literal(title)))
        self.graph.add((self.identifier,DCTERMS['identifier'], Literal(identifier)))
        self.graph.add((self.identifier,FOAF['page'], URIRef(url)))
        #maybe in the future add the logo here:
        #schemaorg:logo <http://www.ebi.ac.uk/rdf/sites/ebi.ac.uk.rdf/files/resize/images/rdf/chembl_service_logo-146x48.gif> .

        #TODO add the licence info
        #FIXME:Temporarily making this in IF statement, can revert after all current resources are updated.
        if(license_url is not None):
            self.graph.add((self.identifier, DCTERMS['license'], URIRef(license_url)))
        else:
            logger.debug('No license provided.')
        if(data_rights is not None):
            self.graph.add((self.identifier, DCTERMS['rights'], Literal(data_rights)))
        else:
            logger.debug('No rights provided.')

        if (description is not None):
            self.graph.add((':'+identifier,DCTERMS['description'],description))
        return

    def load_bindings(self):
        for k in self.core_bindings:
            v=self.core_bindings[k]
            self.graph.bind(k, v)

        for k in self.namespaces.keys():
            v=self.namespaces[k]
            self.graph.bind(k, Namespace(v))

        return


    def setVersion(self,date_issued,version_id=None):
        PAV=Namespace(self.namespaces['pav'])

        #default to setting the version to the date, if no version specified
        if(version_id is None):
            version_id = date_issued
        self.version = URIRef(self.identifier+version_id)
        #todo verify date
        #date in YYYY-MM-DD format

        self.graph.add((self.identifier, DCTERMS['issued'], Literal(date_issued)))
        self.graph.add((self.version, DCTERMS['isVersionOf'], self.identifier))
        self.graph.add((self.version, PAV['version'], Literal(version_id)))

        print("INFO: setting version to",version_id)

        return

    def setFileAccessUrl(self, url):
        DCAT=Namespace(self.namespaces['dcat'])
        self.graph.add((self.identifier, DCAT['accessURL'],URIRef(url)))
        return

    def getGraph(self):
        return self.graph