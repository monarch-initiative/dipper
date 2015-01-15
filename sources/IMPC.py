import csv
from utils import pysed
import os, datetime
from datetime import datetime
from stat import *
from rdflib import Graph, Literal
from rdflib.namespace import RDFS, OWL, RDF

from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from rdflib import Namespace, URIRef
import re
from utils.CurieUtil import CurieUtil


class IMPC(Source):

    files = {
        'impc' : {'file' : 'IMPC_genotype_phenotype.csv.gz',
                  'url' : 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/IMPC_genotype_phenotype.csv.gz'},
        'euro' : {'file' : 'EuroPhenome_genotype_phenotype.csv.gz',
                  'url' : 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/EuroPhenome_genotype_phenotype.csv.gz'},
        'mgd' : {'file' : 'MGP_genotype_phenotype.csv.gz',
                 'url' : 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/MGP_genotype_phenotype.csv.gz'}
    }

    namespaces = {
        'MP': 'http://purl.obolibrary.org/obo/MP_',
        'MA': 'http://purl.obolibrary.org/obo/MA_',
        'MGI' : 'http://www.informatics.jax.org/accession/MGI:'  #IMPC uses MGI identifiers for markers...
    }

    def __init__(self):
        Source.__init__(self, 'impc')

        #assemble all the curie mappings from the imported models
        self.namespaces.update(Assoc.curie_map)
        self.namespaces.update(Genotype.curie_map)

        #update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('impc', 'IMPC', 'http://www.mousephenotype.org')

        #source-specific warnings.  will be cleared when resolved.
        #print("WARN: we are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self):
        #this is fetching the standard files, not from the API/REST service
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],('/').join((self.rawdir,file['file'])))
            self.dataset.setFileAccessUrl(file['url'])
            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        self.dataset.setVersion(filedate)

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:

        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        #pysed.replace("\\\\", '', ('/').join((self.rawdir,self.files['geno']['file'])))

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")



        #self._process_genotype_features(('/').join((self.rawdir,self.files['geno']['file'])), self.outfile, self.graph, limit)
        #self._process_g2p(('/').join((self.rawdir,self.files['pheno']['file'])), self.outfile, self.graph, limit)
        #self._process_pubinfo(('/').join((self.rawdir,self.files['pubs']['file'])), self.outfile, self.graph, limit)

        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Found", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):

        #TODO make this more efficient
        #the problem with this implementation is that it creates many genotypes over and over, if the
        #same genotype has many features (on many rows) then the same genotype is recreated, then it must be
        #merged.  We should probably just create the genotype once, and then find the other
        #items that belong to that genotype.
        print("Processing Genotypes")

        print("INFO: Done with genotypes")
        return



    def _process_g2p(self, raw, out, g, limit=None):
        '''
        This module currently filters for only wild-type environments, which clearly excludes application
        of morpholinos.  Very stringent filter.  To be updated at a later time.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        '''
        print("Processing G2P")




        return



    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


