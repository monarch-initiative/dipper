import csv, os, datetime
from datetime import datetime
from stat import *
import urllib
import gzip,os.path
from utils import fileutils

from sources.Source import Source

from models.D2PAssoc import D2PAssoc
from models.DispositionAssoc import DispositionAssoc
from rdflib import Namespace
from models.Dataset import Dataset
from models.Assoc import Assoc


'''
#down the entry files
cat $DATAFILE | while read fileline ; do
  log $fileline
  # delete file before create it again
  rm -f entry_$fileline.xml
  sleep 0.5
  wget -O entry_$fileline.xml "http://api.omim.org:8000/api/entry?include=all&apiKey=$apiKey&mimNumber=$fileline"
done

'''
class OMIM(Source):

    files = {
        'all' : {
            'file' : 'omim.txt.gz',
            'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/omim.txt.Z'
        },
        'genemap' : {
            'file': 'genemap.txt',
            'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/genemap'
        },
        'genemapkey' : {
            'file': 'genemapkey.txt',
            'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/genemap.key'

        },
        'morbidmap' : {
            'file': 'morbidmap.txt',
            'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/morbidmap'
        }
    }

    disease_prefixes = {
        'OMIM' : 'http://purl.obolibrary.org/obo/OMIM_',
    }

    curie_map = {}

    def __init__(self):
        Source.__init__(self, 'omim')

        self.curie_map.update(D2PAssoc.curie_map)
        self.curie_map.update(DispositionAssoc.curie_map)
        self.curie_map.update(self.disease_prefixes)

        self.load_bindings()

        self.dataset = Dataset('omim', 'Online Mendelian Inheritance in Man', 'http://www.omim.org')

        #data-source specific warnings (will be removed when issues are cleared)
        #print()

        return

    def fetch(self):
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
        For this resource, this currently includes: (Nothing)
        :return: None
        '''
        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.curie_map.keys():
            v=self.curie_map[k]
            self.graph.bind(k, Namespace(v))
        return

    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows")

        print("Parsing files...")

        return


    def verify(self):
        status = True
        #TODO verify some kind of relationship that should be in the file

        return status