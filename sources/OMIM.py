import os
from stat import *
import urllib
from urllib import request
import re
import time
from datetime import datetime
import gzip,os.path
import json
from rdflib import Graph, Literal, URIRef, Namespace
from rdflib.namespace import RDF, RDFS, OWL, DC
from sources.Source import Source

from models.D2PAssoc import D2PAssoc
from models.DispositionAssoc import DispositionAssoc
from models.Dataset import Dataset
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
import config

class OMIM(Source):
    '''
     OMIM is an unusual source.  We can get lots of the disease-gene associations, including allelic variants
     from their ftp site, which is obtainable anonymously.  However, more detailed information is available
     via their API.  So, we pull the basic files from their ftp site, extract the omim identifiers,
     then query their API in batch.  (Note this requires an apiKey, which is not stored in the repo,
     but in a separate conf.json file.)
    '''

    files = {
        'all' : {
            'file' : 'omim.txt.gz',
            'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/omim.txt.Z'
        },
        #'genemap' : {
        #    'file': 'genemap.txt',
        #    'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/genemap'
        #},
        #'genemapkey' : {
        #    'file': 'genemapkey.txt',
        #    'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/genemap.key'
        #},
        #'morbidmap' : {
        #   'file': 'morbidmap.txt',
        #    'url' : 'ftp://anonymous:info%40monarchinitiative.org@ftp.omim.org/OMIM/morbidmap'
        #}
    }

    disease_prefixes = {
        'OMIM' : 'http://purl.obolibrary.org/obo/OMIM_',
    }

    curie_map = {}

    OMIM_API = "http://api.omim.org:8000/api"

    def __init__(self):
        Source.__init__(self, 'omim')

        self.curie_map.update(D2PAssoc.curie_map)
        self.curie_map.update(DispositionAssoc.curie_map)
        self.curie_map.update(self.disease_prefixes)

        self.load_bindings()

        self.dataset = Dataset('omim', 'Online Mendelian Inheritance in Man', 'http://www.omim.org')

        #data-source specific warnings (will be removed when issues are cleared)
        #print()

        #check if config exists; if it doesn't, error out and let user know
        if (not (('keys' in config.get_config()) and ('omim' in config.get_config()['keys']))):
            print("ERROR: not configured with API key.")
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

        omimids = []  #to store the set of omim identifiers
        omimids = self._get_omim_ids()

        omimparams = {
            'format' : 'json',
            'include' : 'all',
        }
        #you will need to add the API key into the conf.json file, like:
        # keys : { 'omim' : '<your api key here>' }
        omimparams.update({'apiKey' : config.get_config()['keys']['omim']})

        #http://api.omim.org/api/entry?mimNumber=100100&include=all

        g = self.graph

        gu = GraphUtils(self.curie_map)
        cu = CurieUtil(self.curie_map)

        it=0  #for counting

        #note that you can only do request batches of 20
        #see info about "Limits" at http://omim.org/help/api
        groupsize=20
        if (limit is not None):
            #just in case the limit is larger than the number of records, max it out
            max = min((limit,omimids.__len__()))

        #max = 10 #for testing

        #TODO write the json to local files - make the assumption that downloads within 24 hrs are the same
        #now, loop through the omim numbers and pull the records as json docs
        while it < max:
            end=min((max,it+groupsize))
            #iterate through the omim ids list, and fetch from the OMIM api in batches of 20
            omimparams.update({'mimNumber' : (',').join(omimids[it:end])})
            p = urllib.parse.urlencode(omimparams)
            url = ('/').join((self.OMIM_API,'entry'))+'?%s' % p
            print ('fetching:',(',').join(omimids[it:end]))
            #print ('fetching:',('/').join((self.OMIM_API,'entry'))+'?%s' % p)
            it+=groupsize

            d = urllib.request.urlopen(url)
            resp = d.read().decode()
            request_time = datetime.now()

            myjson = json.loads(resp)
            entries = myjson['omim']['entryList']

            for e in entries:

                #get the numbers, labels, and descriptions
                omimnum = e['entry']['mimNumber']
                titles = e['entry']['titles']
                label = titles['preferredTitle']
                description = self._get_description(e['entry'])

                if (e['entry']['status'] == 'removed'):
                    n = URIRef(cu.get_uri('OMIM:'+str(omimnum)))
                    #make this a deprecated class
                    g.add((n,RDF['type'],OWL['DeprecatedClass']))
                else:
                    gu.addClassToGraph(g,'OMIM:'+str(omimnum),label,None,description)

                    #check if moved, if so, make it an equivalent class to the other thing
                    if (e['entry']['status'] == 'moved'):
                        gu.addEquivalentClass(g,'OMIM:'+str(omimnum),'OMIM:'+str(e['entry']['movedTo']))


            #can't have more than 4 req per sec,
            #so wait the remaining time, if necessary
            dt=datetime.now()-request_time
            rem=0.25-dt.total_seconds()
            if (rem > 0):
                print("INFO: waiting",str(rem),'s')
                time.sleep(rem/1000)

        ##### Write it out #####
        filewriter = open(self.outfile, 'w')
        self.load_core_bindings()
        self.load_bindings()

        print("Finished parsing files. Writing turtle to",self.outfile)
        print(g.serialize(format="turtle").decode(),file=filewriter)
        filewriter.close()


        return

    def verify(self):
        status = True
        #TODO verify some kind of relationship that should be in the file

        return status

    def _get_description(self,entry):
        d = None
        if entry is not None:
            #print(entry)
            if 'textSectionList' in entry:
                textSectionList = entry['textSectionList']
                for ts in textSectionList:
                    if ts['textSection']['textSectionName'] == 'description':
                        d = ts['textSection']['textSectionContent']
                        #there are internal references to OMIM identifiers in the description, I am
                        #formatting them in our style.
                        d = re.sub('{(\\d+)}','OMIM:\\1',d)
                        break
        return d

    def _get_omim_ids(self):
        omimids = []

        #an omim-specific thing here; from the omim.txt.gz file, get the omim numbers
        #not unzipping the file
        print("INFO: Obtaining OMIM record identifiers")
        line_counter=0
        omimfile=('/').join((self.rawdir,self.files['all']['file']))
        with gzip.open(omimfile, 'rb') as f:
            for line in f:
                line=line.decode().strip()
                if (line=="*FIELD* NO"):
                   line_counter += 1
                   #read the next line
                   number=f.readline().decode().strip()
                   omimids.append(number)

        print("INFO: Done.  I found",omimids.__len__(),"omim ids")
        return omimids