__author__ = 'nicole'

from rdflib import BNode, ConjunctiveGraph, Literal, RDF, Graph
from rdflib.namespace import FOAF, DC, RDFS, OWL

import urllib, csv, os, time
from urllib import request
from datetime import datetime
from stat import *
import hashlib
import subprocess
from subprocess import check_call


core_bindings = {'dc': DC, 'foaf': FOAF, 'rdfs': RDFS}



class Source:
    '''
    Abstract class for any data sources that we'll import and process.
    Each of the subclasses will fetch() the data, scrub() it as necessary,
    then parse() it into a graph.  The graph will then be written out to
    a single self.name().ttl file.
    '''
    def __init__(self,args=[]):
        self.name = args
        self.path = ""
        self.graph = ConjunctiveGraph()
        self.triple_count = 0
        self.outdir = 'out'
        self.rawdir = 'raw'

        self.outfile = ('/').join((self.outdir,self.name + ".ttl"))
        print("Setting outfile to", self.outfile)

        self.datasetfile = ('/').join((self.outdir,self.name + '_dataset.ttl'))
        print("Setting dataset file to", self.datasetfile)

        return

    def load_core_bindings(self):
        self.graph.bind("dc", DC)
        self.graph.bind("foaf", FOAF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind('owl', OWL)

        return

    def load_bindings(self):
        #should each source provide it's own binding dictionary?
        #or should it be a universal one, like the current context file?
        return

    def fetch(self):
        '''
        abstract method to fetch all data from an external resource.
        this should be overridden by subclasses
        :return: None
        '''
        return

    def parse(self):
        '''
        abstract method to parse all data from an external resource, that was fetched in
        fetch()
        this should be overridden by subclasses
        :return: None
        '''
        return

    def write(self, file=None):
        '''
        a basic graph writer (to stdout) for any of the sources, writing raw triples
        :return:
        '''
        if (file is not None):
            filewriter = open(file, 'w')
            print("Writing raw triples to ",file)
            print(self.graph.serialize(format="rdfxml").decode(), file=filewriter)
            filewriter.close()
        else:
            print(self.graph.serialize(format="rdfxml").decode())
        return

    def whoami(self):
        print("I am ",self.name)
        return

    def make_id(self,long_string):
        '''
        a method to create unique identifiers based on very long strings
        currently implemented with md5
        :param long_string:
        :return:
        '''
        #FIXME for now, this will do md5.  probably not the best long-term solution
        #note others available: md5(), sha1(), sha224(), sha256(), sha384(), and sha512()
        byte_string = long_string.encode("utf-8")
        return hashlib.md5(byte_string).hexdigest()


    def checkIfRemoteIsNewer(self,remote,local):
        '''
        Given a remote file location, and the corresponding local file
        this will check the datetime stamp on the files to see if the remote one
        is newer.  This is a convenience method to be used so that we don't have to
        re-fetch files that we already have saved locally
        :param remote: URL of file to fetch from remote server
        :param local: pathname to save file to locally
        :return: True if the remote file is newer and should be downloaded
        '''
        print("Checking if remote file is newer...")
        #check if local file exists
        #if no local file, then remote is newer
        if (not os.path.exists(local)):
            print("File does not exist locally")
            return True
        #get remote file details
        d = urllib.request.urlopen(remote)
        size=d.info()['Content-Length']

        st = os.stat(local)
        print("Local file date:",datetime.utcfromtimestamp(st[ST_CTIME]))

        last_modified = d.info()['Last-Modified']
        if (last_modified is not None):
            #Thu, 07 Aug 2008 16:20:19 GMT
            dt_obj = datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            #get local file details

            #check date on local vs remote file
            if (dt_obj > datetime.utcfromtimestamp(st[ST_CTIME])):
                #check if file size is different
                if (st[ST_SIZE] != size):
                    print("Newer file exists on remote server")
                    return True
                else:
                    print("Remote file has same filesize--will not download")
        elif (st[ST_SIZE] != size):
            print("Object on server is same size as local file; assuming unchanged")
        return False

    def fetch_from_url(self,remotefile,localfile):
        '''
        Given a remote url and a local filename, this will first verify
        if the remote file is newer; if it is, this will pull the remote file
        and save it to the specified localfile, reporting the basic file information
        once it is downloaded
        :param remotefile: URL of remote file to fetch
        :param localfile: pathname of file to save locally
        :return: None
        '''
        if (self.checkIfRemoteIsNewer(remotefile, localfile)):
            print("Fetching from ", remotefile)
            # TODO url verification, etc
            annotation_file = urllib.request
            annotation_file.urlretrieve(remotefile, localfile)
            print("Finished.  Wrote file to", localfile)
        else:
            print("Using existing file.")

        st = os.stat(localfile)
        print("file size:", st[ST_SIZE])
        print("file created:", time.asctime(time.localtime(st[ST_CTIME])))
        return


    def verify(self):
        '''
        abstract method to verify the integrity of the data fetched and turned into triples
        this should be overridden by tests in subclasses
        :return: True if all tests pass
        '''
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        return status

    def _verify(self,f):
        '''
         a simple test to see if the turtle file can be loaded into a graph and parsed
        :param f: file of ttl
        :return:  Boolean status (passed: True; did not pass: False)
        '''
        vg = Graph()
        vg.parse(f, format="turtle")
        print ('Graph nodes:',len(vg))
        return

    def _verifyowl(self,f):
        '''
        test if the ttl can be parsed by owlparser
        this expects owltools to be accessible from commandline
        :param f: file of ttl
        :return: Boolean status (passed: True; did not pass: False)
        '''
        status = check_call(["owltools", f],stderr=subprocess.STDOUT)
        return status

