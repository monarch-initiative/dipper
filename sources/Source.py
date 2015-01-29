import psycopg2

__author__ = 'nicole'

from rdflib import BNode, ConjunctiveGraph, Literal, RDF, Graph, Namespace
from rdflib.namespace import FOAF, DC, RDFS, OWL

import urllib, csv, os, time
from urllib import request
from datetime import datetime
from stat import *
import hashlib
import subprocess
from subprocess import check_call
import curie_map

from utils.GraphUtils import GraphUtils

core_bindings = {'dc': DC, 'foaf': FOAF, 'rdfs': RDFS}



class Source:
    '''
    Abstract class for any data sources that we'll import and process.
    Each of the subclasses will fetch() the data, scrub() it as necessary,
    then parse() it into a graph.  The graph will then be written out to
    a single self.name().ttl file.
    '''

    namespaces = {}
    files = {}

    def __init__(self,args=[]):
        self.name = args
        self.path = ""
        self.graph = ConjunctiveGraph()
        self.triple_count = 0
        self.outdir = 'out'
        self.rawdir = 'raw/'+self.name

        #if raw data dir doesn't exist, create it
        if not os.path.exists(self.rawdir):
            print("INFO: creating raw directory for resource "+self.name)
            os.makedirs(self.rawdir)

        self.outfile = ('/').join((self.outdir,self.name + ".ttl"))
        print("INFO: Setting outfile to", self.outfile)

        self.datasetfile = ('/').join((self.outdir,self.name + '_dataset.ttl'))
        print("INFO: Setting dataset file to", self.datasetfile)

        return

    def load_core_bindings(self):
        self.graph.bind("dc", DC)
        self.graph.bind("foaf", FOAF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind('owl', OWL)

        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.namespaces.keys():
            v = self.namespaces[k]
            self.graph.bind(k, Namespace(v))

        for k in curie_map.get().keys():
            v = curie_map.get()[k]
            self.graph.bind(k, Namespace(v))
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

    def write(self, format='rdfxml', stream=None):
        '''
        This convenience method will write out all of the graphs associated with the source.
        Right now these are hardcoded to be a single "graph" and a "dataset".
        If you do not supply stream='stdout' it will default write these to files
        :return: None
        '''
        format_to_xtn = {
            'rdfxml' : 'xml', 'turtle' : 'ttl'
        }

        #make the regular graph output file
        file=('/').join((self.outdir,self.name))
        if (format in format_to_xtn):
            file=('.').join((file,format_to_xtn.get(format)))
        else:
            file=('.').join((file,format))

        #make the datasetfile name
        datasetfile=('/').join((self.outdir,self.name+'_dataset'))
        if (format in format_to_xtn):
            datasetfile=('.').join((datasetfile,format_to_xtn.get(format)))
        else:
            datasetfile=('.').join((datasetfile,format))

        graphs = [
            {'g':self.graph,'file':file},
            {'g':self.dataset.getGraph(),'file':datasetfile}
        ]
        gu = GraphUtils(None)
        #loop through each of the graphs and print them out
        for g in graphs:
            if (stream is None):
                gu.write(g['g'],format,file=g['file'])
            elif(stream.lowercase().strip() != 'stdout'):
                gu.write(g['g'],format)
            else:
                print("ERROR: I don't understand your stream.")
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
        return (':').join(('MONARCH',hashlib.md5(byte_string).hexdigest()))


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
        print("INFO: Checking if remote file is newer...")
        #check if local file exists
        #if no local file, then remote is newer
        if (not os.path.exists(local)):
            print("INFO: File does not exist locally")
            return True
        #get remote file details
        d = urllib.request.urlopen(remote)
        size=d.info()['Content-Length']

        st = os.stat(local)
        print("INFO: Local file date:",datetime.utcfromtimestamp(st[ST_CTIME]))

        last_modified = d.info()['Last-Modified']
        if (last_modified is not None):
            #Thu, 07 Aug 2008 16:20:19 GMT
            dt_obj = datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            #get local file details

            #check date on local vs remote file
            if (dt_obj > datetime.utcfromtimestamp(st[ST_CTIME])):
                #check if file size is different
                if (st[ST_SIZE] != size):
                    print("INFO: Newer file exists on remote server")
                    return True
                else:
                    print("INFO: Remote file has same filesize--will not download")
        elif (st[ST_SIZE] != size):
            print("INFO: Object on server is same size as local file; assuming unchanged")
        return False

    def get_files(self):
        '''
        Given a set of files for this source, it will go fetch them, and add
        set a default version by date.  If you need to set the version number
        by another method, then it can be set again.
        :return:
        '''

        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],('/').join((self.rawdir,file['file'])))
            self.dataset.setFileAccessUrl(file['url'])

            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        self.dataset.setVersion(filedate)

        return

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
            print("INFO: Fetching from ", remotefile)
            # TODO url verification, etc
            annotation_file = urllib.request
            annotation_file.urlretrieve(remotefile, localfile)
            print("INFO: Finished.  Wrote file to", localfile)
        else:
            print("INFO: Using existing file", localfile)

        st = os.stat(localfile)
        print("INFO: file size:", st[ST_SIZE])
        print("INFO: file created:", time.asctime(time.localtime(st[ST_CTIME])))
        return

    def fetch_files(self, files=dict):
        """
        Given a files dict containing a key pointing to a remote url
        and local filename downloads and stores the file on the local os
        :param files: dict containing the filename and url
            Example: files = {'interactions': {'file': 'filename.csv',
                                               'url': 'http://ctdbase.org/'
                                                      'filename.csv'}
                     }
        :return: date as a string
        """
        for f in files.keys():
            file = files.get(f)
            self.fetch_from_url(file['url'], '/'.join((self.rawdir, file['file'])))
            self.dataset.setFileAccessUrl(file['url'])

            st = os.stat('/'.join((self.rawdir, file['file'])))

        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        return filedate

    def fetch_from_pgdb(self,tables,cxn,limit=None):
        '''
        Will fetch all Postgres tables fom the specified database in the cxn connection parameters.
        This will save them to a local file named the same as the table, in tab-delimited format, including a header.
        :param tables: Names of tables to fetch
        :param cxn: database connection details
        :param limit: A max row count to fetch for each table
        :return: None
        '''
        try:
            con = None
            con = psycopg2.connect(host=cxn['host'],database=cxn['database'],port=cxn['port'], user=cxn['user'], password=cxn['password'])
            cur = con.cursor()
            for t in self.tables:
                print("INFO: Fetching data from table",t)
                self._getcols(cur,t)
                query=(' ').join(("SELECT * FROM",t))
                if (limit is not None):
                    query=(' ').join((query,"LIMIT",str(limit)))
                print("COMMAND:",query)
                cur.execute(query)

                outfile=('/').join((self.rawdir,t))

                #check local copy.  assume that if the # rows are the same, that the table is the same
                #TODO may want to fix this assumption
                filerowcount=-1
                if os.path.exists(outfile):
                    #get rows in the file
                    filerowcount=self.file_len(outfile)
                    print("INFO: rows in local file: ",filerowcount)

                #get rows in the table
                tablerowcount=cur.rowcount
                if (filerowcount < 0 or (filerowcount-1) != tablerowcount):  #rowcount-1 because there's a header
                    print("INFO: local data different from remote; fetching.")
                    #download the file
                    outputquery = "COPY ({0}) TO STDOUT WITH DELIMITER AS '\t' CSV HEADER".format(query)
                    with open(outfile, 'w') as f:
                        cur.copy_expert(outputquery, f)
                else:
                    print("INFO: local data same as remote; reusing.")

        finally:
            if con:
                con.close()
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

