import psycopg2

__author__ = 'nicole'

from rdflib import BNode, ConjunctiveGraph, Literal, RDF, Graph, Namespace
from rdflib.namespace import FOAF, DC, RDFS, OWL

import urllib, csv, os, time, logging
from urllib import request
from datetime import datetime
from stat import *
import hashlib
import subprocess
from subprocess import check_call
import curie_map

from utils.GraphUtils import GraphUtils

logger = logging.getLogger(__name__)
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
            logger.info("creating raw directory for resource %s"+self.name)
            os.makedirs(self.rawdir)

        self.outfile = ('/').join((self.outdir,self.name + ".ttl"))
        logger.info("Setting outfile to %s", self.outfile)

        self.datasetfile = ('/').join((self.outdir,self.name + '_dataset.ttl'))
        logger.info("Setting dataset file to %s", self.datasetfile)

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
                logger.error("I don't understand your stream.")
        return

    def whoami(self):
        logger.info("I am %s", self.name)
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
        logger.info("Checking if remote file is newer...")
        #check if local file exists
        #if no local file, then remote is newer
        if (not os.path.exists(local)):
            logger.info("File does not exist locally")
            return True
        #get remote file details
        d = urllib.request.urlopen(remote)
        size=d.info()['Content-Length']

        st = os.stat(local)
        logger.info("Local file date: %s",datetime.utcfromtimestamp(st[ST_CTIME]))

        last_modified = d.info()['Last-Modified']
        if (last_modified is not None):
            #Thu, 07 Aug 2008 16:20:19 GMT
            dt_obj = datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            #get local file details

            #check date on local vs remote file
            if (dt_obj > datetime.utcfromtimestamp(st[ST_CTIME])):
                #check if file size is different
                if (st[ST_SIZE] != size):
                    logger.info("Newer file exists on remote server")
                    return True
                else:
                    logger.info("Remote file has same filesize--will not download")
        elif (st[ST_SIZE] != size):
            logger.info("Object on server is same size as local file; assuming unchanged")
        return False

    def get_files(self, is_dl_forced):
        '''
        Given a set of files for this source, it will go fetch them, and add
        set a default version by date.  If you need to set the version number
        by another method, then it can be set again.
        :return:
        '''

        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                ('/').join((self.rawdir,file['file'])),
                                is_dl_forced)
            self.dataset.setFileAccessUrl(file['url'])

            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        self.dataset.setVersion(filedate)

        return

    def fetch_from_url(self, remotefile, localfile, is_dl_forced):
        '''
        Given a remote url and a local filename, this will first verify
        if the remote file is newer; if it is, this will pull the remote file
        and save it to the specified localfile, reporting the basic file information
        once it is downloaded
        :param remotefile: URL of remote file to fetch
        :param localfile: pathname of file to save locally
        :return: None
        '''
        if ((is_dl_forced is True) or
           (self.checkIfRemoteIsNewer(remotefile, localfile))):
            logger.info("Fetching from ", remotefile)
            # TODO url verification, etc
            annotation_file = urllib.request
            annotation_file.urlretrieve(remotefile, localfile)
            logger.info("Finished.  Wrote file to %s", localfile)
        else:
            logger.info("Using existing file %s", localfile)

        st = os.stat(localfile)
        logger.info("file size: %s", st[ST_SIZE])
        logger.info("file created: %s", time.asctime(time.localtime(st[ST_CTIME])))
        return

    def fetch_from_pgdb(self,tables,cxn,limit=None):
        '''
        Will fetch all Postgres tables from the specified database in the cxn connection parameters.
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
                logger.info("Fetching data from table %s",t)
                self._getcols(cur,t)
                query=(' ').join(("SELECT * FROM",t))
                countquery=(' ').join(("SELECT COUNT(*) FROM",t))
                if (limit is not None):
                    query=(' ').join((query,"LIMIT",str(limit)))
                    countquery=(' ').join((countquery,"LIMIT",str(limit)))

                outfile=('/').join((self.rawdir,t))

                #check local copy.  assume that if the # rows are the same, that the table is the same
                #TODO may want to fix this assumption
                filerowcount=-1
                if os.path.exists(outfile):
                    #get rows in the file
                    filerowcount=self.file_len(outfile)
                    logger.info("rows in local file: %s",filerowcount)

                #get rows in the table
                #tablerowcount=cur.rowcount
                cur.execute(countquery)
                tablerowcount=cur.fetchone()[0]
                if (filerowcount < 0 or (filerowcount-1) != tablerowcount):  #rowcount-1 because there's a header
                    logger.info("local (%s) different from remote (%s); fetching.", filerowcount, tablerowcount)
                    #download the file
                    logger.info("COMMAND:%s", query)
                    outputquery = "COPY ({0}) TO STDOUT WITH DELIMITER AS '\t' CSV HEADER".format(query)
                    with open(outfile, 'w') as f:
                        cur.copy_expert(outputquery, f)
                else:
                    logger.info("local data same as remote; reusing.")

        finally:
            if con:
                con.close()
        return

    def fetch_query_from_pgdb(self,qname,query,con,cxn,limit=None):
        """
        Supply either an already established connection, or connection parameters.
        The supplied connection will override any separate cxn parameter
        :param qname:  The name of the query to save the output to
        :param query:  The SQL query itself
        :param con:  The already-established connection
        :param cxn: The postgres connection information
        :param limit: If you only want a subset of rows from the query
        :return:
        """
        if (con is None and cxn is None):
            print("ERROR: you need to supply connection information")
            return
        if (con is None and cxn is not None):
            con = psycopg2.connect(host=cxn['host'],database=cxn['database'],port=cxn['port'], user=cxn['user'], password=cxn['password'])

        outfile=('/').join((self.rawdir,qname))
        cur = con.cursor()
        countquery=(' ').join(("SELECT COUNT(*) FROM (",query,") x"))  #wrap the query to get the count
        if (limit is not None):
            countquery=(' ').join((countquery,"LIMIT",str(limit)))

        #check local copy.  assume that if the # rows are the same, that the table is the same
        filerowcount=-1
        if os.path.exists(outfile):
            #get rows in the file
            filerowcount=self.file_len(outfile)
            print("INFO: rows in local file: ",filerowcount)

        #get rows in the table
        #tablerowcount=cur.rowcount
        cur.execute(countquery)
        tablerowcount=cur.fetchone()[0]

        if (filerowcount < 0 or (filerowcount-1) != tablerowcount):  #rowcount-1 because there's a header
            print("INFO: local (",filerowcount,") different from remote (",tablerowcount,"); fetching.")
            #download the file
            print("COMMAND:",query)
            outputquery = "COPY ({0}) TO STDOUT WITH DELIMITER AS '\t' CSV HEADER".format(query)
            with open(outfile, 'w') as f:
                cur.copy_expert(outputquery, f)
        else:
            print("INFO: local data same as remote; reusing.")

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
        logger.info('Found %s graph nodes',len(vg))
        if len(vg) > 0:
            return True
        return False

    def _verifyowl(self,f):
        '''
        test if the ttl can be parsed by owlparser
        this expects owltools to be accessible from commandline
        :param f: file of ttl
        :return: Boolean status (passed: True; did not pass: False)
        '''
        status = check_call(["owltools", f],stderr=subprocess.STDOUT)
        #returns zero is success!
        if (status != 0):
            logger.error('finished verifying with owltools with status %s',status)
            return False
        else:
            return True

    def _check_list_len(self, row, length):
        """
        Sanity check for csv parser
        :param row
        :param length
        :return:None
        """
        if len(row) != length:
            raise Exception("row length does not match expected length of " +
                            str(length)+"\nrow: "+str(row))