
import re
import hashlib
import urllib
import os
import time
import logging
from datetime import datetime
from stat import ST_CTIME, ST_SIZE
from rdflib import ConjunctiveGraph, Namespace
from rdflib.namespace import FOAF, DC, RDFS, OWL
from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils

__author__ = 'nicole'

logger = logging.getLogger(__name__)
core_bindings = {'dc': DC, 'foaf': FOAF, 'rdfs': RDFS}


class Source:
    """
    Abstract class for any data sources that we'll import and process.
    Each of the subclasses will fetch() the data, scrub() it as necessary,
    then parse() it into a graph.  The graph will then be written out to
    a single self.name().ttl file.
    """

    namespaces = {}
    files = {}

    def __init__(self, name=None):
        if name is not None:
            logger.info("Processing Source \"%s\"", name)
        self.testOnly = False
        self.name = name
        self.path = ""
        self.graph = ConjunctiveGraph()
        # to be used to store a subset of data for testing downstream.
        self.testgraph = ConjunctiveGraph()
        self.triple_count = 0
        self.outdir = 'out'
        self.testdir = 'tests'
        self.rawdir = 'raw'
        self.dataset = None
        # set to True if you want to materialze identifiers for BNodes
        self.nobnodes = False
        if self.name is not None:
            self.rawdir = '/'.join((self.rawdir, self.name))
            self.outfile = '/'.join((self.outdir, self.name + ".ttl"))
            logger.info("Setting outfile to %s", self.outfile)

            self.datasetfile = '/'.join((self.outdir,
                                         self.name + '_dataset.ttl'))
            logger.info("Setting dataset file to %s", self.datasetfile)

            self.testfile = '/'.join((self.outdir, self.name + "_test.ttl"))
            logger.info("Setting testfile to %s", self.testfile)

        # if raw data dir doesn't exist, create it
        if not os.path.exists(self.rawdir):
            os.makedirs(self.rawdir)
            p = os.path.abspath(self.rawdir)
            logger.info("creating raw directory for %s at %s", self.name, p)

        # if output dir doesn't exist, create it
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            p = os.path.abspath(self.outdir)
            logger.info("created output directory %s", p)

        # will be set to True if the intention is
        # to only process and write the test data
        self.testOnly = False
        self.testMode = False

        for g in [self.graph, self.testgraph]:
            self.declareAsOntology(g)

        return

    def load_core_bindings(self):

        for g in [self.graph, self.testgraph]:
            g.bind("dc", DC)
            g.bind("foaf", FOAF)
            g.bind("rdfs", RDFS)
            g.bind('owl', OWL)

        return

    def load_bindings(self):
        self.load_core_bindings()
        for g in [self.graph, self.testgraph]:

            for k in self.namespaces.keys():
                v = self.namespaces[k]
                g.bind(k, Namespace(v))

            for k in curie_map.get().keys():
                v = curie_map.get()[k]
                g.bind(k, Namespace(v))
        return

    def fetch(self, is_dl_forced=False):
        """
        abstract method to fetch all data from an external resource.
        this should be overridden by subclasses
        :return: None

        """
        return

    def parse(self):
        """
        abstract method to parse all data from an external resource,
        that was fetched in fetch() this should be overridden by subclasses
        :return: None

        """
        return

    def write(self, format='rdfxml', stream=None):
        """
        This convenience method will write out all of the graphs
        associated with the source.
        Right now these are hardcoded to be a single "graph" and a "dataset".
        If you do not supply stream='stdout'
        it will default write these to files.

        In addition, if the version number isn't yet set in the dataset,
        it will be set to the date on file.
        :return: None

        """
        format_to_xtn = {
            'rdfxml': 'xml', 'turtle': 'ttl'
        }

        # make the regular graph output file
        file = None
        if self.name is not None:
            file = '/'.join((self.outdir, self.name))
            if format in format_to_xtn:
                file = '.'.join((file, format_to_xtn.get(format)))
            else:
                file = '.'.join((file, format))
            # make the datasetfile name
            datasetfile = '/'.join((self.outdir, self.name+'_dataset'))
            if format in format_to_xtn:
                datasetfile = '.'.join((datasetfile, format_to_xtn.get(format)))
            else:
                datasetfile = '.'.join((datasetfile, format))

            logger.info("No version set for this datasource; setting to date issued.")

            if self.dataset is not None and self.dataset.version is None:
                self.dataset.set_version_by_date()
        else:
            logger.warning("No output file set. Using stdout")
            stream = 'stdout'

        # start off with only the dataset descriptions
        graphs = [
            {'g': self.dataset.getGraph(), 'file': datasetfile},
        ]

        # add the other graphs to the set to write, if not in the test mode
        if self.testMode:
            graphs += [{'g': self.testgraph, 'file': self.testfile}]
        else:
            graphs += [{'g': self.graph, 'file': file}]

        gu = GraphUtils(None)
        # loop through each of the graphs and print them out

        for g in graphs:
            f = None
            if stream is None:
                f = g['file']
            elif stream.lower().strip() == 'stdout':
                f = None
            else:
                logger.error("I don't understand your stream.")
                return
            if format == 'raw':
                gu.write_raw_triples(g['g'], file=f)
            else:
                gu.write(g['g'], format, file=f)

        return

    def whoami(self):
        logger.info("I am %s", self.name)
        return

    def make_id(self, long_string):
        """
        a method to create unique identifiers based on very long strings
        currently implemented with md5
        :param long_string:
        :return:
        """
        # FIXME for now, this will do md5.  probably not the best long-term solution
        # note others available: md5(), sha1(), sha224(), sha256(), sha384(), and sha512()

        byte_string = long_string.encode("utf-8")

        return ':'.join(('MONARCH', hashlib.md5(byte_string).hexdigest()))


    def checkIfRemoteIsNewer(self, remote, local, headers):
        """
        Given a remote file location, and the corresponding local file this will
        check the datetime stamp on the files to see if the remote one is newer.
        This is a convenience method to be used so that we don't have to
        re-fetch files that we already have saved locally
        :param remote: URL of file to fetch from remote server
        :param local: pathname to save file to locally
        :return: True if the remote file is newer and should be downloaded
        """
        logger.info("Checking if remote file (%s) is newer than local (%s)...", remote, local)

        # check if local file exists
        # if no local file, then remote is newer
        if not os.path.exists(local):
            logger.info("File does not exist locally")
            return True
        # get remote file details
        if headers is not None:
            req = urllib.request.Request(remote, headers=headers)
        else:
            req = urllib.request.Request(remote)

        logger.info("Request header: %s", str(req.header_items()))
        response = urllib.request.urlopen(req)

        try:
            resp_header = response.getheaders()
            size = resp_header.get('Content-length')
            last_modified = resp_header.get('last-modified')  # check me
        except OSError as e:  #URLError?
            resp_header = None
            size = 0
            last_modified = None
            logger.error(e)

        st = os.stat(local)
        logger.info("Local file date: %s",
                    datetime.utcfromtimestamp(st[ST_CTIME]))

        if last_modified is not None:
            # Thu, 07 Aug 2008 16:20:19 GMT
            dt_obj = datetime.strptime(last_modified,
                                       "%a, %d %b %Y %H:%M:%S %Z")
            # get local file details

            # check date on local vs remote file
            if dt_obj > datetime.utcfromtimestamp(st[ST_CTIME]):
                # check if file size is different
                if st[ST_SIZE] != size:
                    logger.info("Newer file exists on remote server")
                    return True
                else:
                    logger.info("Remote file has same filesize--will not download")
        elif st[ST_SIZE] != size:
            logger.info("Object on server is difference size in comparison to local file")
            return True

        return False

    def get_files(self, is_dl_forced):
        """
        Given a set of files for this source, it will go fetch them, and add
        set a default version by date.  If you need to set the version number
        by another method, then it can be set again.
        :return:
        """

        st = None
        for f in self.files.keys():
            logger.info("Getting %s", f)
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                '/'.join((self.rawdir, file['file'])),
                                is_dl_forced, file.get('headers'))
            self.dataset.setFileAccessUrl(file['url'])

            st = os.stat('/'.join((self.rawdir, file['file'])))

        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        # FIXME change this so the date is attached only to each file,
        # not the entire dataset
        self.dataset.set_date_issued(filedate)

        return

    def fetch_from_url(self, remotefile, localfile, is_dl_forced, headers=None):
        """
        Given a remote url and a local filename, this will first verify
        if the remote file is newer; if it is,
        this will pull the remote file and save it to the specified localfile,
        reporting the basic file information once it is downloaded
        :param remotefile: URL of remote file to fetch
        :param localfile: pathname of file to save locally
        :return: None

        """

        CHUNK = 16 * 1024
        if ((is_dl_forced is True) or
                (self.checkIfRemoteIsNewer(remotefile, localfile, headers))):
            logger.info("Fetching from %s", remotefile)
            # TODO url verification, etc
            if headers is not None:
                r = urllib.request.Request(remotefile, headers=headers)
            else:
                r = urllib.request.Request(remotefile)

            response = urllib.request.urlopen(r)

            with open(localfile, 'wb') as fd:
                while True:
                    chunk = response.read(CHUNK)
                    if not chunk:
                        break
                    fd.write(chunk)

            logger.info("Finished.  Wrote file to %s", localfile)
            if self.compare_local_remote_bytes(remotefile, localfile):
                logger.debug("local file is same size as remote after download")
            else:
                raise Exception("Error when downloading files: "+\
                                "local file size does not match"+\
                                " remote file size")
        else:
            logger.info("Using existing file %s", localfile)

        st = os.stat(localfile)
        logger.info("file size: %s", st[ST_SIZE])
        logger.info("file created: %s",
                    time.asctime(time.localtime(st[ST_CTIME])))
        return

    def process_xml_table(self, elem, table_name, processing_function, limit):
        """
        This is a convenience function to process the elements of
        an xml document, when the xml is used as an alternative way
        of distributing sql-like tables.  In this case, the "elem" is akin to an
        sql table, with it's name of ```table_name```.
        It will then process each ```row```
            given the ```processing_function``` supplied.

        :param elem: The element data
        :param table_name: The name of the table to process
        :param processing_function: The row processing function
        :param limit:
        :return:

        """

        line_counter = 0
        table_data = elem.find("[@name='"+table_name+"']")
        if table_data is not None:
            logger.info("Processing "+table_name)
            row = {}
            for r in table_data.findall('row'):
                for f in r.findall('field'):
                    ats = f.attrib
                    row[ats['name']] = f.text
                processing_function(row)
                line_counter += 1
                if self.testMode and limit is not None and line_counter > limit:
                    continue

            elem.clear()  # discard the element

        return

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

        return

    def get_file_md5(self, directory, file, blocksize=2**20):
        # reference: http://stackoverflow.com/questions/
        #            1131220/get-md5-hash-of-big-files-in-python

        md5 = hashlib.md5()
        with open(os.path.join(directory, file), "rb") as f:
            while True:
                buffer = f.read(blocksize)
                if not buffer:
                    break
                md5.update(buffer)

        return md5.hexdigest()

    def get_remote_content_len(self, remote, headers=None):
        """
        :param remote:
        :return: size of remote file
        """

        if headers is not None:
            req = urllib.request.Request(remote, headers=headers)
        else:
            req = urllib.request.Request(remote)

        try:
            response = urllib.request.urlopen(req)

            resp_header = response.getheaders()
            byte_size = resp_header.get('Content-length')
        except OSError as e:
            byte_size = None
            logger.error(e)

        return byte_size

    def get_local_file_size(self, localfile):
        """
        :param localfile:
        :return: size of file
        """
        byte_size = os.stat(localfile)
        return byte_size[ST_SIZE]

    def compare_local_remote_bytes(self, remotefile, localfile):
        """
        test to see if fetched file is the same size as the remote file
        using information in the content-length field in the HTTP header
        :return: True or False
        """
        is_equal = True
        remote_size = self.get_remote_content_len(remotefile)
        local_size = self.get_local_file_size(localfile)
        if remote_size is not None and local_size != int(remote_size):
            is_equal = False
            logger.error('local file and remote file different sizes\n'
                         '%s has size %s, %s has size %s', localfile,
                         local_size, remotefile, remote_size)
        return is_equal

    def file_len(self, fname):
        with open(fname) as f:
            l = sum(1 for line in f)
        return l

    def settestonly(self, testonly):
        """
        Set that this source should only be processed in testMode
        :param testOnly:
        :return: None
        """

        self.testOnly = testonly

        return

    def settestmode(self, mode):
        """
        Set testMode to (mode).
            True: run the Source in testMode;
            False: run it in full mode
        :param mode:
        :return: None
        """

        self.testMode = mode

        return

    def getTestSuite(self):
        """
        An abstract method that should be overwritten with
        tests appropriate for the specific source.
        :return:
        """
        return None

    def setnobnodes(self, materialize_bnodes):
        """
        If materialze_bnodes is True,
        then all usages of BNodes will be materialized
        by putting the BNodes into the BASE space,
        and prefixing the numeric portion (after the colon) with an underscore.
        :param materialize_bnodes:
        :return:

        """

        self.nobnodes = materialize_bnodes

        return

    def declareAsOntology(self, graph):
        """
        The file we output needs to be declared as an ontology,
        including it's version information.
        Further information will be augmented in the dataset object.
        :param version:
        :return:
        """
        # <http://data.monarchinitiative.org/ttl/biogrid.ttl> a owl:Ontology ;
        # owl:versionInfo <http://archive.monarchinitiative.org/ttl/biogrid-YYYY-MM-DD.ttl>

        gu = GraphUtils(curie_map.get())

        ontology_file_id = 'MonarchData:'+self.name+".ttl"
        gu.addOntologyDeclaration(graph, ontology_file_id)

        # add timestamp as version info

        t = datetime.now()
        t_string = t.strftime("%Y-%m-%d-%H-%M")
        ontology_version = self.name+'-'+t_string
        archive_url = 'MonarchArchive:'+ontology_version+'.ttl'
        gu.addOWLVersionIRI(graph, ontology_file_id, archive_url)
        gu.addOWLVersionInfo(graph, ontology_file_id, ontology_version)

        # TODO make sure this is synced with the Dataset class

        return

    @staticmethod
    def remove_backslash_r(filename, encoding):
        """
        A helpful utility to remove '\r' from any file.
        This will read a file into memory,
        and overwrite the contents of the original file.
        :param filename:
        :return:
        """

        f = open(filename, 'r', encoding=encoding, newline='\n')
        contents = f.read()
        f.close()
        contents = re.sub(r'\r', '', contents)
        with open(filename, "w") as f:
            f.truncate()
            f.write(contents)

        return
