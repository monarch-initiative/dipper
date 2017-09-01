import re
import hashlib
import os
import time
import logging
import urllib
import csv
import yaml
from datetime import datetime
from stat import ST_CTIME, ST_SIZE
from dipper.graph.RDFGraph import RDFGraph
from dipper.graph.StreamedGraph import StreamedGraph
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Model import Model

logger = logging.getLogger(__name__)
CHUNK = 16 * 1024


class Source:
    """
    Abstract class for any data sources that we'll import and process.
    Each of the subclasses will fetch() the data, scrub() it as necessary,
    then parse() it into a graph.  The graph will then be written out to
    a single self.name().ttl file.
    """

    namespaces = {}
    files = {}

    def __init__(self, graph_type, are_bnodes_skized=False, name=None):

        self.graph_type = graph_type
        self.are_bnodes_skized = are_bnodes_skized

        if name is not None:
            logger.info("Processing Source \"%s\"", name)
        self.testOnly = False
        self.name = name
        self.path = ""
        # to be used to store a subset of data for testing downstream.
        self.triple_count = 0
        self.outdir = 'out'
        self.testdir = 'tests'
        self.rawdir = 'raw'
        self.dataset = None
        # set to True if you want to materialze identifiers for BNodes
        if self.name is not None:
            self.rawdir = '/'.join((self.rawdir, self.name))
            self.outfile = '/'.join((self.outdir, self.name + ".ttl"))
            logger.info("Setting outfile to %s", self.outfile)

            self.datasetfile = '/'.join(
                (self.outdir, self.name + '_dataset.ttl'))
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

        if graph_type == 'rdf_graph':
            self.graph = RDFGraph(are_bnodes_skized)
            self.testgraph = RDFGraph(True)
        elif graph_type == 'streamed_graph':
            source_file = open(self.outfile.replace(".ttl", ".nt"), 'w')
            test_file = open(self.testfile.replace(".ttl", ".nt"), 'w')
            self.graph = StreamedGraph(are_bnodes_skized, source_file)
            self.testgraph = StreamedGraph(are_bnodes_skized, test_file)
        else:
            logger.error("{} graph type not supported\n"
                         "valid types: rdf_graph, streamed_graph"
                         .format(graph_type))

        # will be set to True if the intention is
        # to only process and write the test data
        self.testOnly = False
        self.testMode = False

        for g in [self.graph, self.testgraph]:
            self.declareAsOntology(g)

        return

    def fetch(self, is_dl_forced=False):
        """
        abstract method to fetch all data from an external resource.
        this should be overridden by subclasses
        :return: None

        """
        raise NotImplementedError

    def parse(self, limit):
        """
        abstract method to parse all data from an external resource,
        that was fetched in fetch() this should be overridden by subclasses
        :return: None

        """
        raise NotImplementedError

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
                datasetfile = '.'.join((datasetfile,
                                        format_to_xtn.get(format)))
            else:
                datasetfile = '.'.join((datasetfile, format))

            logger.info(
                "No version set for this datasource; setting to date issued.")

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
            gu.write(g['g'], format, file=f)

        return

    def whoami(self):
        logger.info("I am %s", self.name)
        return

    @staticmethod
    def make_id(long_string, prefix='MONARCH'):
        """
        a method to create DETERMINISTIC identifiers
        based on a string's digest. currently implemented with sha1
        :param long_string:
        :return:

        """
        return ':'.join((prefix, Source.hash_id(long_string)))

    @staticmethod
    def hash_id(long_string):
        """
        prepend 'b' to avoid leading with digit
        truncate to 64bit sized word
        return sha1 hash of string
        :param long_string: str string to be hashed
        :return: str hash of id
        """
        return r'b' + hashlib.sha1(
            long_string.encode('utf-8')).hexdigest()[0:15]

    def checkIfRemoteIsNewer(self, remote, local, headers):
        """
        Given a remote file location, and the corresponding local file
        this will check the datetime stamp on the files to see if the remote
        one is newer.
        This is a convenience method to be used so that we don't have to
        re-fetch files that we already have saved locally
        :param remote: URL of file to fetch from remote server
        :param local: pathname to save file to locally
        :return: True if the remote file is newer and should be downloaded

        """
        logger.info(
            "Checking if remote file (%s) is newer than local (%s)...",
            remote, local)

        # check if local file exists
        # if no local file, then remote is newer
        if os.path.exists(local):
            logger.info("File does exist locally")
        else:
            logger.info("File does not exist locally")
            return True

        # get remote file details
        if headers is not None and headers != []:
            req = urllib.request.Request(remote, headers=headers)
        else:
            req = urllib.request.Request(remote)

        logger.info("Request header: %s", str(req.header_items()))

        response = urllib.request.urlopen(req)

        try:
            resp_headers = response.info()
            size = resp_headers.get('Content-length')
            last_modified = resp_headers.get('last-modified')  # check me
        except OSError as e:  # URLError?
            resp_headers = None
            size = 0
            last_modified = None
            logger.error(e)

        st = os.stat(local)
        logger.info(
            "Local file date: %s",
            datetime.utcfromtimestamp(st[ST_CTIME]))

        if last_modified is not None:
            # Thu, 07 Aug 2008 16:20:19 GMT
            dt_obj = datetime.strptime(
                last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            # get local file details

            # check date on local vs remote file
            if dt_obj > datetime.utcfromtimestamp(st[ST_CTIME]):
                # check if file size is different
                if st[ST_SIZE] != size:
                    logger.info("Newer file exists on remote server")
                    return True
                else:
                    logger.info(
                        "Remote file has same filesize--will not download")
        elif st[ST_SIZE] != size:
            logger.info(
                "Object on server is difference size to local file")
            return True

        return False

    def get_files(self, is_dl_forced, files=None):
        """
        Given a set of files for this source, it will go fetch them, and
        set a default version by date.  If you need to set the version number
        by another method, then it can be set again.
        :param is_dl_forced - boolean
        :param files dict - override instance files dict
        :return: None
        """

        st = None
        if files is None:
            files = self.files
        for fname in files.keys():
            logger.info("Getting %s", fname)
            filesource = files.get(fname)
            self.fetch_from_url(
                filesource['url'], '/'.join((self.rawdir, filesource['file'])),
                is_dl_forced, filesource.get('headers'))
            self.dataset.setFileAccessUrl(filesource['url'])

            st = os.stat('/'.join((self.rawdir, filesource['file'])))

        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        # FIXME change this so the date is attached only to each file,
        # not the entire dataset
        self.dataset.set_date_issued(filedate)

        return

    def fetch_from_url(
            self, remotefile, localfile=None, is_dl_forced=False,
            headers=None):
        """
        Given a remote url and a local filename, this will first verify
        if the remote file is newer; if it is,
        this will pull the remote file and save it to the specified localfile,
        reporting the basic file information once it is downloaded
        :param remotefile: URL of remote file to fetch
        :param localfile: pathname of file to save locally
        :return: None

        """
        response = None
        if ((is_dl_forced is True) or localfile is None or
                (self.checkIfRemoteIsNewer(remotefile, localfile, headers))):
            logger.info("Fetching from %s", remotefile)
            # TODO url verification, etc
            if headers is not None:
                request = urllib.request.Request(remotefile, headers=headers)
            else:
                request = urllib.request.Request(remotefile)

            response = urllib.request.urlopen(request)

            if localfile is not None:
                with open(localfile, 'wb') as fd:
                    while True:
                        chunk = response.read(CHUNK)
                        if not chunk:
                            break
                        fd.write(chunk)

                logger.info("Finished.  Wrote file to %s", localfile)
                if self.compare_local_remote_bytes(remotefile, localfile):
                    logger.debug(
                        "local file is same size as remote after download")
                else:
                    raise Exception(
                        "Error when downloading files: local file size " +
                        "does not match remote file size")
                st = os.stat(localfile)
                logger.info("file size: %s", st[ST_SIZE])
                logger.info(
                    "file created: %s",
                    time.asctime(time.localtime(st[ST_CTIME])))
        else:
            logger.info("Using existing file %s", localfile)

        return response

    def process_xml_table(self, elem, table_name, processing_function, limit):
        """
        This is a convenience function to process the elements of an
        xml document, when the xml is used as an alternative way of
        distributing sql-like tables.  In this case, the "elem" is akin to an
        sql table, with it's name of ```table_name```.
        It will then process each ```row``` given the ```processing_function```
        supplied.

        :param elem: The element data
        :param table_name: The name of the table to process
        :param processing_function: The row processing function
        :param limit:

        Appears to be making calls to the elementTree library
        although it not explicitly imported here.

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
                if self.testMode \
                        and limit is not None and line_counter > limit:
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
            raise Exception(
                "row length does not match expected length of " +
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
            resp_header = response.info()
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
            logger.error(
                'local file and remote file different sizes\n'
                '%s has size %s, %s has size %s',
                localfile, local_size, remotefile, remote_size)
        return is_equal

    def file_len(self, fname):
        with open(fname) as f:
            l = sum(1 for line in f)
        return l

    @staticmethod
    def get_eco_map(url):
        """
        To conver the three column file to
        a hashmap we join primary and secondary keys,
        for example
        IEA	GO_REF:0000002	ECO:0000256
        IEA	GO_REF:0000003	ECO:0000501
        IEA	Default	ECO:0000501

        becomes
        IEA-GO_REF:0000002: ECO:0000256
        IEA-GO_REF:0000003: ECO:0000501
        IEA: ECO:0000501

        :return: dict
        """
        eco_map = {}
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)

        for line in response:
            line = line.decode('utf-8').rstrip()
            if re.match(r'^#', line):
                continue
            (code, go_ref, eco_curie) = line.split('\t')
            if go_ref != 'Default':
                eco_map["{}-{}".format(code, go_ref)] = eco_curie
            else:
                eco_map[code] = eco_curie

        return eco_map

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
        - True: run the Source in testMode;
        - False: run it in full mode
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

    def declareAsOntology(self, graph):
        """
        The file we output needs to be declared as an ontology,
        including it's version information.

        TEC: I am not convinced dipper reformating external data as RDF triples
        makes an OWL ontology (nor that it should be considered a goal).

        Proper ontologies are built by ontologists. Dipper reformats data
        and anotates/decorates it with a minimal set of carefully arranged
        terms drawn from from multiple proper ontologies.
        Which allows the whole (dipper's RDF triples and parent ontologies)
        to function as a single ontology we can reason over when combined
        in a store such as SciGraph.

        Including more than the minimal ontological terms in dipper's RDF
        output constitutes a liability as it allows divergence.

        Further information will be augmented in the dataset object.
        :param version:
        :return:

        """

        # <http://data.monarchinitiative.org/ttl/biogrid.ttl> a owl:Ontology ;
        # owl:versionInfo
        # <https://archive.monarchinitiative.org/YYYYMM/ttl/biogrid.ttl>

        model = Model(graph)

        ontology_file_id = 'MonarchData:' + self.name + ".ttl"
        model.addOntologyDeclaration(ontology_file_id)

        # add timestamp as version info

        t = datetime.now()
        t_string = t.strftime("%Y-%m-%d")
        ontology_version = t_string
        # TEC this means the MonarchArchive IRI needs the release updated
        # maybe extract the version info from there

        archive_url = 'MonarchArchive:' + 'ttl/' + self.name + '.ttl'
        model.addOWLVersionIRI(ontology_file_id, archive_url)
        model.addOWLVersionInfo(ontology_file_id, ontology_version)

        # TODO make sure this is synced with the Dataset class

        return

    @staticmethod
    def remove_backslash_r(filename, encoding):
        """
        A helpful utility to remove Carriage Return from any file.
        This will read a file into memory,
        and overwrite the contents of the original file.

        TODO: This function may be a liability

        :param filename:

        :return:

        """

        f = open(filename, 'r', encoding=encoding, newline=r'\n')
        contents = f.read()
        f.close()
        contents = re.sub(r'\r', '', contents)
        with open(filename, "w") as f:
            f.truncate()
            f.write(contents)

        return

    @staticmethod
    def open_and_parse_yaml(file):
        """
        :param file: String, path to file containing label-id mappings in
                             the first two columns of each row
        :return: dict where keys are labels and values are ids
        """
        map = dict()
        if os.path.exists(os.path.join(os.path.dirname(__file__), file)):
            map_file = open(os.path.join(os.path.dirname(__file__), file), 'r')
            map = yaml.load(map_file)
            map_file.close()
        else:
            logger.warn("file: {0} not found".format(file))

        return map

    @staticmethod
    def parse_mapping_file(file):
        """
        :param file: String, path to file containing label-id mappings in
                             the first two columns of each row
        :return: dict where keys are labels and values are ids
        """
        id_map = {}
        if os.path.exists(os.path.join(os.path.dirname(__file__), file)):
            with open(
                    os.path.join(os.path.dirname(__file__), file)) as tsvfile:
                reader = csv.reader(tsvfile, delimiter="\t")
                for row in reader:
                    label = row[0]
                    id = row[1]
                    id_map[label] = id

        return id_map
