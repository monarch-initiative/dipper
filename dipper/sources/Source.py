import re
import hashlib
import os
import time
import logging
import urllib
import csv
from datetime import datetime
from stat import ST_CTIME, ST_SIZE
from inspect import getdoc
from rdflib import XSD, Literal

import yaml
from dipper.graph.RDFGraph import RDFGraph
from dipper.graph.StreamedGraph import StreamedGraph
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset

LOG = logging.getLogger(__name__)
CHUNK = 16 * 1024  # read remote urls of unknown size in 16k chunks
USER_AGENT = "The Monarch Initiative (https://monarchinitiative.org/; " \
             "info@monarchinitiative.org)"


class Source:
    """
    Abstract class for any data sources that we'll import and process.
    Each of the subclasses will fetch() the data, scrub() it as necessary,
    then parse() it into a graph.  The graph will then be written out to
    a single self.name().<dest_fmt>  file.

    Also provides a means to marshal metadata in a consistent fashion

    Houses the global translation table (from ontology label to ontology term)
    so it may as well be used everywhere.

    """

    namespaces = {}
    files = {}

    def __init__(
            self,
            graph_type='rdf_graph',     # or streamed_graph
            are_bnodes_skized=False,    # typically True
            data_release_version=None,
            name=None,                  # identifier; make an IRI for nquads
            ingest_title=None,
            ingest_url=None,
            ingest_logo=None,     # this should be the name of file on 'MonarchLogoRepo'
            ingest_description=None,
            license_url=None,           # only if it is _our_ lic
            data_rights=None,           # external page that points to their current lic
            file_handle=None
    ):

        # pull in the common test identifiers
        self.all_test_ids = self.open_and_parse_yaml('../../resources/test_ids.yaml')

        self.graph_type = graph_type
        self.are_bnodes_skized = are_bnodes_skized
        self.data_release_version = data_release_version
        self.ingest_title = ingest_title
        self.ingest_url = ingest_url
        self.ingest_logo = ingest_logo
        self.ingest_description = ingest_description
        self.license_url = license_url
        self.data_rights = data_rights
        self.localtt = self.load_local_translationtable(name)

        if name is not None:
            self.name = name.lower()
        elif self.whoami() is not None:
            self.name = self.whoami().lower()

        LOG.info("Processing Source \"%s\"", self.name)
        self.test_only = False
        self.path = ""
        # to be used to store a subset of data for testing downstream.
        self.triple_count = 0
        self.outdir = 'out'
        self.testdir = 'tests'
        self.rawdir = 'raw'
        self.rawdir = '/'.join((self.rawdir, self.name))
        self.testname = name + "_test"
        self.testfile = '/'.join((self.outdir, self.testname + ".ttl"))
        self.datasetfile = None

        # if raw data dir doesn't exist, create it
        if not os.path.exists(self.rawdir):
            os.makedirs(self.rawdir)
            pth = os.path.abspath(self.rawdir)
            LOG.info("creating raw directory for %s at %s", self.name, pth)

        # if output dir doesn't exist, create it
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            pth = os.path.abspath(self.outdir)
            LOG.info("created output directory %s", pth)

        LOG.info("Creating Test graph %s", self.testname)
        # note: tools such as protoge need skolemized blank nodes
        self.testgraph = RDFGraph(True, self.testname)

        if graph_type == 'rdf_graph':
            graph_id = ':MONARCH_' + str(self.name) + "_" + \
                datetime.now().isoformat(' ').split()[0]

            LOG.info("Creating graph  %s", graph_id)
            self.graph = RDFGraph(are_bnodes_skized, graph_id)

        elif graph_type == 'streamed_graph':
            # need to expand on export formats
            dest_file = open(pth + '/' + name + '.nt', 'w')    # where is the close?
            self.graph = StreamedGraph(are_bnodes_skized, dest_file)
            # leave test files as turtle (better human readibility)
        else:
            LOG.error(
                "%s graph type not supported\n"
                "valid types: rdf_graph, streamed_graph", graph_type)

        # pull in global ontology mapping datastructures
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid

        self.curie_map = self.graph.curie_map
        # self.prefix_base = {v: k for k, v in self.curie_map.items()}

        # will be set to True if the intention is
        # to only process and write the test data
        self.test_only = False
        self.test_mode = False

        if self.ingest_description and getdoc(self) is not None:
            self.ingest_description = getdoc(self)

        self.dataset = Dataset(
            identifier=self.curie_map.get("DatasetBase") + self.name,
            data_release_version=self.data_release_version,
            ingest_name=self.name,
            ingest_title=self.ingest_title,
            ingest_url=self.ingest_url,
            ingest_logo=self.ingest_logo,
            ingest_description=self.ingest_description,   # description
            license_url=self.license_url,    # only _OUR_ lic
            data_rights=self.data_rights,    # tries to point to others lics
            graph_type=graph_type,
            file_handle=file_handle
        )

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

    def write(self, fmt='turtle', stream=None, write_metadata_in_main_graph=True):
        """
        This convenience method will write out all of the graphs
        associated with the source.
        Right now these are hardcoded to be a single "graph"
        and a "src_dataset.ttl" and a "src_test.ttl"
        If you do not supply stream='stdout'
        it will default write these to files.

        In addition, if the version number isn't yet set in the dataset,
        it will be set to the date on file.
        :return: None

        """
        fmt_ext = {
            'rdfxml': 'xml',
            'turtle': 'ttl',
            'nt': 'nt',         # ntriples
            'nquads':  'nq',
            'n3': 'n3'          # notation3
        }

        # make the regular graph output file
        dest = None
        if self.name is not None:
            dest = '/'.join((self.outdir, self.name))
            if fmt in fmt_ext:
                dest = '.'.join((dest, fmt_ext.get(fmt)))
            else:
                dest = '.'.join((dest, fmt))
            LOG.info("Setting outfile to %s", dest)

            # make the dataset_file name, always format as turtle
            self.datasetfile = '/'.join(
                (self.outdir, self.name + '_dataset.ttl'))
            LOG.info("Setting dataset file to %s", self.datasetfile)
        else:
            LOG.warning("No output file set. Using stdout")
            stream = 'stdout'

        graph_util = GraphUtils(None)

        # the  _dataset description is always turtle
        graph_util.write(self.dataset.get_graph(), 'turtle', filename=self.datasetfile)

        if self.test_mode:
            # unless we stop hardcoding, the test dataset is always turtle
            LOG.info("Setting testfile to %s", self.testfile)
            graph_util.write(self.testgraph, 'turtle', filename=self.testfile)

        if write_metadata_in_main_graph:
            self.graph = self.graph + self.dataset.get_graph()

        # print graph out
        if stream is None:
            outfile = dest
        elif stream.lower().strip() == 'stdout':
            outfile = None
        else:
            LOG.error("I don't understand our stream.")
            return

        graph_util.write(self.graph, fmt, filename=outfile)

    def whoami(self):
        '''
            pointless convieniance
        '''
        LOG.info("Ingest is %s", self.name)

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
    def hash_id(wordage):  # same as graph/GraphUtils.digest_id(wordage)
        """
        prepend 'b' to avoid leading with digit
        truncate to a 20 char sized word with a leading 'b'
        return truncated sha1 hash of string.

        by the birthday paradox;
            expect 50% chance of collision after 69 billion invocations
            however these are only hoped to be unique within a single file

        Consider reducing to 17 hex chars to fit in a 64 bit word
        16 discounting a leading constant
        gives a 50% chance of collision at about 4.3b billion unique input strings
        (currently _many_ orders of magnitude below that)

        :param long_string: str string to be hashed
        :return: str hash of id
        """
        return 'b' + hashlib.sha1(wordage.encode('utf-8')).hexdigest()[1:20]

    def check_if_remote_is_newer(self, remote, local, headers):
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
        is_remote_newer = False
        LOG.info("Checking if remote file is newer than local \n(%s)", local)

        # check if local file exists
        # if no local file, then remote is newer
        if os.path.exists(local):
            LOG.info("Local File exists as %s", local)
        else:
            LOG.info("Local File does NOT exist as %s", local)
            return True

        # get remote file details
        if headers is None:
            headers = self._get_default_request_headers()

        req = urllib.request.Request(remote, headers=headers)
        LOG.info("Request header: %s", str(req.header_items()))

        response = urllib.request.urlopen(req)

        try:
            resp_headers = response.info()
            size = resp_headers.get('Content-Length')
            last_modified = resp_headers.get('Last-Modified')
        except urllib.error.URLError as err:
            resp_headers = None
            size = 0
            last_modified = None
            LOG.error(err)

        if size is not None and size != '':
            size = int(size)
        else:
            size = 0

        fstat = os.stat(local)
        LOG.info(
            "Local File date: %s",
            datetime.utcfromtimestamp(fstat[ST_CTIME]))

        if last_modified is not None:
            # Thu, 07 Aug 2008 16:20:19 GMT
            dt_obj = datetime.strptime(
                last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            # get local file details

            # check date on local vs remote file
            if dt_obj > datetime.utcfromtimestamp(fstat[ST_CTIME]):
                # check if file size is different
                if fstat[ST_SIZE] < size:
                    LOG.info("New Remote File exists")
                    is_remote_newer = True
                if fstat[ST_SIZE] > size:
                    LOG.warning("New Remote File exists but it is SMALLER")
                    is_remote_newer = True
                # filesize is a fairly imperfect metric here
                LOG.info("New Remote File has same filesize--will not download")
        elif fstat[ST_SIZE] != size:
            LOG.info(
                "Remote File is %i  \t Local File is %i", size, fstat[ST_SIZE])
            is_remote_newer = True

        response.close()

        return is_remote_newer

    def get_files(self, is_dl_forced, files=None, delay=0):
        """
        Given a set of files for this source, it will go fetch them, and
        set a default version by date.  If you need to set the version number
        by another method, then it can be set again.
        :param is_dl_forced - boolean
        :param files dict - override instance files dict
        :return: None
        """

        fstat = None
        if files is None:
            files = self.files
        for fname in files:
            headers = None
            filesource = files[fname]
            if 'headers' in filesource:
                headers = filesource['headers']
            LOG.info("Getting %s", fname)
            # if the key 'clean' exists in the sources `files` dict
            # expose that instead of the longer url
            if 'clean' in filesource and filesource['clean'] is not None:
                self.dataset.set_ingest_source(filesource['clean'])
            else:
                self.dataset.set_ingest_source(filesource['url'])
                LOG.info('Fetching %s in %i seconds', filesource['url'], delay)

            time.sleep(delay)

            self.fetch_from_url(
                filesource['url'], '/'.join((self.rawdir, filesource['file'])),
                is_dl_forced, headers)

            fstat = os.stat('/'.join((self.rawdir, filesource['file'])))
            self.dataset.graph.addTriple(self.dataset.version_level_curie,
                                         self.globaltt["Source (dct)"],
                                         filesource['url']
                                         )
            filedate = Literal(
                datetime.utcfromtimestamp(fstat[ST_CTIME]).strftime("%Y%m%d"),
                datatype=XSD.date)
            self.dataset.graph.addTriple(filesource['url'],
                                         self.globaltt['retrieved_on'],
                                         filedate)

    def fetch_from_url(
            self, remotefile, localfile=None, is_dl_forced=False, headers=None):
        """
        Given a remote url and a local filename, attempt to determine
        if the remote file is newer; if it is,
        fetch the remote file and save it to the specified localfile,
        reporting the basic file information once it is downloaded
        :param remotefile: URL of remote file to fetch
        :param localfile: pathname of file to save locally
        :return: None

        """

        response = None
        if ((is_dl_forced is True) or localfile is None or
                (self.check_if_remote_is_newer(remotefile, localfile, headers))):
            # TODO url verification, etc
            if headers is None:
                headers = self._get_default_request_headers()

            request = urllib.request.Request(remotefile, headers=headers)
            response = urllib.request.urlopen(request)

            if localfile is not None:
                with open(localfile, 'wb') as binwrite:
                    while True:
                        chunk = response.read(CHUNK)
                        if not chunk:
                            break
                        binwrite.write(chunk)

                LOG.info("Finished.  Wrote file to %s", localfile)
                if self.compare_local_remote_bytes(remotefile, localfile, headers):
                    LOG.debug("local file is same size as remote after download")
                else:
                    raise Exception(
                        "Error downloading file: local file size  != remote file size")

                fstat = os.stat(localfile)
                LOG.info("file size: %s", fstat[ST_SIZE])
                LOG.info(
                    "file created: %s", time.asctime(time.localtime(fstat[ST_CTIME])))
                response.close()
            else:
                LOG.error('Local filename is required')
                exit(-1)
        else:
            LOG.info("Using existing file %s", localfile)

    # TODO: rephrase as mysql-dump-xml specific format
    def process_xml_table(self, elem, table_name, processing_function, limit):
        """
        This is a convenience function to process the elements of an xml dump of
        a mysql relational database.
        The "elem" is akin to a mysql table, with it's name of ```table_name```.
        It will process each ```row``` given the ```processing_function``` supplied.
        :param elem: The element data
        :param table_name: The name of the table to process
        :param processing_function: The row processing function
        :param limit:

        Appears to be making calls to the elementTree library
        although it not explicitly imported here.

        :return:

        """

        line_counter = 0
        table_data = elem.find("[@name='" + table_name + "']")
        if table_data is not None:
            LOG.info("Processing " + table_name)
            row = {}
            for line in table_data.findall('row'):
                for field in line.findall('field'):
                    atts = dict(field.attrib)
                    row[atts['name']] = field.text
                processing_function(row)
                line_counter += 1
                if self.test_mode and limit is not None and line_counter > limit:
                    continue

            elem.clear()  # discard the element

    @staticmethod
    def _check_list_len(row, length):
        """
        Sanity check for csv parser
        :param row
        :param length
        :return:None
        """
        if len(row) != length:
            raise Exception(
                "row length does not match expected length of " +
                str(length) + "\nrow: " + str(row))

    @staticmethod
    def get_file_md5(directory, filename, blocksize=2**20):
        # reference:
        # http://stackoverflow.com/questions/1131220/get-md5-hash-of-big-files-in-python

        md5 = hashlib.md5()
        with open(os.path.join(directory, filename), "rb") as bin_reader:
            while True:
                buff = bin_reader.read(blocksize)
                if not buff:
                    break
                md5.update(buff)

        return md5.hexdigest()

    def get_remote_content_len(self, remote, headers=None):
        """
        :param remote:
        :return: size of remote file
        """

        if headers is None:
            headers = self._get_default_request_headers()

        req = urllib.request.Request(remote, headers=headers)

        try:
            response = urllib.request.urlopen(req)
            resp_header = response.info()
            byte_size = resp_header.get('Content-length')
        except OSError as err:
            byte_size = None
            LOG.error(err)

        return byte_size

    @staticmethod
    def get_local_file_size(localfile):
        """
        :param localfile:
        :return: size of file
        """
        byte_size = os.stat(localfile)
        return byte_size[ST_SIZE]

    def compare_local_remote_bytes(self, remotefile, localfile, remote_headers=None):
        """
        test to see if fetched file is the same size as the remote file
        using information in the content-length field in the HTTP header
        :return: True or False
        """
        is_equal = True
        remote_size = self.get_remote_content_len(remotefile, remote_headers)
        local_size = self.get_local_file_size(localfile)
        if remote_size is not None and local_size != int(remote_size):
            is_equal = False
            LOG.error(
                'local file and remote file different sizes\n'
                '%s has size %s, %s has size %s',
                localfile, local_size, remotefile, remote_size)
        return is_equal

    @staticmethod
    def file_len(fname):
        with open(fname) as lines:
            length = sum(1 for line in lines)
        return length

    @staticmethod
    def get_eco_map(url):
        """
        To convert the three column file to
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
        # this would go in a translation table but it is generated dynamically
        # maybe when we move to a make driven system
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

        self.test_only = testonly

    def settestmode(self, mode):
        """
        Set testMode to (mode).
        - True: run the Source in testMode;
        - False: run it in full mode
        :param mode:
        :return: None

        """

        self.test_mode = mode

    def getTestSuite(self):
        """
        An abstract method that should be overwritten with
        tests appropriate for the specific source.
        :return:

        """
        return None

    # TODO: pramaterising the release date

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

        with open(filename, 'r', encoding=encoding, newline=r'\n') as filereader:
            contents = filereader.read()
        contents = re.sub(r'\r', '', contents)
        with open(filename, "w") as filewriter:
            filewriter.truncate()
            filewriter.write(contents)

    @staticmethod
    def open_and_parse_yaml(yamlfile):
        """
        :param file: String, path to file containing label-id mappings in
                             the first two columns of each row
        :return: dict where keys are labels and values are ids
        """

        # ??? what if the yaml file does not contain a dict datastructure?
        mapping = dict()
        if os.path.exists(os.path.join(os.path.dirname(__file__), yamlfile)):
            map_file = open(os.path.join(os.path.dirname(__file__), yamlfile), 'r')
            mapping = yaml.safe_load(map_file)
            map_file.close()
        else:
            LOG.warning("file: %s not found", yamlfile)

        return mapping

    @staticmethod
    def parse_mapping_file(file):
        """
        :param file: String, path to file containing label-id mappings
                in the first two columns of each row
        :return: dict where keys are labels and values are ids
        """
        id_map = {}
        if os.path.exists(os.path.join(os.path.dirname(__file__), file)):
            with open(os.path.join(os.path.dirname(__file__), file)) as tsvfile:
                reader = csv.reader(tsvfile, delimiter="\t")
                for row in reader:
                    key = row[0]
                    value = row[1]
                    id_map[key] = value

        return id_map

    @staticmethod
    def _get_default_request_headers():
        return {
            'User-Agent': USER_AGENT
        }

    # @staticmethod
    # def getTestSuite(ingest):  # WIP
    #    '''
    #    try to avoid having one of these per ingest
    #    '''
    #    import unittest
    #    testcase = ingest + 'TestCase'
    #    # construct import names ... how
    #    from tests.test_ + ingest import testcase
    #    return unittest.TestLoader().loadTestsFromTestCase(testcase)

    def load_local_translationtable(self, name):
        '''
        Load "ingest specific" translation from whatever they called something
        to the ontology label we need to map it to.
        To facilitate seeing more ontology labels in dipper ingests
        a reverse mapping from ontology labels to external strings is also generated
        and available as a dict localtcid

        '---\n# %s.yaml\n"": ""  # example'
        '''

        localtt_file = '../../translationtable/' + name + '.yaml'

        try:
            with open(os.path.join(os.path.dirname(__file__), localtt_file)):
                pass
        except IOError:
            # write a stub file as a place holder if none exists
            with open(os.path.join(os.path.dirname(__file__),
                                   localtt_file), 'w') as write_yaml:
                print('---\n# %s.yaml\n"": ""  # example' % name, file=write_yaml)
        finally:
            with open(os.path.join(os.path.dirname(__file__),
                                   localtt_file), 'r') as read_yaml:
                localtt = yaml.safe_load(read_yaml)

        # inverse local translation.
        # note: keeping this invertable will be work.
        # Useful to not litter an ingest with external syntax
        self.localtcid = {v: k for k, v in localtt.items()}

        return localtt

    def resolve(self, word, mandatory=True, default=None):
        '''
        composite mapping
        given f(x) and g(x)
        here: localtt & globaltt respectivly
        return g(f(x))|g(x)||f(x)|x in order of preference
        returns x on fall through if finding a mapping
        is not mandatory (by default finding is mandatory).

        This may be specialized further from any mapping
        to a global mapping only; if need be.

        :param word:  the srting to find as a key in translation tables
        :param mandatory: boolean to cauae failure when no key exists
        :param default: string to return if nothing is found (& not manandatory)
        :return
            value from global translation table,
            or value from local translation table,
            or the query key if finding a value is not mandatory (in this order)

        '''

        assert word is not None

        # we may not agree with a remote sources use of a global term we have
        # this provides opportunity for us to override
        if word in self.localtt:
            label = self.localtt[word]
            if label in self.globaltt:
                term_id = self.globaltt[label]
            else:
                logging.info(
                    "Translated to '%s' but no global term_id for: '%s'", label, word)
                term_id = label
        elif word in self.globaltt:
            term_id = self.globaltt[word]
        else:
            if mandatory:
                raise KeyError("Mapping required for: ", word)
            logging.warning("We have no translation for: '%s'", word)

            if default is not None:
                term_id = default
            else:
                term_id = word
        return term_id

    @staticmethod
    def check_fileheader(expected, received):
        '''
        Compare file headers received versus file headers expected
        if the expected headers are a subset (proper or not)
        of received headers report suscess (warn if proper subset)

            param:  expected  list
            param:  received  list

            return: truthyness
        '''
        exp = set(expected)
        got = set(received)
        if expected != received:
            LOG.error('\nExpected header: %s\nRecieved header: %s', expected, received)

            # pass reordering and adding new columns (after protesting)
            # hard fail on missing expected columns (temper with mandatory cols?)
            if exp - got != set():
                LOG.error('Missing: %s', exp - got)
                raise AssertionError('Incomming headers are missing expected column.')

            if got - exp != set():
                LOG.warning('Addtional new columns: %s', got - exp)
            else:
                LOG.warning('Check columns order')

        return (exp ^ got) & exp == set()
