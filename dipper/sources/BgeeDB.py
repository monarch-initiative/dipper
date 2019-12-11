
import logging
import os
import time
import gzip

import sqlite3
from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.assoc.Association import Assoc
from dipper.graph.StreamedGraph import StreamedGraph


LOG = logging.getLogger(__name__)
BGEE_FTP = 'ftp.bgee.org'
CACHE = 'https://archive.monarchinitiative.org/cache'


class BgeeDB(Source):
    """
    Bgee is a database to retrieve and compare gene expression
    patterns between species.

    Bgee first maps heterogeneous expression data (currently RNA-Seq,
    Affymetrix, in situ hybridization, and EST data) to anatomy and
    development of different species.

    Then, in order to perform automated cross species comparisons,
    homology relationships across anatomies, and comparison criteria
    between developmental stages, are designed.
    """

    default_species = [
        "Cavia porcellus",                  # guinea pig
        "Mus musculus",                     # mouse
        "Rattus norvegicus",                # rat
        "Monodelphis domestica",            # gray short-tailed opossum
        "Anolis carolinensis",              # lizard
        "Caenorhabditis elegans",           # worm
        "Drosophila melanogaster",          # fly
        "Danio rerio",                      # zebrafish
        "Xenopus (Silurana) tropicalis",    # frog
        "Gallus gallus",                    # chicken
        "Ornithorhynchus anatinus",         # platypus
        "Erinaceus europaeus",              # hedgehog
        "Macaca mulatta",                   # monkey (rhesus macaque)
        "Gorilla gorilla",                  # gorilla
        "Pan paniscus",                     # bonobo
        "Pan troglodytes",                  # chimp
        "Homo sapiens",                     # corporal people
        "Canis lupus familiaris",           # dog  "Canis familiaris"
        "Felis catus",                      # cat
        "Equus caballus",                   # horse
        "Sus scrofa",                       # pig
        "Bos taurus",                       # cow
        "Oryctolagus cuniculus",            # rabbit
        # 7217	Drosophila_ananassae
        # 7230	Drosophila_mojavensis
        # 7237	Drosophila_pseudoobscura
        # 7240	Drosophila_simulans
        # 7244	Drosophila_virilis
        # 7245	Drosophila_yakuba
    ]
    files = {
        'bgeedb_lite': {
            # 'url': 'ftp://ftp.bgee.org/current/sql_lite_dump.tar.gz', # pre fetched
            # 'file': 'sql_lite_dump.tar.gz',   # un processed
            # zcat sql_lite_dump.tar.gz |
            #   scripts/mysql2sqlite -  |
            #   sqlite3 bgee.sqlite3


            'url':  CACHE + '/bgee.sqlite.gz',  # processed
            'file': 'bgee.sqlite.gz',

            'tables': {  # see resources/bgee for generator hints
                'anatEntity': {
                    'columns':  [
                        'anatEntityId',
                        'anatEntityName',
                        'anatEntityDescription',
                    ]
                },
                'gene': {
                    'columns': [
                        'bgeeGeneId',
                        'geneId',
                        'geneName',
                        'geneDescription',
                        'speciesId',
                    ]
                },
                'species': {
                    'columns': [
                        'speciesId',
                        'genus',
                        'species',
                        'speciesCommonName',
                        'genomeVersion',
                        'genomeSpeciesId',
                    ]
                },
                'globalCond': {
                    'columns': [
                        'globalConditionId',
                        'anatEntityId',
                        'stageId',
                        'speciesId',
                    ]
                },
                'stage': {
                    'columns': [
                        'stageId',
                        'stageName',
                        'stageDescription',
                    ]
                },
                'globalExpression': {
                    'columns': [
                        'globalExpressionId',
                        'bgeeGeneId',
                        'globalConditionId',
                        'summaryQuality',
                    ]
                },
            }
        }
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None,
                 tax_ids=None,
                 version=None):
        """
        :param tax_ids: [str,], List of NCBI taxon  identifiers
        :return:
        """
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='bgee',
            ingest_title='Bgee Gene expression data in animals',
            ingest_url='http://bgee.org/',
            ingest_logo='source-bgee.png',
            # license_url=None,
            data_rights='https://bgee.org/?page=about'
            # file_handle=None
        )
        self.default_taxa = {
            x: self.globaltt[x].split(':')[-1] for x in self.default_species}
        # names for logging
        self.txid_name = {v: k for k, v in self.default_taxa.items()}

        if tax_ids is None:
            tax_ids = self.default_taxa.values()
        self.tax_ids = [str(x) for x in tax_ids]  # incase they were passed in
        LOG.info(
            "Filtering on tax_ids %s",
            [{t: self.txid_name[t]} for t in self.tax_ids])

        if version is None:
            self.version = 'current'
        else:
            self.version = version

        # Database connection and cursor
        self.con = None  # sqlite3.connect('raw/bgee/bgee.sqlite3?mode=ro')
        self.cur = None  # self.con.cursor()

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced: boolean, force download

        """
        # get the remote copy unless the local copy is even newer
        self.get_files('files'=None, 'is_dl_forced'=is_dl_forced)

        # if the unzipped version is olere than the zipped version then unzip it
        bgeedb_gzip =  '/'.join((self.rawdir, self.files['bgeedb_lite']['file'])



        self.con = sqlite3.connect(bgeedb + '?mode=ro')
        self.con.isolation_level = None  # None for autocommit
        self.cur = self.con.cursor()

        # zcat omia.sql.gz | ./mysql2sqlite - > omia_sqlite.sql
        with open('bgee.sqlite3', 'r') as sql:
            self.cur.executescript(sql.read())

        # update internal metadata (incase it does not come with)
        self.cur.execute('vacuum;')  # a couple of minutes
        self.cur.execute('analyze;') # a few seconds



    #######################################################
        for dlname in files_to_download:
            localfile = '/'.join((self.rawdir, dlname))
            info = ftp.sendcmd("MLST {}".format(dlname))   # fetch remote file stats
            LOG.info(
                '%s\n'
                'Remote File Size: %i\n'
                'Remote timestamp: %s',
                dlname, int(info['size']),
                self._convert_ftp_time_to_iso(info['modify']))

            if not os.path.exists(localfile) or is_dl_forced or \
                    self.check_if_remote_is_newer(
                            localfile, int(info['size']), info['modify']):

                LOG.info("Fetching %s", dlname)
                LOG.info("Writing to %s", localfile)



    def parse(self, limit=None):
        """
        Given the input taxa, expects files in the raw directory
        with the name {tax_id}_anat_entity_all_data_Pan_troglodytes.tsv.zip

        :param limit: int Limit to top ranked anatomy associations per group
        :return: None
        """

        files_to_download, ftp = self._get_file_list(
            self.files['anat_entity']['path'],
            self.files['anat_entity']['pattern'], None)
        for dlname in files_to_download:
            localfile = '/'.join((self.rawdir, dlname))
            with gzip.open(localfile, 'rt', encoding='ISO-8859-1') as fh:
                LOG.info("Processing %s", localfile)
                self._parse_gene_anatomy(fh, limit)
