import logging
import re

from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper import config
from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.GenomicFeature import Feature, makeChromID


logger = logging.getLogger(__name__)


class FlyBase(Source):
    """
    This is the [Drosophila Genetics](http://www.flybase.org/) resource,
    from which we process genotype and phenotype data about fruitfly.
    Genotypes leverage the GENO genotype model.

    Here, we connect to their public database, and download a subset of tables/views to get specifically at the
    geno-pheno data, then iterate over the tables.  We end up effectively performing joins when adding nodes
    to the graph.
    We connect using the [Direct Chado Access](http://gmod.org/wiki/Public_Chado_Databases#Direct_Chado_Access)

    """


    tables = [
        'genotype',
        'feature_genotype',
        'pub', 'feature_pub',
        'pub_dbxref', 'feature_dbxref',
        'feature_relationship',
        'cvterm',
        'stock_genotype',
        'stock',
        'organism',
        'environment',
        'phenotype',
        'phenstatement',
        'dbxref',
        'db',
        'phenotype_cvterm',
        # 'phendesc',
        # 'strain',
        # 'environment_cvterm',

    ]

    def __init__(self):
        Source.__init__(self, 'flybase')
        self.namespaces.update(curie_map.get())

        # update the dataset object with details about this resource
        self.dataset = Dataset('flybase', 'FlyBase', 'http://www.flybase.org/', None,
                               None, 'http://flybase.org/wiki/FlyBase:FilesOverview')

        # check if config exists; if it doesn't, error out and let user know
        if 'dbauth' not in config.get_config() and 'mgi' not in config.get_config()['dbauth']:
            logger.error("not configured with PG user/password.")

        # source-specific warnings.  will be cleared when resolved.
        logger.warn("we are ignoring normal phenotypes for now")

        # so that we don't have to deal with BNodes, we will create hash lookups for the internal identifiers
        # the hash will hold the type-specific-object-keys to FB public identifiers.  then, subsequent
        # views of the table will lookup the identifiers in the hash.  this allows us to do the 'joining' on the
        # fly
        self.idhash = {'allele': {}, 'marker': {}, 'publication': {}, 'strain': {},
                       'genotype': {}, 'annot': {}, 'notes': {}}
        self.markers = {'classes': [], 'indiv': []}  # to store if a marker is a class or indiv

        self.label_hash = {}  # use this to store internally generated labels for various features
        self.geno_bkgd = {}  # use this to store the genotype strain ids for building genotype labels

        self.wildtype_alleles = set()

        return

    def fetch(self, is_dl_forced=False):
        """
        For the MGI resource, we connect to the remote database, and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        """

        # create the connection details for MGI
        cxn = {'host': 'flybase.org', 'database': 'flybase', 'port': 5432, 'user': 'flybase',
               'password': 'no password'}

        self.dataset.setFileAccessUrl(''.join(('jdbc:postgresql://', cxn['host'], ':',
                                               str(cxn['port']), '/', cxn['database'])))

        # process the tables
        # self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables, cxn, None, is_dl_forced)

        return

    def parse(self, limit=None):
        """
        We process each of the postgres tables in turn.  The order of processing is important here, as we build
        up a hashmap of internal vs external identifers (unique keys by type to FB id).  These include
        allele, marker (gene), publication, strain, genotype, annotation (association), and descriptive notes.
        :param limit: Only parse this many lines of each table
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %d rows of each file", limit)
        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        # the following will provide us the hash-lookups
        # These must be processed in a specific order


        logger.info("Finished parsing.")

        self.load_bindings()
        for g in [self.graph, self.testgraph]:
            Assoc(self.name).load_all_properties(g)
            gu = GraphUtils(curie_map.get())
            gu.loadAllProperties(g)

        logger.info("Loaded %d nodes", len(self.graph))
        return

    @staticmethod
    def _makeInternalIdentifier(prefix, key):
        """
        This is a special Flybase-to-MONARCH-ism.  Flybase tables have unique keys that we use here, but don't want
        to necessarily re-distribute those internal identifiers.  Therefore, we make them into keys in a consistent
        way here.
        :param prefix: the object type to prefix the key with, since the numbers themselves are not unique across tables
        :param key: the number (unique key)
        :return:
        """

        return '_fb'+prefix+'key'+key


    # def getTestSuite(self):
    #     import unittest
    #     from tests.test_mgi import MGITestCase
    #     # TODO test genotypes
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(MGITestCase)
    #
    #     return test_suite