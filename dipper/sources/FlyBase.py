import logging
import re
import csv
import gzip
import io
import os

from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.Model import Model
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.Environment import Environment
from dipper.utils.DipperUtil import DipperUtil


LOG = logging.getLogger(__name__)


class FlyBase(PostgreSQLSource):
    """
    This is the [Drosophila Genetics](http://www.flybase.org/) resource,
    from which we process genotype and phenotype data about the fruit fly.

    Here, we connect to their public database and download preprocessed files

    Queries from the relational db
    1. allele-phenotype data: ../../sources/sql/fb/allele_phenotype.sql
    2. gene dbxrefs: ../../resources/sql/fb/gene_xref.sql

    Downloads:
    1. allele_human_disease_model_data_fb_*.tsv.gz - models of disease
    2. species.ab.gz - species prefix mappings
    3. fbal_to_fbgn_fb_*.tsv.gz - allele to gene

    We connect using the
    [Direct Chado Access](http://gmod.org/wiki/Public_Chado_Databases#Direct_Chado_Access)

    When running the whole set,
    it performs best by dumping raw triples using the flag ```--format nt```.

    Note that this script underwent a major revision after commit bd5f555

    """

    FLYFTP = 'ftp://ftp.flybase.net/releases/current/precomputed_files/'

    resources = [
        {
            'query': '../../resources/sql/fb/allele_phenotype.sql',
            'outfile': 'allele_phenotype'
        },
        {
            'query': '../../resources/sql/fb/gene_xref.sql',
            'outfile': 'gene_xref'
        }
    ]

    files = {
        'disease_model': {
            'file': 'allele_human_disease_model_data.tsv.gz',
            'url':  FLYFTP +
                    'human_disease/allele_human_disease_model_data_fb_*.tsv.gz'
        },
        'species_map': {
            'file': 'species.ab.gz',
            'url': FLYFTP + 'species/species.ab.gz'
        },
        'allele_gene': {
            'file': 'fbal_to_fbgn_fb.tsv.gz',
            'url': FLYFTP + 'alleles/fbal_to_fbgn_fb_*.tsv.gz'
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'flybase',
            ingest_title='FlyBase',
            ingest_url='http://www.flybase.org/',
            license_url=None,
            data_rights='https://wiki.flybase.org/wiki/FlyBase_Wiki:General_disclaimer',
            file_handle=None
            )

        return

    def fetch(self, is_dl_forced=False):
        """
        :return: None
        """

        # create the connection details for Flybase
        cxn = {
            'host': 'chado.flybase.org', 'database': 'flybase', 'port': 5432,
            'user': 'flybase', 'password': 'no password'}

        self.dataset.setFileAccessUrl(
            ''.join(('jdbc:postgresql://', cxn['host'], ':', str(cxn['port']),
                     '/', cxn['database'])), is_object_literal=True)

        for query_map in self.resources:
            query_fh = open(os.path.join(
                os.path.dirname(__file__), query_map['query']), 'r')
            query = query_fh.read()
            self.fetch_query_from_pgdb(
                query_map['outfile'], query, None, cxn)

        self.get_files(False)

    def parse(self, limit=None):
        """
        :param limit: Only parse this many lines of each table
        :return: None

        """
        if limit is not None:
            LOG.info("Only parsing first %d rows of each file", limit)
        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self._process_disease_model(limit)

        LOG.info("Finished parsing.")
        LOG.info("Loaded %d nodes", len(self.graph))


    def _process_disease_model(self, limit):
        """
        Here we make associations between a disease and the supplied "model".
        In this case it's an allele.
        :param limit:
        :return: None

        """

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        raw = '/'.join((self.rawdir, self.files['disease_models']['file']))
        LOG.info("processing disease models")

        line_counter = 0

        with gzip.open(raw, 'rb') as f:
            filereader = csv.reader(
                io.TextIOWrapper(f, newline=""),
                delimiter='\t', quotechar='\"')
            for line in filereader:
                # skip comments
                if re.match(r'#', ''.join(line)) or ''.join(line) == '':
                    continue
                (allele_id, allele_symbol, qualifier, doid_label, doid_id,
                 evidence_or_interacting_allele, pub_id) = line
                line_counter += 1

                allele_id = 'FlyBase:' + allele_id
                if qualifier == 'model of':
                    relation = self.globaltt['is model of']
                else:
                    # TODO amelorates, exacerbates, and DOES NOT *
                    continue

                assoc = G2PAssoc(graph, self.name, allele_id, doid_id, relation)
                if pub_id != '':
                    pub_id = 'FlyBase:' + pub_id
                    assoc.add_source(pub_id)
                if evidence_or_interacting_allele == 'inferred from mutant phenotype':
                    evidence_id = self.globaltt['mutant phenotype evidence']
                    assoc.add_evidence(evidence_id)
                else:
                    assoc.set_description(evidence_or_interacting_allele)

                assoc.add_association_to_graph()

                if not self.test_mode and limit is not None and line_counter > limit:
                    break
