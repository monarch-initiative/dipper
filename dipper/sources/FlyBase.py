import logging
import re
import csv
import gzip
import os
from ftplib import FTP
from typing import Dict

from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.assoc.G2PAssoc import G2PAssoc


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
    4. fbrf_pmid_pmcid_doi_fb_*.tsv.gz -  flybase ref to pmid

    We connect using the
    [Direct Chado Access](http://gmod.org/wiki/Public_Chado_Databases#Direct_Chado_Access)

    When running the whole set,
    it performs best by dumping raw triples using the flag ```--format nt```.

    Note that this script underwent a major revision after commit bd5f555
    in which genotypes, stocks, and environments were removed

    """
    FLYFTP = 'ftp.flybase.net'
    FLYFILES = '/releases/current/precomputed_files'

    resources = [
        {
            'query': '../../resources/sql/fb/allele_phenotype.sql',
            'outfile': 'allele_phenotype.tsv',
            'columns': [
                'allele_id',
                'pheno_desc',
                'pub_id',
                'pub_title',
                'pmid_id',
            ]
        },
        {
            'query': '../../resources/sql/fb/gene_xref.sql',
            'outfile': 'gene_xref.tsv',
            'columns': [
                'gene_id',
                'xref_id',
                'xref_source',
            ]
        }
    ]

    files = {
        'disease_model': {
            'file': 'allele_human_disease_model_data.tsv.gz',
            'url':  r'human_disease/allele_human_disease_model_data_fb_.*\.tsv\.gz',
            'columns': [
                'FBal_ID',
                'AlleleSymbol',
                'DOID_qualifier',
                'DOID_term',
                'DOID_ID',
                'Evidence/interacting_alleles',
                'Reference_FBid',
            ]
        },
        'species_map': {
            'file': 'species.ab.gz',
            'url': r'species/species\.ab\.gz',
            'columns': [
                'internal_id',
                'taxgroup',
                'abbreviation',
                'genus',
                'species name',
                'common name',
                'comment',
                'ncbi-taxon-id',
            ]
        },
        'allele_gene': {
            'file': 'fbal_to_fbgn_fb.tsv.gz',
            'url': r'alleles/fbal_to_fbgn_fb_(.*)\.tsv\.gz',
            'columns': [
                'AlleleID',
                'AlleleSymbol',
                'GeneID',
                'GeneSymbol',
            ]
        },
        'ref_pubmed': {
            'file': 'fbrf_pmid_pmcid_doi_fb.tsv.gz',
            'url': r'references/fbrf_pmid_pmcid_doi_fb_.*\.tsv\.gz',
            'columns': [
                'FBrf',
                'PMID',
                'PMCID',
                'DOI',
                'pub_type',
                'miniref',
                'pmid_added',
            ]
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

    def fetch(self, is_dl_forced=False):
        """
        Fetch flat files and sql queries

        :param is_dl_forced: force download
        :return: None

        """

        # create the connection details for Flybase
        cxn = {
            'host': 'chado.flybase.org', 'database': 'flybase', 'port': 5432,
            'user': 'flybase', 'password': 'no password'}

        self.dataset.setFileAccessUrl(
            ''.join(('jdbc:postgresql://', cxn['host'], ':', str(cxn['port']),
                     '/', cxn['database'])), is_object_literal=True)

        # Get flat files
        ftp = FTP(FlyBase.FLYFTP)
        ftp.login("anonymous", "info@monarchinitiative.org")

        for src_key, file in self.files.items():
            filename = self._resolve_filename(file['url'], ftp)
            # prepend ftp:// since this gets added to dataset rdf model
            self.files[src_key]['url'] = "ftp://" + self.FLYFTP + filename

        ftp.close()

        self.get_files(is_dl_forced)

        # Get data from remote db
        # Each query takes 2 minutes or so
        for query_map in self.resources:
            query_fh = open(os.path.join(
                os.path.dirname(__file__), query_map['query']), 'r')
            query = query_fh.read()
            self.fetch_query_from_pgdb(
                query_map['outfile'], query, None, cxn)

    def parse(self, limit=None):
        """
        Parse flybase files and add to graph

        :param limit: number of rows to process
        :return: None

        """
        if limit is not None:
            LOG.info("Only parsing first %d rows of each file", limit)
        LOG.info("Parsing files...")

        self._process_disease_model(limit)

        LOG.info("Finished parsing.")
        LOG.info("Loaded %d nodes", len(self.graph))

    def _flyref_to_pmid(self) -> Dict[str, str]:
        """
        Generates a dictionary of flybase reference and PMID curie mappings;
        "FBrf0241315":"PMID:30328653"

        :raises TypeError: If len(row) is different from len(columns)
        :raises ValueError: If headers differ from FlyBase.files[src_key]['columns']
        :return: Dict with FBrf ids as keys and PMID curies as values
        """
        pub_map = {}
        src_key = 'ref_pubmed'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("creating map of flybase ref ids and pmids")

        col = self.files[src_key]['columns']

        with gzip.open(raw, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            # skip first four lines
            for _ in range(0, 2):
                # file name, note
                next(reader)

            row = next(reader)  # headers
            row[0] = row[0][1:]  # uncomment

            if self.check_fileheader(col, row):
                pass

            next(reader)  # Go to next line

            for row in reader:
                # File ends with a final comment
                if ''.join(row).startswith('#'):
                    continue

                pmid = row[col.index('PMID')]
                fly_ref = row[col.index('FBrf')]

                pub_map[fly_ref] = 'PMID:' + pmid

            return pub_map

    def _process_disease_model(self, limit):
        """
        Make associations between a disease and fly alleles
        Adds triples to self.graph

        Pulls FBrf to pmid eqs from _flyref_to_pmid
        for 2019_02 maps all but FlyBase:FBrf0211649

        Right now only model_of is processed from DOID_qualifier
        As of 2019_02 release there are also:
        ameliorates
        exacerbates
        DOES NOT ameliorate
        DOES NOT exacerbate
        DOES NOT model
        DOID_qualifier

        :raises TypeError: If len(row) is different from len(columns)
        :raises ValueError: If headers differ from FlyBase.files[src_key]['columns']
        :param limit: number of rows to process
        :return: None

        """
        graph = self.graph
        pub_map = self._flyref_to_pmid()
        src_key = 'disease_model'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("processing disease models")

        col = self.files[src_key]['columns']

        with gzip.open(raw, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            # skip first four lines
            for _ in range(0, 4):
                # file name, generated datetime, db info, blank line
                next(reader)

            row = next(reader)  # headers
            row[0] = row[0][2:]  # uncomment

            if self.check_fileheader(col, row):
                pass

            next(reader) # Go to next line

            for row in reader:
                # File ends with a blank line and a final comment
                if ''.join(row).startswith('#') or not ''.join(row).strip():
                    continue

                allele_id = row[col.index('FBal_ID')]
                flybase_ref = row[col.index('Reference_FBid')]
                evidence_or_allele = row[col.index('Evidence/interacting_alleles')]
                doid_id = row[col.index('DOID_ID')]
                doid_qualifier = row[col.index('DOID_qualifier')]

                allele_curie = 'FlyBase:' + allele_id
                if doid_qualifier == 'model of':
                    relation = self.globaltt['is model of']
                else:
                    # amelorates, exacerbates, and DOES NOT *
                    continue

                assoc = G2PAssoc(graph, self.name, allele_curie, doid_id, relation)
                if flybase_ref != '':
                    ref_curie = None
                    try:
                        ref_curie = pub_map[flybase_ref]
                    except KeyError:
                        ref_curie = 'FlyBase:' + flybase_ref

                    assoc.add_source(ref_curie)
                if evidence_or_allele == 'inferred from mutant phenotype':
                    evidence_id = self.globaltt['mutant phenotype evidence']
                    assoc.add_evidence(evidence_id)
                else:
                    assoc.set_description(evidence_or_allele)

                assoc.add_association_to_graph()

                if limit is not None and reader.line_num > limit:
                    break

    @staticmethod
    def _resolve_filename(filename: str, ftp: FTP) -> str:
        """
        Resolve a file name from ftp server given a regex
        :return: str, filepath on ftp server
        """

        # Represent file path as a list of directories
        dir_path = FlyBase.FLYFILES.split('/')
        # Also split on the filename to get prepending dirs
        file_path = filename.split('/')
        file_regex = file_path.pop()
        working_dir = "/".join(dir_path + file_path)

        LOG.info('Looking for remote files in %s', working_dir)

        ftp.cwd(working_dir)
        remote_files = ftp.nlst()

        files_to_download = [dnload for dnload in remote_files
                             if re.match(file_regex, dnload)]

        if len(files_to_download) > 1:
            raise ValueError("Could not resolve filename from regex, "
                             "too many matches for {}, matched: {}"
                             .format(file_regex, files_to_download))
        if not files_to_download:
            raise ValueError("Could not resolve filename from regex, "
                             "no matches for {}".format(file_regex))

        return working_dir + '/' + files_to_download[0]
