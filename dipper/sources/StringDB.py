import logging
import gzip

import pandas as pd
from dipper.sources.Source import Source, USER_AGENT
from dipper.sources.Ensembl import Ensembl

LOG = logging.getLogger(__name__)

STRING = 'https://string-db.org'
STRING_DWN = STRING + "/download"
STRING_MAP = STRING + '/mapping_files/entrez'
DEFAULT_TAXA = ['9606', '10090', '7955', '7227', '6239', '4932', '10116']
# https://string-db.org/cgi/access.pl?footer_active_subpage=archive
# TODO automate it with:
# curl https://string-db.org/api/tsv-no-header/version
#   11.0        https://version-11-0.string-db.org
VERSION = '11.0'
YEAR = '2018'


class StringDB(Source):
    """
    STRING is a database of known and predicted protein-protein interactions.
    The interactions include direct (physical) and indirect
    (functional) associations; they stem from computational prediction,
    from knowledge transfer between organisms, and from interactions
    aggregated from other (primary) databases.
    From: http://string-db.org/cgi/about.pl?footer_active_subpage=content

    STRING uses one protein per gene. If there is more than one isoform
    per gene, we usually select the longest isoform, unless we have
    information that suggest that other isoform regarded as
    cannonical (e.g., proteins in the CCDS database).
    From: http://string-db.org/cgi/help.pl
    """

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None,
                 tax_ids=None,
                 version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='string',
            ingest_title='Known and predicted protein-protein interactions',
            ingest_url=STRING,
            ingest_logo='source-string.png',
            license_url=None,
            data_rights=STRING+'/cgi/access.pl?footer_active_subpage=licensing'
            # file_handle=None
        )

        if tax_ids is None:
            self.tax_ids = DEFAULT_TAXA
        else:
            LOG.info("Filtering on taxa %s", tax_ids)
            self.tax_ids = [str(tax_id) for tax_id in tax_ids]

        if version is None:
            version = VERSION
        self.version = version

        self.files = {
            'protein_links': {
                # https://stringdb-static.org/download/protein.links.detailed.v11.0/
                'path': '{}/protein.links.detailed.v{}/'.format(
                    STRING_DWN, self.version),
                # 9606.protein.links.detailed.v11.0.txt.gz
                'pattern': 'protein.links.detailed.v{}.txt.gz'.format(self.version)
            }
        }

        self.id_map_files = {
            '9606': {
                'url': STRING_MAP + '/human.entrez_2_string.2018.tsv.gz',
                'file': 'human.entrez_2_string.'+YEAR+'.tsv.gz',
                'headers': {'User-Agent': USER_AGENT}
            },
            '10090': {
                'url': STRING_MAP + '/mouse.entrez_2_string.2018.tsv.gz',
                'file': 'mouse.entrez_2_string.'+YEAR+'.tsv.gz',
                'headers': {'User-Agent': USER_AGENT}
            },
            '6239': {
                'url': STRING_MAP + '/celegans.entrez_2_string.2018.tsv.gz',
                'file': 'celegans.entrez_2_string.'+YEAR+'.tsv.gz',
                'headers': {'User-Agent': USER_AGENT}
            },
            '7227': {
                'url': STRING_MAP + '/fly.entrez_2_string.2018.tsv.gz',
                'file': 'fly.entrez_2_string.'+YEAR+'.tsv.gz',
                'headers': {'User-Agent': USER_AGENT}
            },
            '7955': {
                'url': STRING_MAP + '/zebrafish.entrez_2_string.2018.tsv.gz',
                'file': 'zebrafish.entrez_2_string.'+YEAR+'.tsv.gz',
                'headers': {'User-Agent': USER_AGENT}
            },
            '4932': {
                'url': STRING_MAP + '/yeast.entrez_2_string.2018.tsv.gz',
                'file': 'yeast.entrez_2_string.'+YEAR+'.tsv.gz',
                'headers': {'User-Agent': USER_AGENT}
            },
            #'10116': {  # rat is not special enough to get its own mapping file
            #    'url': STRING_MAP + '/rat.entrez_2_string.2018.tsv.gz',
            #    'file': 'rat.entrez_2_string.'+YEAR+'.tsv.gz',
            #    'headers': {'User-Agent': USER_AGENT}
            #},
        }

    def fetch(self, is_dl_forced=False):
        """
        Override Source.fetch()
        Fetches resources from String

        We also fetch ensembl to determine if protein pairs are from
        the same species
        Args:
            :param is_dl_forced (bool): Force download
        Returns:
            :return None
        """

        file_paths = self._get_file_paths(self.tax_ids, 'protein_links')
        self.get_files(is_dl_forced, file_paths, delay=5)
        self.get_files(is_dl_forced, self.id_map_files, delay=5)

    def parse(self, limit=None):
        """
        Override Source.parse()
        Args:
            :param limit (int, optional) limit the number of rows processed
        Returns:
            :return None
        """
        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        protein_paths = self._get_file_paths(self.tax_ids, 'protein_links')
        col = ['NCBI taxid', 'entrez', 'STRING']
        for taxon in protein_paths:
            ensembl = Ensembl(self.graph_type, self.are_bnodes_skized)
            string_file_path = '/'.join((
                self.rawdir, protein_paths[taxon]['file']))
            p2gene_map = dict()
            with gzip.open(string_file_path, 'rb') as reader:
                dataframe = pd.read_csv(reader, sep=r'\s+')

            if taxon in self.id_map_files:
                LOG.info("Using string provided id_map files")
                map_file = '/'.join((self.rawdir, self.id_map_files[taxon]['file']))

                with gzip.open(map_file, 'rt') as reader:
                    line = next(reader)
                    row = line[2:-2].split(' / ')
                    if not self.check_fileheader(col, row):
                        pass
                    for line in reader.readlines():
                        row = line.rstrip('\n').split('\t')
                        # tax = row[col.index(''NCBI taxid')].strip()
                        gene = row[col.index('entrez')].strip()
                        prot = row[col.index('STRING')].strip()

                        genes = gene.split('|')
                        p2gene_map[prot.replace(taxon + '.', '')] = [
                            "NCBIGene:" + entrez_id for entrez_id in genes]
            else:
                LOG.info("Fetching ensembl protein_gene dict for NCBITaxon:%s", taxon)
                p2gene_map = ensembl.fetch_protein_gene_map(taxon)
                p2gene_map.update({k: ['ENSEMBL:' + p2gene_map[k]] for k in p2gene_map})

            LOG.info(
                "Finished fetching ENSP ID mappings, fetched %i proteins",
                len(p2gene_map))

            LOG.info(
                "Fetching protein protein interactions for taxon %s", taxon)

            self._process_protein_links(dataframe, p2gene_map, taxon, limit)

    def _process_protein_links(
            self, dataframe, p2gene_map, taxon, limit=None, rank_min=700):
        filtered_df = dataframe[dataframe['combined_score'] > rank_min]
        filtered_out_count = 0
        for index, row in filtered_df.iterrows():
            # Check if proteins are in same species
            protein1 = row['protein1'].replace('{}.'.format(taxon), '')
            protein2 = row['protein2'].replace('{}.'.format(taxon), '')
            gene1_curies = None
            gene2_curies = None
            try:
                # Keep orientation the same since RO!"interacts with" is symmetric
                # TEC: symeteric expansion is the job of post processing not ingest
                if protein1 >= protein2:
                    gene1_curies = p2gene_map[protein1]
                    gene2_curies = p2gene_map[protein2]
                else:
                    gene1_curies = p2gene_map[protein2]
                    gene2_curies = p2gene_map[protein1]
            except KeyError:
                filtered_out_count += 1

            if gene1_curies is not None and gene2_curies is not None:
                for gene1 in gene1_curies:
                    for gene2 in gene2_curies:
                        self.graph.addTriple(
                            gene1, self.globaltt['interacts with'], gene2)
                if limit is not None and index >= limit:
                    break

        LOG.info(
            "Finished parsing p-p interactions for %s, "
            "%i rows filtered out based on checking ensembl proteins",
            taxon, filtered_out_count)

    def _get_file_paths(self, tax_ids, file_type):
        """
        Assemble file paths from tax ids
        Args:
            :param tax_ids (list) list of taxa
        Returns:
            :return file dict
        """
        file_paths = dict()
        if file_type not in self.files:
            raise KeyError("file type {} not configured".format(file_type))
        for taxon in tax_ids:
            file_paths[taxon] = {
                'file': "{}.{}".format(taxon, self.files[file_type]['pattern']),
                'url': "{}{}.{}".format(
                    self.files[file_type]['path'], taxon,
                    self.files[file_type]['pattern']),
                'headers': {'User-Agent': USER_AGENT}
            }
        return file_paths
