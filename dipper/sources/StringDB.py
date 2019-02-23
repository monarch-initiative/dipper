import logging
import gzip

import pandas as pd
from dipper.sources.Source import Source
from dipper.sources.Ensembl import Ensembl

LOG = logging.getLogger(__name__)

STRING_BASE = "http://string-db.org/download/"
DEFAULT_TAXA = ['9606', '10090', '7955', '7227', '6239']


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

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None, version=None):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'string',
            ingest_title='Known and predicted protein-protein interactions',
            ingest_url='https://string-db.org/',
            license_url=None,
            data_rights='https://string-db.org/cgi/access.pl?footer_active_subpage=licensing'
            # file_handle=None
        )

        if tax_ids is None:
            self.tax_ids = DEFAULT_TAXA
        else:
            LOG.info("Filtering on taxa %s", tax_ids)
            self.tax_ids = [str(tax_id) for tax_id in tax_ids]

        if version is None:
            self.version = 'v10.5'
        elif not version.startswith('v'):
            self.version = 'v' + version
        else:
            self.version = version

        self.files = {
            'protein_links': {
                'path': '{}protein.links.detailed.{}/'.format(
                    STRING_BASE, self.version),
                'pattern': 'protein.links.detailed.{}.txt.gz'.format(self.version)
            }
        }

        self.id_map_files = {
            '9606': {
                'url': 'https://string-db.org/mapping_files/entrez/'
                       'human.entrez_2_string.2018.tsv.gz',
                'file': 'human.entrez_2_string.2018.tsv.gz'
            },
            '10090': {
                'url': 'https://string-db.org/mapping_files/entrez/'
                       'mouse.entrez_2_string.2018.tsv.gz',
                'file': 'mouse.entrez_2_string.2018.tsv.gz'
            },
            '6239': {
                'url': 'https://string-db.org/mapping_files/entrez/'
                       'celegans.entrez_2_string.2018.tsv.gz',
                'file': 'celegans.entrez_2_string.2018.tsv.gz'
            },
            '7227': {
                'url': 'https://string-db.org/mapping_files/entrez/'
                       'fly.entrez_2_string.2018.tsv.gz',
                'file': 'fly.entrez_2_string.2018.tsv.gz'
            },
            '7955': {
                'url': 'https://string-db.org/mapping_files/entrez/'
                       'zebrafish.entrez_2_string.2018.tsv.gz',
                'file': 'zebrafish.entrez_2_string.2018.tsv.gz'
            },
            '4932': {
                'url': 'https://string-db.org/mapping_files/entrez/'
                       'yeast.entrez_2_string.2018.tsv.gz',
                'file': 'yeast.entrez_2_string.2018.tsv.gz'
            }
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
        self.get_files(is_dl_forced, file_paths)
        self.get_files(is_dl_forced, self.id_map_files)

        return

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

        for taxon in protein_paths:
            ensembl = Ensembl(self.graph_type, self.are_bnodes_skized)
            string_file_path = '/'.join((
                self.rawdir, protein_paths[taxon]['file']))

            fh = gzip.open(string_file_path, 'rb')
            dataframe = pd.read_csv(fh, sep='\s+')
            fh.close()
            p2gene_map = dict()

            if taxon in self.id_map_files:
                LOG.info("Using string provided id_map files")
                map_file = '/'.join((self.rawdir, self.id_map_files[taxon]['file']))

                with gzip.open(map_file, 'rt') as reader:
                    for line in reader.readlines():
                        if line.startswith('#'): continue
                        tax, gene, prot = line.rstrip("\n").split("\t")
                        genes = gene.split('|')
                        p2gene_map[prot.replace(taxon + '.', '')] = ["NCBIGene:"+ entrez_id
                                                                     for entrez_id in genes]
            else:
                LOG.info("Fetching ensembl proteins for taxon %s", taxon)
                p2gene_map = ensembl.fetch_protein_gene_map(taxon)
                for key in p2gene_map.keys():
                    for i, gene in enumerate(p2gene_map[key]):
                        p2gene_map[key][i] = "ENSEMBL:{}".format(gene)

            LOG.info(
                "Finished fetching ENSP ID mappings, fetched %i proteins",
                len(p2gene_map))

            LOG.info(
                "Fetching protein protein interactions for taxon %s", taxon)

            self._process_protein_links(dataframe, p2gene_map, taxon, limit)

    def _process_protein_links(self, dataframe, p2gene_map, taxon,
                               limit=None, rank_min=700):
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
                if protein1 > protein1:
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
            "Finished parsing p-p interactions for %s, " +
            "%i rows filtered out based on checking ensembl proteins",
            taxon, filtered_out_count)
        return

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
                    self.files[file_type]['pattern'])
            }
        return file_paths
