from dipper.sources.Source import Source
from dipper.sources.Ensembl import Ensembl
from dipper.models.Dataset import Dataset
import logging
import gzip
import pandas as pd

logger = logging.getLogger(__name__)


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
    STRING_BASE = "http://string-db.org/download/"
    DEFAULT_TAXA = [9606, 10090, 7955, 7227, 6239]

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None, version=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'string')
        self.dataset = Dataset(
            'string', 'String', 'http://string-db.org/', None,
            'http://string-db.org/cgi/access.pl?footer_active_subpage=licensing')

        if tax_ids is None:
            self.tax_ids = StringDB.DEFAULT_TAXA
        else:
            logger.info("Filtering on taxa {}".format(tax_ids))
            self.tax_ids = tax_ids

        if version is None:
            self.version = 'v10.5'

        self.files = {
            'protein_links': {
                'path': '{}protein.links.detailed.{}/'.format(
                    StringDB.STRING_BASE, self.version),
                'pattern': 'protein.links.detailed.{}.txt.gz'.format(self.version)
            }
        }

        self.id_map_files = {
            9606: {
                'url': 'https://string-db.org/mapping_files/entrez_mappings/'
                       'entrez_gene_id.vs.string.v10.28042015.tsv',
                'file': 'entrez_gene_id.vs.string.v10.28042015.tsv'
            },
            10090: {
                'url': 'https://data.monarchinitiative.org/dipper/'
                       'cache/10090.string2mgi.tsv',
                'file': '10090.string2mgi.tsv'
            },
            6239: {
                'url': 'https://data.monarchinitiative.org/dipper/'
                       'cache/6239.string2ensembl_gene.tsv',
                'file': '6239.string2ensembl_gene.tsv'
            },
            7227: {
                'url': 'https://data.monarchinitiative.org/dipper/'
                       'cache/7227.string2ensembl_gene.tsv',
                'file': '7227.string2ensembl_gene.tsv'
            },
            7955: {
                'url': 'https://data.monarchinitiative.org/dipper/'
                       'cache/7955.string2zfin.tsv',
                'file': '7955.string2zfin.tsv'
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
            logger.info("Only parsing first %d rows", limit)

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
                map_file = '/'.join((
                    self.rawdir, self.id_map_files[taxon]['file']))

                mfile_handle = open(map_file, 'r')
                if taxon == 9606:
                    for line in mfile_handle.readlines():
                        gene, prot = line.rstrip("\n").split("\t")
                        p2gene_map[prot.replace('9606.', '')] \
                            = "NCBIGene:{}".format(gene)
                else:
                    for line in mfile_handle.readlines():
                        prot, gene = line.rstrip("\n").split("\t")
                        p2gene_map[prot] = gene
                mfile_handle.close()
            else:
                logger.info("Fetching ensembl proteins "
                            "for taxon {}".format(taxon))
                p2gene_map = ensembl.fetch_protein_gene_map(taxon)
                for key in p2gene_map.keys():
                    p2gene_map[key] = "ENSEMBL:{}".format(p2gene_map[key])
            if taxon == 9606:
                temp_map = ensembl.fetch_protein_gene_map(taxon)
                for key in temp_map:
                    if key not in p2gene_map:
                        p2gene_map[key] = "ENSEMBL:{}".format(temp_map[key])

            logger.info("Finished fetching ENSP ID mappings, "
                        "fetched {} proteins".format(len(p2gene_map)))

            logger.info("Fetching protein protein interactions "
                        "for taxon {}".format(taxon))

            self._process_protein_links(dataframe, p2gene_map, taxon, limit)

    def _process_protein_links(self, dataframe, p2gene_map, taxon,
                               limit=None, rank_min=700):
        filtered_df = dataframe[dataframe['combined_score'] > rank_min]
        filtered_out_count = 0
        for index, row in filtered_df.iterrows():
            # Check if proteins are in same species
            protein1 = row['protein1'].replace('{}.'.format(str(taxon)), '')
            protein2 = row['protein2'].replace('{}.'.format(str(taxon)), '')

            gene1_curie = None
            gene2_curie = None

            try:
                # Keep orientation the same since RO:0002434 is symmetric
                if protein1 > protein1:
                    gene1_curie = p2gene_map[protein1]
                    gene2_curie = p2gene_map[protein2]
                else:
                    gene1_curie = p2gene_map[protein2]
                    gene2_curie = p2gene_map[protein1]
            except KeyError:
                filtered_out_count += 1

            if gene1_curie is not None and gene2_curie is not None:
                # RO:0002434 ! interacts_with
                interacts_with = 'RO:0002434'
                self.graph.addTriple(gene1_curie, interacts_with, gene2_curie)
                if limit is not None and index >= limit:
                    break

        logger.info("Finished parsing p-p interactions for {},"
                    " {} rows filtered out based on checking"
                    " ensembl proteins".format(taxon, filtered_out_count))
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
                'file': "{}.{}".format(str(taxon), self.files[file_type]['pattern']),
                'url': "{}{}.{}".format(
                    self.files[file_type]['path'], str(taxon),
                    self.files[file_type]['pattern'])
            }
        return file_paths
