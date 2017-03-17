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
    STRING_VERSION = "v10"
    DEFAULT_TAXA = [9606, 10090, 7955, 7227, 6239]

    files = {
        'protein_links': {
            'path': '{}protein.links.detailed.{}/'.format(STRING_BASE, STRING_VERSION),
            'pattern': 'protein.links.detailed.{}.txt.gz'.format(STRING_VERSION)
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'string')
        self.dataset = Dataset(
            'string', 'String', 'http://string-db.org/', None,
            'http://string-db.org/cgi/access.pl?footer_active_subpage=licensing')

        if tax_ids is None:
            self.tax_ids = StringDB.DEFAULT_TAXA
        else:
            logger.info("Filtering on taxa {}".format(tax_ids))
            self.tax_ids = tax_ids

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
        ensembl = Ensembl(self.graph_type, self.are_bnodes_skized)

        for taxon in protein_paths:
            string_file_path = '/'.join((
                self.rawdir, protein_paths[taxon]['file']))

            fh = gzip.open(string_file_path, 'rb')
            dataframe = pd.read_csv(fh, sep='\s+')
            protein_list = ensembl.fetch_protein_list(taxon)

            logger.info("Fetching protein protein interactions "
                        "for taxon {}".format(taxon))

            self._process_protein_links(dataframe, protein_list, taxon, limit)

    def _process_protein_links(self, dataframe, prot_list, taxon,
                               limit=None, rank_min=10):
        for index, row in dataframe.iterrows():
            # Check if proteins are in same species
            protein1 = row['protein1'].replace('{}.'.format(str(taxon)), '')
            protein2 = row['protein2'].replace('{}.'.format(str(taxon)), '')
            filtered_out_count = 0

            if protein1 not in prot_list \
                    or protein2 not in prot_list:
                filtered_out_count += 1
            elif row['experimental'] > rank_min:
                protein1_curie = "ENSEMBL:{}".format(protein1)
                protein2_curie = "ENSEMBL:{}".format(protein2)
                # RO:0002434 ! interacts_with
                interacts_with = 'RO:0002434'
                self.graph.addTriple(protein1_curie, interacts_with, protein2_curie)
                # Until we have symmetry working
                self.graph.addTriple(protein2_curie, interacts_with, protein1_curie)
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
