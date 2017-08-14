from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Pathway import Pathway
from dipper.models.Dataset import Dataset
from ontobio.io.gafparser import GafParser
from pprint import pprint

import logging

__author__ = 'timputman'

logger = logging.getLogger(__name__)


class RGD(Source):
    RGD_BASE = 'ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/'
    files = {
        'rat_gene2mammalian_phenotype': {
            'file': 'rattus_genes_mp',
            'url': RGD_BASE + 'rattus_genes_mp'},
    }

    map_files = {
        'eco_map': 'http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt',
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'rat_genome_database')
        self.dataset = Dataset(
            'rat_genome_database', 'Rat_Genome_Database', 'http://rgd.mcw.edu/', None,
            None)

    def fetch(self, is_dl_forced=False):
        """
        Override Source.fetch()
        Fetches resources from rat_genome_database using the rat_genome_database ftp site
        Args:
            :param is_dl_forced (bool): Force download
        Returns:
            :return None
        """
        print(self.rawdir)

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

        rgd_file = '/'.join((self.rawdir, self.files['rat_gene2mammalian_phenotype']['file']))
        p = GafParser()
        assocs = p.parse(open(rgd_file, "r"))

        return


example_assoc = {'aspect': 'N',
                 'date': '2006-10-17',
                 'evidence': {'has_supporting_reference': ['RGD:1581602', 'PMID:16368876'],
                              'type': 'IAGP',
                              'with_support_from': []},
                 'negated': False,
                 'object': {'id': 'MP:0005211', 'taxon': 'NCBITaxon:10116'},
                 'provided_by': 'RGD',
                 'qualifiers': [],
                 'relation': {'id': None},
                 'source_line': 'RGD\t621503\tKcnq1\t\tMP:0005211\tRGD:1581602|PMID:16368876\t'
                                'IAGP\t\tN\tpotassium voltage-gated channel subfamily Q member '
                                '1\t\tgene\ttaxon:10116\t20061017\tRGD\t\t\n',
                 'subject': {'fullname': 'potassium voltage-gated channel subfamily Q member 1',
                             'id': 'RGD:621503',
                             'label': 'Kcnq1',
                             'synonyms': [],
                             'taxon': {'id': 'NCBITaxon:10116'},
                             'type': 'gene'},
                 'subject_extensions': [{'filler': '\n', 'property': 'isoform'}]}
