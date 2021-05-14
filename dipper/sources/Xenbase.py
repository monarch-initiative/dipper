import csv
import logging
from typing import Dict, List

from dipper.sources.Source import Source
from dipper.models.Genotype import Genotype
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Model import Model


LOG = logging.getLogger(__name__)

DIPPER_CACHE = 'https://archive.monarchinitiative.org/DipperCache/xenbase/'
XENBASE_FTP = 'http://ftp.xenbase.org/'

OBAN_COLS = [
    'SUBJECT',
    'SUBJECT_LABEL',
    'SUBJECT_TAXON',
    'SUBJECT_TAXON_LABEL',
    'OBJECT',
    'OBJECT_LABEL',
    'RELATION',
    'RELATION_LABEL',
    'EVIDENCE',
    'EVIDENCE_LABEL',
    'SOURCE',
    'IS_DEFINED_BY',
    'QUALIFIER',
]


class Xenbase(Source):
    """
    Xenbase is a web-accessible resource that integrates all the diverse biological,
    genomic, genotype and phenotype data available from Xenopus research.
    """

    files = {
        'g2p_assertions': {
            'file': 'xb_xpo_spo_v20210511b.csv',
            'url': DIPPER_CACHE + 'xb_xpo_spo_v20210511b.csv',
            'columns': OBAN_COLS
        },
        'orthologs': {  # Not ingesting this yet since it appears to be from another source
                        # and not custom to Xenbase
            'file': 'xb_ortho_spo_v_v20210318a.csv',
            'url': DIPPER_CACHE + 'xb_ortho_spo_v_v20210318a.csv',
            'columns': OBAN_COLS
        },
        'gene_literature': {
            'file': 'LiteratureMatchedGenesByPaper.txt',
            'url': XENBASE_FTP + '/pub/GenePageReports/LiteratureMatchedGenesByPaper.txt',
            'columns': [
                'xb_article',
                'pmid',
                'gene_pages'
            ]
        },
        'genepage2gene': {
            'file': 'XenbaseGenepageToGeneIdMapping.txt',
            'url': XENBASE_FTP + '/pub/GenePageReports/XenbaseGenepageToGeneIdMapping.txt',
            'columns': [
                'gene_page_id',
                'gene_page_label',
                'tropicalis_id',
                'tropicalis_label',
                'laevis_l_id',
                'laevis_l_label',
                'laevis_s_id',
                'laevis_s_label',
            ]
        }
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='xenbase',
            ingest_title='Xenbase: The Xenopus model organism knowledgebase',
            ingest_url='http://www.xenbase.org',
            ingest_logo = 'source-xenbase.png',
            license_url=None,
            data_rights='http://www.xenbase.org/other/static-xenbase/citingMOD.jsp',
            file_handle=None
        )

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows fo each file", str(limit))

        LOG.info("Parsing files...")

        self._parse_g2p_file(limit)
        genepage2gene = self._parse_genepage2gene(limit)
        self._parse_gene_literature(limit, genepage2gene)

    def _parse_g2p_file(self, limit=None):
        """
        Parse gene to XPO file, currently custom for Monarch
        :param limit:
        :return:
        """
        src_key = 'g2p_assertions'
        geno = Genotype(self.graph)
        model = Model(self.graph)

        columns = self.files[src_key]['columns']
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))

        LOG.info("Processing Gene to XPO associations")

        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile)

            # File has headers
            row = next(reader)
            if not self.check_fileheader(columns, row):
                pass

            for row in reader:

                gene = row[columns.index('SUBJECT')]
                gene_label = row[columns.index('SUBJECT_LABEL')]
                gene_taxon = row[columns.index('SUBJECT_TAXON')]
                #gene_taxon_label = row[columns.index('SUBJECT_TAXON_LABEL')]
                phenotype_curie = row[columns.index('OBJECT')]
                #phenotype_label = row[columns.index('OBJECT_LABEL')]
                relation = row[columns.index('RELATION')]
                #relation_label = row[columns.index('RELATION_LABEL')]
                evidence = row[columns.index('EVIDENCE')]
                #evidence_label = row[columns.index('EVIDENCE_LABEL')]
                source = row[columns.index('SOURCE')]
                #is_defined_by = row[columns.index('IS_DEFINED_BY')]
                #qualifier = row[columns.index('QUALIFIER')]

                gene_curie = 'Xenbase:' + gene
                relation_curie = relation.replace('_', ':')

                geno.addGene(gene_curie, gene_label)
                geno.addTaxon(gene_taxon, gene_curie)

                assoc = G2PAssoc(
                    self.graph,
                    self.name,
                    entity_id=gene_curie,
                    phenotype_id=phenotype_curie,
                    rel=relation_curie
                )

                if evidence:
                    assoc.add_evidence(evidence)

                if source:
                    model.addType(source, self.globaltt['journal article'])
                    assoc.add_source(source)

                assoc.add_association_to_graph()

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

    def _parse_genepage2gene(self, limit) -> Dict[str, List[str]]:
        """
        :return:
        """
        src_key = 'genepage2gene'
        columns = self.files[src_key]['columns']
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))

        geno = Genotype(self.graph)
        genepage2gene = {}

        LOG.info("Processing GenePage to Gene file")

        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')

            for row in reader:

                gene_page = row[columns.index('gene_page_id')]
                # gene_page_label = row[columns.index('gene_page_label')]
                tropicalis_id = row[columns.index('tropicalis_id')]
                tropicalis_label = row[columns.index('tropicalis_label')]
                laevis_l_id = row[columns.index('laevis_l_id')]
                laevis_l_label = row[columns.index('laevis_l_label')]
                laevis_s_id = row[columns.index('laevis_s_id')]
                laevis_s_label = row[columns.index('laevis_s_label')]

                tropicalis_curie = 'Xenbase:' + tropicalis_id
                laevis_l_curie = 'Xenbase:' + laevis_l_id
                laevis_s_curie = 'Xenbase:' + laevis_s_id

                genepage2gene[gene_page] = [tropicalis_curie, laevis_l_curie, laevis_s_curie]

                geno.addGene(tropicalis_curie, tropicalis_label)
                geno.addGene(laevis_l_curie, laevis_l_label)
                geno.addGene(laevis_s_curie, laevis_s_label)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        return genepage2gene

    def _parse_gene_literature(self, limit, genepage2gene):
        """
        :return:
        """
        src_key = 'gene_literature'
        columns = self.files[src_key]['columns']
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))

        model = Model(self.graph)

        LOG.info("Processing gene page to literature")

        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')

            for row in reader:

                gene_pages = row[columns.index('gene_pages')]
                pmid = row[columns.index('pmid')]
                # xb_article = row[columns.index('xb_article')]

                pmid_curie = 'PMID:' + pmid

                for gene_page in gene_pages.split(','):
                    gene_page_id = gene_page.split(' ')[0]
                    try:
                        gene_ids = genepage2gene[gene_page_id]
                    except KeyError:
                        LOG.debug("Could not locate genepage_id: %s in row %s", gene_page_id, row)
                        continue
                    for gene in gene_ids:
                        model.addTriple(pmid_curie, self.globaltt['mentions'], gene)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break
