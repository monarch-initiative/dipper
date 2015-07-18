import logging
import csv
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Genotype import Genotype
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper import config
from dipper.models.GenomicFeature import Feature, makeChromID


logger = logging.getLogger(__name__)


class HGNC(Source):
    """
    This is the processing module for HGNC.

    We create equivalences between HGNC identifiers and ENSEMBL and NCBIGene.
    We also add the links to cytogenic locations for the gene features.

    """

    files = {
        'genes': {'file': 'hgnc_complete_set.txt',
                  'url': 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt'},
    }

    def __init__(self, tax_ids=None, gene_ids=None):
        Source.__init__(self, 'hgnc')

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids
        self.load_bindings()

        self.dataset = Dataset('hgnc', 'HGNC', 'http://www.genenames.org', None)

        self.gene_ids = []
        if 'test_ids' not in config.get_config() or 'gene' not in config.get_config()['test_ids']:
            logger.warn("not configured with gene test ids.")
        else:
            self.gene_ids = config.get_config()['test_ids']['gene']

        self.properties = Feature.properties

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        if self.testOnly:
            self.testMode = True

        logger.info("Parsing files...")

        self._process_genes(limit)

        self.load_core_bindings()
        self.load_bindings()

        logger.info("Done parsing files.")

        logger.info("Found %d nodes in graph", len(self.graph))
        logger.info("Found %d nodes in testgraph", len(self.testgraph))

        return

    def _process_genes(self, limit=None):
        gu = GraphUtils(curie_map.get())

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)

        raw = '/'.join((self.rawdir, self.files['genes']['file']))
        line_counter = 0
        logger.info("Processing HGNC genes")

        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (hgnc_id, symbol, name, locus_group, locus_type, status, location, location_sortable,
                 alias_symbol, alias_name, prev_symbol, prev_name, gene_family, gene_family_id,
                 date_approved_reserved, date_symbol_changed, date_name_changed, date_modified,
                 entrez_id, ensembl_gene_id, vega_id, ucsc_id, ena, refseq_accession, ccds_id,
                 uniprot_ids, pubmed_id, mgd_id, rgd_id, lsdb, cosmic, omim_id, mirbase, homeodb,
                 snornabase, bioparadigms_slc, orphanet, pseudogene_org, horde_id, merops, imgt,
                 iuphar, kznf_gene_catalog, mamit_trnadb, cd, lncrnadb, enzyme_id,
                 intermediate_filament_db) = row

                line_counter += 1

                # skip header
                if line_counter <= 1:
                    continue

                if self.testMode and entrez_id != '' and int(entrez_id) not in self.gene_ids:
                    continue

                if name == '':
                    name = None
                gene_type_id = self._get_gene_type(locus_type)
                gu.addClassToGraph(g, hgnc_id, symbol, gene_type_id, name)
                if locus_type == 'withdrawn':
                    gu.addDeprecatedClass(g, hgnc_id)
                if entrez_id != '':
                    gu.addEquivalentClass(g, hgnc_id, 'NCBIGene:' + entrez_id)
                if ensembl_gene_id != '':
                    gu.addEquivalentClass(g, hgnc_id, 'ENSEMBL:' + ensembl_gene_id)
                geno.addTaxon('NCBITaxon:9606', hgnc_id)

                # add pubs as "is about"
                if pubmed_id != '':
                    for p in re.split('\|', pubmed_id.strip()):
                        if str(p) != '':
                            gu.addTriple(g, 'PMID:' + str(p.strip()), gu.object_properties['is_about'], hgnc_id)

                # add chr location
                # sometimes two are listed, like: 10p11.2 or 17q25  -- there are only 2 of these FRA10A and MPFD
                # sometimes listed like "1 not on reference assembly"
                # sometimes listed like 10q24.1-q24.3
                # sometimes like 11q11 alternate reference locus
                band = chrom = None
                chr_pattern = '(\d+|X|Y|Z|W|MT)[pq$]'
                chr_match = re.match(chr_pattern, location)
                if chr_match is not None and len(chr_match.groups()) > 0:
                    chrom = chr_match.group(1)
                    chrom_id = makeChromID(chrom, 'NCBITaxon:9606')
                    band_pattern = '([pq][A-H\d]?\d?(?:\.\d+)?)'
                    band_match = re.search(band_pattern, location)
                    f = Feature(hgnc_id, None, None)
                    if band_match is not None and len(band_match.groups()) > 0:
                        band = band_match.group(1)
                        band = chrom + band
                        # add the chr band as the parent to this gene as a feature
                        # but assume that the band is created as a class with properties elsewhere in Monochrom
                        band_id = makeChromID(band, 'NCBITaxon:9606')
                        gu.addClassToGraph(g, band_id, None)
                        f.addSubsequenceOfFeature(g, band_id)
                    else:
                        gu.addClassToGraph(g, chrom_id, None)
                        f.addSubsequenceOfFeature(g, chrom_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

            # end loop through file

        gu.loadProperties(g, Feature.object_properties, gu.OBJPROP)
        gu.loadProperties(g, Feature.data_properties, gu.DATAPROP)
        gu.loadProperties(g, Genotype.object_properties, gu.OBJPROP)
        gu.loadAllProperties(g)

        return

    @staticmethod
    def _get_gene_type(locus_type):
        """
        Given the locus_type string supplied by HGNC, we map them to relevant SO terms.  We use the
        mappings listed in the "?" popup on any HGNC gene page.
        :param locus_type:
        :return:
        """

        locus_type_map = {
            'RNA, Y': 'SO:0000405',
            'RNA, cluster': 'SO:0000655',  # ??? ncRNA
            'RNA, long non-coding': 'SO:0001877',  # lncRNA
            'RNA, micro': 'SO:0001265',  # miRNA
            'RNA, misc': 'SO:0000655',  # ??? ncRNA
            'RNA, ribosomal': 'SO:0001637',  # rRNA
            'RNA, small cytoplasmic': 'SO:0001266',  # small cytoplasmic is synonym of scRNA_primary_transcript
            'RNA, small nuclear': 'SO:0001268',  # snRNA
            'RNA, small nucleolar': 'SO:0001267',  # snoRNA
            'RNA, transfer': 'SO:0001272',  # tRNA
            'RNA, vault': 'SO:0000404',  # vault_RNA
            'T cell receptor gene': 'SO:0000460',  # vertebrate_immunoglobulin_T_cell_receptor_segment
            'T cell receptor pseudogene': 'SO:0000336',  # pseudogene
            'complex locus constituent': 'SO:0001217',  # protein_coding - no way to say member of cluster
            'endogenous retrovirus': 'SO:0000100',  # endogenous_retroviral_gene
            'fragile site': 'SO:0001411',  # biological_region TODO https://sourceforge.net/p/song/term-tracker/433/
            'gene with protein product': 'SO:0001217',  # protein_coding_gene
            'immunoglobulin gene': 'SO:0000460',  # immunoglobulin_region
            'immunoglobulin pseudogene': 'SO:0000336',  # pseudogene  FIXME
            'phenotype only': 'SO:0001500',  # heritable_phenotypic_marker
            'protocadherin': 'SO:0000110',  # FIXME
            'pseudogene': 'SO:0000336',  # pseudogene
            'readthrough': 'SO:0000110',  # FIXME
            'region': 'SO:0000001',  # region
            'transposable element': 'SO:0000101',  # transposable element
            'unknown': 'SO:0000110',  # sequence_feature
            'virus integration site': 'SO:0000110',  # FIXME
            'withdrawn': None  # TODO
        }
        t = None
        if locus_type in locus_type_map:
            t = locus_type_map[locus_type]
        else:
            logger.error("Unknown/unmapped locus type: %s", locus_type)

        return t

    def getTestSuite(self):
        import unittest
        from tests.test_hgnc import HGNCTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(HGNCTestCase)

        return test_suite