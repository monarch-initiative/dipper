import logging
import csv
import re

from dipper.sources.OMIMSource import OMIMSource
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class HGNC(OMIMSource):
    """
    This is the processing module for HGNC.

    - HGNC is authoratatitive for human gene symbols.
    - We adopt HGNC's reporting of
        - OMIM, Ensembl & Entrez equivalences of HGNC identifiers.

    - We include HGNC's reporting of publications "about" a human gene.

    - We also add the links to cytogenic locations for the gene features.
        - although HGNC is not authoratatitive for genomic locations

    """

    EBIFTP = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/'
    files = {
        'genes': {
            'file': 'hgnc_complete_set.txt',
            'url': EBIFTP + 'new/tsv/hgnc_complete_set.txt',
            'columns': [
                'hgnc_id',
                'symbol',
                'name',
                'locus_group',
                'locus_type',
                'status',
                'location',
                'location_sortable',
                'alias_symbol',
                'alias_name',
                'prev_symbol',
                'prev_name',
                'gene_family',
                'gene_family_id',
                'date_approved_reserved',
                'date_symbol_changed',
                'date_name_changed',
                'date_modified',
                'entrez_id',
                'ensembl_gene_id',
                'vega_id',
                'ucsc_id',
                'ena',
                'refseq_accession',
                'ccds_id',
                'uniprot_ids',
                'pubmed_id',
                'mgd_id',
                'rgd_id',
                'lsdb',
                'cosmic',
                'omim_id',
                'mirbase',
                'homeodb',
                'snornabase',
                'bioparadigms_slc',
                'orphanet',
                'pseudogene.org',
                'horde_id',
                'merops',
                'imgt',
                'iuphar',
                'kznf_gene_catalog',
                'mamit-trnadb',
                'cd',
                'lncrnadb',
                'enzyme_id',
                'intermediate_filament_db',
                'rna_central_ids',
                'lncipedia',
                'gtrnadb'
            ],
        },
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None,
                 tax_ids=None, gene_ids=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skolemized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='hgnc',
            ingest_title='HGNC',
            ingest_url='https://www.genenames.org/',
            ingest_logo='source-hgnc.png',
            license_url=None,
            data_rights='ftp://ftp.ebi.ac.uk/pub/databases/genenames/README.txt',
            # file_handle=None
        )

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids

        self.gene_ids = []

        if 'gene' not in self.all_test_ids:
            LOG.warning("not configured with gene test ids.")
        else:
            self.gene_ids = self.all_test_ids['gene']
        self.hs_txid = self.globaltt['Homo sapiens']

        # to help detect obsolete usages in other ingests (someday)
        self.withdrawn = {}

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)
        if self.test_only:
            self.test_mode = True
        LOG.info("Parsing files...")
        self._process_genes(limit)
        LOG.info("Done parsing files.")

    def _process_genes(self, limit=None):

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        geno = Genotype(graph)
        model = Model(graph)

        raw = '/'.join((self.rawdir, self.files['genes']['file']))
        col = self.files['genes']['columns']
        LOG.info("Processing HGNC genes")

        chr_pattern = re.compile(r'(\d+|X|Y|Z|W|MT)[pq$]')
        band_pattern = re.compile(r'([pq][A-H\d]?\d?(?:\.\d+)?)')

        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')

            row = next(filereader)
            if not self.check_fileheader(col, row):
                pass

            for row in filereader:
                # To generate:
                # head -1 hgnc_complete_set.txt.1 | tr '\t' '\n' |
                # sed "s/\(.*\)/\1 = row[col.index(\'\1\')]/g"

                hgnc_id = row[col.index('hgnc_id')].strip()
                symbol = row[col.index('symbol')].strip()
                name = row[col.index('name')].strip()
                # locus_group = row[col.index('locus_group')]
                locus_type = row[col.index('locus_type')].strip()
                status = row[col.index('status')].strip()
                # 41622 Approved  & 1752 Entry Withdrawn
                location = row[col.index('location')].strip()
                # location_sortable = row[col.index('location_sortable')]
                # alias_symbol = row[col.index('alias_symbol')]
                # alias_name = row[col.index('alias_name')]
                # prev_symbol = row[col.index('prev_symbol')]
                # prev_name = row[col.index('prev_name')]
                # gene_family = row[col.index('gene_family')]
                # gene_family_id = row[col.index('gene_family_id')]
                # date_approved_reserved = row[col.index('date_approved_reserved')]
                # date_symbol_changed = row[col.index('date_symbol_changed')]
                # date_name_changed = row[col.index('date_name_changed')]
                # date_modified = row[col.index('date_modified')]
                entrez_id = row[col.index('entrez_id')].strip()
                ensembl_gene_id = row[col.index('ensembl_gene_id')].strip()
                # vega_id = row[col.index('vega_id')]
                # ucsc_id = row[col.index('ucsc_id')]
                # ena = row[col.index('ena')]
                # refseq_accession = row[col.index('refseq_accession')]
                # ccds_id = row[col.index('ccds_id')]
                # uniprot_ids = row[col.index('uniprot_ids')]
                pubmed_ids = row[col.index('pubmed_id')].strip()  # pipe separated!
                # mgd_id = row[col.index('mgd_id')]
                # rgd_id = row[col.index('rgd_id')]
                # lsdb = row[col.index('lsdb')]
                # cosmic = row[col.index('cosmic')]
                omim_ids = row[col.index('omim_id')].strip()  # pipe separated!
                # mirbase = row[col.index('mirbase')]
                # homeodb = row[col.index('homeodb')]
                # snornabase = row[col.index('snornabase')]
                # bioparadigms_slc = row[col.index('bioparadigms_slc')]
                # orphanet = row[col.index('orphanet')]
                # pseudogene.org = row[col.index('pseudogene.org')]
                # horde_id = row[col.index('horde_id')]
                # merops = row[col.index('merops')]
                # imgt = row[col.index('imgt')]
                # iuphar = row[col.index('iuphar')]
                # kznf_gene_catalog = row[col.index('kznf_gene_catalog')]
                # mamit_trnadb = row[col.index('mamit-trnadb')]
                # cd = row[col.index('cd')]
                # lncrnadb = row[col.index('lncrnadb')]
                # enzyme_id = row[col.index('enzyme_id')]
                # intermediate_filament_db = row[col.index('intermediate_filament_db')]
                # rna_central_ids = row[col.index('rna_central_ids')]
                # lncipedia = row[col.index('lncipedia')]
                # gtrnadb = row[col.index('gtrnadb')]

                if status != 'Approved':
                    self.withdrawn[hgnc_id]=symbol
                    continue

                if (self.test_mode and entrez_id != '' and
                        entrez_id not in self.gene_ids):
                    continue

                if name == '':
                    name = None

                if locus_type == 'withdrawn':
                    model.addDeprecatedClass(hgnc_id,
                                             old_id_category=blv.terms.Gene.value)
                elif symbol[-1] == '@':  # 10)  region (HOX), RNA cluster, gene (PCDH)
                    continue

                else:
                    gene_type_id = self.resolve(locus_type, mandatory=False)
                    if gene_type_id != locus_type:
                        model.addClassToGraph(hgnc_id, symbol, gene_type_id, name,
                                              class_category=blv.terms.Gene.value)
                    model.makeLeader(hgnc_id,
                                     node_category=blv.terms.Gene.value)

                if entrez_id != '':
                    model.addEquivalentClass(hgnc_id, 'NCBIGene:' + entrez_id,
                                             subject_category=blv.terms.Gene.value,
                                             object_category=blv.terms.Gene.value)

                if ensembl_gene_id != '':
                    model.addEquivalentClass(hgnc_id, 'ENSEMBL:' + ensembl_gene_id,
                                             subject_category=blv.terms.Gene.value,
                                             object_category=blv.terms.Gene.value)

                for omim_id in omim_ids.split('|'):
                    if omim_id in self.omim_replaced:
                        repl = self.omim_replaced[omim_id]
                        LOG.warning('%s is replaced with %s', omim_id, repl)
                        for omim in repl:
                            if self.omim_type[omim] == self.globaltt['gene']:
                                omim_id = omim

                    if omim_id in self.omim_type and \
                            self.omim_type[omim_id] == self.globaltt['gene']:
                        model.addEquivalentClass(hgnc_id, 'OMIM:' + omim_id,
                                                 subject_category=blv.terms.Gene.value,
                                                 object_category=blv.terms.Gene.value)

                geno.addTaxon(self.hs_txid, hgnc_id,
                              genopart_category=blv.terms.Gene.value)

                # add pubs as "is about"
                for pubmed_id in pubmed_ids.split('|'):
                    graph.addTriple(
                        'PMID:' + pubmed_id, self.globaltt['is_about'], hgnc_id,
                        subject_category=blv.terms.Publication.value,
                        object_category=blv.terms.Gene.value)

                # add chr location
                # sometimes two are listed, like: 10p11.2 or 17q25
                # -- there are only 2 of these FRA10A and MPFD
                # sometimes listed like "1 not on reference assembly"
                # sometimes listed like 10q24.1-q24.3
                # sometimes like 11q11 alternate reference locus
                band = chrom = None
                chr_match = chr_pattern.match(location)
                if chr_match is not None and chr_match.groups():
                    chrom = chr_match.group(1)
                    chrom_id = makeChromID(chrom, self.hs_txid, 'CHR')
                    band_match = band_pattern.search(location)
                    feat = Feature(graph, hgnc_id, None, None,
                                   feature_category=blv.terms.Gene.value)
                    if band_match is not None and band_match.groups():
                        band = band_match.group(1)
                        band = chrom + band
                        # add the chr band as the parent to this gene
                        # as a feature but assume that the band is created
                        # as a class with properties elsewhere in Monochrom
                        band_id = makeChromID(band, self.hs_txid, 'CHR')
                        model.addClassToGraph(band_id, None,
                                              class_category=
                                              blv.terms.GenomicSequenceLocalization.value)
                        feat.addSubsequenceOfFeature(band_id)
                    else:
                        model.addClassToGraph(chrom_id, None,
                                              blv.terms.GenomicEntity.value)
                        feat.addSubsequenceOfFeature(chrom_id)

                if not self.test_mode and limit is not None and \
                        filereader.line_num > limit:
                    break

    def getTestSuite(self):
        import unittest
        from tests.test_hgnc import HGNCTestCase
        test_suite = unittest.TestLoader().loadTestsFromTestCase(HGNCTestCase)
        return test_suite
