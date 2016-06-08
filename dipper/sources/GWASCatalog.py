import logging
import csv
import re

from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Dataset import Dataset
from dipper import config
from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper.models.Genotype import Genotype
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Reference import Reference
from dipper.models.GenomicFeature import Feature, makeChromID

logger = logging.getLogger(__name__)


class GWASCatalog(Source):
    """
    The NHGRI-EBI Catalog of published genome-wide association studies.

    We link the variants recorded here to the curated EFO-classes using a
    "contributes_to" linkage because the only thing we know is that the SNPs
    are associated with the trait/disease,
    but we don't know if it is actually causative.

    Description of the GWAS catalog is here:
    http://www.ebi.ac.uk/gwas/docs/fileheaders#_file_headers_for_catalog_version_1_0_1

    Status:  IN PROGRESS

    """

    terms = {
        'cell_line_repository': 'CLO:0000008',
        'race': 'SIO:001015',
        'ethnic_group': 'EFO:0001799',
        'age': 'EFO:0000246',
        'sampling_time': 'EFO:0000689',
        'collection': 'ERO:0002190'
    }

    files = {
        'catalog': {
            'file': 'gwas_catalog.tsv',
            'url': 'http://www.ebi.ac.uk/gwas/api/search/downloads/alternative'}
    }

    def __init__(self):
        Source.__init__(self, 'gwascatalog')

        self.load_bindings()

        self.dataset = Dataset(
            'gwascatalog', 'GWAS Catalog', 'http://www.ebi.ac.uk/gwas/',
            'The NHGRI-EBI Catalog of published genome-wide association studies',
            'http://creativecommons.org/licenses/by/3.0/', None)
        # 'http://www.ebi.ac.uk/gwas/docs/about'  # TODO add this

        if 'test_ids' not in config.get_config() or \
                'gene' not in config.get_config()['test_ids']:
            logger.warning("not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']

        return

    def fetch(self, is_dl_forced=False):
        """

        :param is_dl_forced:
        :return:
        """
        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self.process_catalog(limit)

        logger.info("Finished parsing.")

        self.load_bindings()

        logger.info("Found %d nodes in graph", len(self.graph))
        logger.info("Found %d nodes in testgraph", len(self.testgraph))

        return

    def process_catalog(self, limit=None):
        """
        :param limit:
        :return:

        """
        raw = '/'.join((self.rawdir, self.files['catalog']['file']))
        logger.info("Processing Data from %s", raw)
        gu = GraphUtils(curie_map.get())

        if self.testMode:      # set the graph to build
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        geno = Genotype(g)

        gu.loadProperties(g, geno.object_properties, gu.OBJPROP)
        gu.loadAllProperties(g)

        tax_id = 'NCBITaxon:9606'  # hardcode
        genome_version = 'GRCh38'  # hardcode

        # build a hashmap of genomic location to identifiers,
        # to try to get the equivalences

        loc_to_id_hash = {}

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                if not row:
                    pass
                else:
                    line_counter += 1
                    (date_added_to_catalog, pubmed_num, first_author, pub_date,
                     journal, link, study_name, disease_or_trait,
                     initial_sample_description, replicate_sample_description,
                     region, chrom_num, chrom_pos, reported_gene_nums,
                     mapped_gene, upstream_gene_num, downstream_gene_num,
                     snp_gene_nums, upstream_gene_distance,
                     downstream_gene_distance, strongest_snp_risk_allele, snps,
                     merged, snp_id_current, context, intergenic_flag,
                     risk_allele_frequency, pvalue, pvalue_mlog, pvalue_text,
                     or_or_beta, confidence_interval_95,
                     platform_with_snps_passing_qc, cnv_flag, mapped_trait,
                     mapped_trait_uri) = row

                    intersect = \
                        list(set([str(i) for i in self.test_ids['gene']]) &
                             set(re.split(r',', snp_gene_nums)))
                    # skip if no matches found in test set
                    if self.testMode and len(intersect) == 0:
                        continue

# 06-May-2015	25917933	Zai CC	20-Nov-2014	J Psychiatr Res	http://europepmc.org/abstract/MED/25917933
# A genome-wide association study of suicide severity scores in bipolar disorder.
# Suicide in bipolar disorder
# 959 European ancestry individuals	NA
# 10p11.22	10	32704340	C10orf68, CCDC7, ITGB1	CCDC7
# rs7079041-A	rs7079041	0	7079041	intron	0		2E-6	5.698970
                    if chrom_num != '' and chrom_pos != '':
                        loc = 'chr'+str(chrom_num)+':'+str(chrom_pos)
                        if loc not in loc_to_id_hash:
                            loc_to_id_hash[loc] = set()
                    else:
                        loc = None

                    if re.search(r' x ', strongest_snp_risk_allele) \
                            or re.search(r',', strongest_snp_risk_allele):
                        # TODO deal with haplotypes
                        logger.warning(
                            "We can't deal with haplotypes yet: %s",
                            strongest_snp_risk_allele)
                        continue
                    elif re.match(r'rs', strongest_snp_risk_allele):
                        rs_id = 'dbSNP:'+strongest_snp_risk_allele.strip()
                        # remove the alteration
                    elif re.match(r'kgp', strongest_snp_risk_allele):
                        # FIXME this isn't correct
                        rs_id = 'dbSNP:'+strongest_snp_risk_allele.strip()
                        # http://www.1000genomes.org/faq/what-are-kgp-identifiers
                        # for some information
                        # They were created by Illumina for their genotyping
                        # platform before some variants identified during the
                        # pilot phase of the project had been assigned
                        # rs numbers.
                    elif re.match(r'chr', strongest_snp_risk_allele):
                        # like: chr10:106180121-G
                        rs_id = ':gwas-' + \
                            re.sub(
                                r':', '-', strongest_snp_risk_allele.strip())
                    elif strongest_snp_risk_allele.strip() == '':
                        # logger.debug(
                        #    "No strongest SNP risk allele for %s:\n%s",
                        #    pubmed_num, str(row))
                        # FIXME still consider adding in the EFO terms
                        # for what the study measured?
                        continue
                    else:
                        logger.warning(
                            "There's a snp id i can't manage: %s",
                            strongest_snp_risk_allele)
                        continue

                    alteration = re.search(r'-(.*)$', rs_id)
                    if alteration is not None \
                            and re.match(r'[ATGC]', alteration.group(1)):
                        # add variation to snp
                        pass  # TODO
                    rs_id = re.sub(r'-.*$', '', rs_id).strip()
                    if loc is not None:
                        loc_to_id_hash[loc].add(rs_id)

                    pubmed_id = 'PMID:'+pubmed_num

                    r = Reference(
                        pubmed_id, Reference.ref_types['journal_article'])
                    r.addRefToGraph(g)

                    # create the chromosome
                    chrom_id = makeChromID(chrom_num, genome_version, 'CHR')

                    # add the feature to the graph
                    snp_description = None
                    if risk_allele_frequency != '' and \
                            risk_allele_frequency != 'NR':
                        snp_description = \
                            str(risk_allele_frequency) + \
                            ' [risk allele frequency]'

                    f = Feature(
                        rs_id, strongest_snp_risk_allele.strip(),
                        Feature.types[r'SNP'], snp_description)
                    if chrom_num != '' and chrom_pos != '':
                        f.addFeatureStartLocation(chrom_pos, chrom_id)
                        f.addFeatureEndLocation(chrom_pos, chrom_id)
                    f.addFeatureToGraph(g)
                    f.addTaxonToFeature(g, tax_id)
                    # TODO consider adding allele frequency as property;
                    # but would need background info to do that

                    # also want to add other descriptive info about
                    # the variant from the context
                    for c in re.split(r';', context):
                        cid = self._map_variant_type(c.strip())
                        if cid is not None:
                            gu.addType(g, rs_id, cid)

                    # add deprecation information
                    if merged == 1 and str(snp_id_current.strip()) != '':
                        # get the current rs_id
                        current_rs_id = 'dbSNP:'
                        if not re.match(r'rs', snp_id_current):
                            current_rs_id += 'rs'
                        if loc is not None:
                            loc_to_id_hash[loc].append(current_rs_id)
                        current_rs_id += str(snp_id_current)
                        gu.addDeprecatedIndividual(g, rs_id, current_rs_id)
                        # TODO check on this
                        # should we add the annotations to the current
                        # or orig?
                        gu.makeLeader(g, current_rs_id)
                    else:
                        gu.makeLeader(g, rs_id)

                    # add the feature as a sequence alteration
                    # affecting various genes
                    # note that intronic variations don't necessarily list
                    # the genes such as for rs10448080  FIXME
                    if snp_gene_nums != '':
                        for s in re.split(r',', snp_gene_nums):
                            s = s.strip()
                            # still have to test for this,
                            # because sometimes there's a leading comma
                            if s != '':
                                gene_id = 'NCBIGene:'+s
                                geno.addAlleleOfGene(rs_id, gene_id)

                    # add the up and downstream genes if they are available
                    if upstream_gene_num != '':
                        downstream_gene_id = 'NCBIGene:'+downstream_gene_num
                        gu.addTriple(
                            g, rs_id,
                            Feature.object_properties[
                                r'upstream_of_sequence_of'],
                            downstream_gene_id)
                    if downstream_gene_num != '':
                        upstream_gene_id = 'NCBIGene:'+upstream_gene_num
                        gu.addTriple(
                            g, rs_id,
                            Feature.object_properties[
                                'downstream_of_sequence_of'],
                            upstream_gene_id)

                    description = 'A study of ' + disease_or_trait + \
                        ' in ' + initial_sample_description
                    if replicate_sample_description != '':
                        description = \
                            ' '.join(
                                (description, 'with',
                                 replicate_sample_description))
                    if platform_with_snps_passing_qc != '':
                        description = ' '.join(
                            (description, 'on platform',
                             platform_with_snps_passing_qc))
                    description = ' '.join((description, '(p='+pvalue+')'))

                    # make associations to the EFO terms; there can be >1
                    if mapped_trait_uri.strip() != '':
                        for t in re.split(r',', mapped_trait_uri):
                            t = t.strip()

                            cu = CurieUtil(curie_map.get())
                            tid = cu.get_curie(t)

                            assoc = G2PAssoc(
                                self.name, rs_id, tid,
                                gu.object_properties['contributes_to'])
                            assoc.add_source(pubmed_id)
                            # combinatorial evidence
                            # used in automatic assertion
                            eco_id = 'ECO:0000213'
                            assoc.add_evidence(eco_id)

                            # assoc.set_description(description)
                            # FIXME score should get added to provenance/study
                            # assoc.set_score(pvalue)
                            assoc.add_association_to_graph(g)
                    else:
                        print('foo')

                    if not self.testMode and\
                            (limit is not None and line_counter > limit):
                        break

            Assoc(self.name).load_all_properties(g)

        # loop through the location hash,
        # and make all snps at that location equivalent
        for l in loc_to_id_hash:
            snp_ids = loc_to_id_hash[l]
            if len(snp_ids) > 1:
                logger.info("%s has >1 snp id: %s", l, str(snp_ids))
        return

    @staticmethod
    def _map_variant_type(sample_type):
        ctype = None
        type_map = {
            'stop_gained': 'SO:0001587',       # stop-gain variant
            'intron_variant': 'SO:0001627',  # intron variant
            '3_prime_UTR_variant': 'SO:0001624',           # 3'utr variant
            '5_prime_UTR_variant': 'SO:0001623',           # 5'UTR variant
            'synonymous_variant': 'SO:0001819',       # synonymous variant
            'frameshift_variant': 'SO:0001589',      # frameshift
            'intergenic_variant': 'SO:0001628',     # intergenic_variant
            'non_coding_transcript_exon_variant': 'SO:0001619', # noncoding transcript variant
            'splice_acceptor_variant': 'SO:0001574',        # splice acceptor variant
            'splice_donor_variant': 'SO:0001575',        # splice donor variant
            'missense_variant': 'SO:0001583',       # missense variant
            'downstream_gene_variant': 'SO:0001634',      # 500B_downstream_variant
            'upstream_gene_variant': 'SO:0001636',      # 2KB_upstream_variant
            'coding_sequence_variant': 'SO:0001580',  #coding_sequence_variant
            'non_coding_exon_variant ': 'SO:0001792',
            'regulatory_region_variant': 'SO:0001566',
            'splice_region_variant': 'SO:0001630',
            'stop_lost': 'SO:0001578',
            'TF_binding_site_variant': 'SO:0001782'
        }
        if sample_type.strip() in type_map:
            ctype = type_map.get(sample_type)
        elif sample_type.strip() != '':
            logger.error("Variant type not mapped: %s", sample_type)

        return ctype

    def getTestSuite(self):
        import unittest
        from tests.test_gwascatalog import GWASCatalogTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(GWASCatalogTestCase)

        return test_suite
