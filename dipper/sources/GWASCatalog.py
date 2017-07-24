import logging
import csv
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper import config
from dipper.utils.CurieUtil import CurieUtil
from dipper.utils.DipperUtil import DipperUtil
from dipper.models.Model import Model
from dipper import curie_map
from dipper.models.Genotype import Genotype
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Reference import Reference
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.graph.RDFGraph import RDFGraph

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

    GWAS also pulishes Owl files described here
    http://www.ebi.ac.uk/gwas/docs/ontology


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

    GWASFTP = 'ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest'
    GWASFILE = 'gwas-catalog-associations_ontology-annotated.tsv'
    files = {
        'catalog': {
            'file': GWASFILE,
            'url': GWASFTP + '/' + GWASFILE},
        'efo': {
            'file': 'efo.owl',
            'url': 'http://www.ebi.ac.uk/efo/efo.owl'},
        'so': {
            'file': 'so.owl',
            'url': 'http://purl.obolibrary.org/obo/so.owl'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'gwascatalog')

        if graph_type != 'rdf_graph':
            raise ValueError("UDP requires a rdf_graph")

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

        # build a dictionary of genomic location to identifiers,
        # to try to get the equivalences
        self.id_location_map = dict()

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
        return

    def process_catalog(self, limit=None):
        """
        :param limit:
        :return:

        """
        raw = '/'.join((self.rawdir, self.files['catalog']['file']))
        logger.info("Processing Data from %s", raw)
        efo_ontology = RDFGraph()
        logger.info("Loading EFO ontology in separate rdf graph")
        efo_ontology.parse(self.files['efo']['url'], format='xml')
        efo_ontology.bind_all_namespaces()
        logger.info("Finished loading EFO ontology")

        so_ontology = RDFGraph()
        logger.info("Loading SO ontology in separate rdf graph")
        so_ontology.parse(self.files['so']['url'], format='xml')
        so_ontology.bind_all_namespaces()
        logger.info("Finished loading SO ontology")

        line_counter = 0

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t')
            header = next(filereader, None)  # the header row
            header_len = len(header)
            logger.info('header length:\t %i', header_len)

            for row in filereader:
                if not row:
                    pass
                else:
                    line_counter += 1
                    if header_len != len(row):
                        logger.error(
                            'BadRow: %i has %i columns', line_counter, row)
                        pass
                    (date_added_to_catalog, pubmed_num, first_author,
                     pub_date, journal, link, study_name, disease_or_trait,
                     initial_sample_description, replicate_sample_description,
                     region, chrom_num, chrom_pos, reported_gene_nums,
                     mapped_gene, upstream_gene_num, downstream_gene_num,
                     snp_gene_nums, upstream_gene_distance,
                     downstream_gene_distance, strongest_snp_risk_allele, snps,
                     merged, snp_id_current, context, intergenic_flag,
                     risk_allele_frequency, pvalue, pvalue_mlog, pvalue_text,
                     or_or_beta, confidence_interval_95,
                     platform_with_snps_passing_qc, cnv_flag, mapped_trait,
                     mapped_trait_uri, study_accession) = row

                    intersect = list(
                        set([str(i) for i in self.test_ids['gene']]) &
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

                    variant_curie, variant_type = \
                        self._get_curie_and_type_from_id(
                            strongest_snp_risk_allele)

                    if strongest_snp_risk_allele.strip() == '':
                        logger.debug(
                            "No strongest SNP risk allele for %s:\n%s",
                            pubmed_num, str(row))
                        # still consider adding in the EFO terms
                        # for what the study measured?
                        continue

                    if variant_type == 'snp':
                        self._add_snp_to_graph(
                            variant_curie, strongest_snp_risk_allele,
                            chrom_num, chrom_pos,
                            context, risk_allele_frequency)

                        self._add_deprecated_snp(variant_curie, snp_id_current,
                                                 merged, chrom_num, chrom_pos)

                        self._add_snp_gene_relation(
                            variant_curie, snp_gene_nums, upstream_gene_num,
                            downstream_gene_num)
                    elif variant_type == 'haplotype':
                        self._process_haplotype(
                            variant_curie, strongest_snp_risk_allele,
                            chrom_num, chrom_pos, context,
                            risk_allele_frequency, mapped_gene, so_ontology)
                    elif variant_type is None:
                        logger.warning(
                            "There's a snp id i can't manage: %s",
                            strongest_snp_risk_allele)
                        continue

                    description = self._make_description(
                        disease_or_trait, initial_sample_description,
                        replicate_sample_description,
                        platform_with_snps_passing_qc, pvalue)

                    self._add_variant_trait_association(
                        variant_curie, mapped_trait_uri, efo_ontology,
                        pubmed_num, description)

                    if not self.testMode and\
                            (limit is not None and line_counter > limit):
                        break

        # TODO loop through the location hash,
        # and make all snps at that location equivalent
        for l in self.id_location_map:
            snp_ids = self.id_location_map[l]
            if len(snp_ids) > 1:
                logger.info("%s has >1 snp id: %s", l, str(snp_ids))
        return

    def _process_haplotype(
            self, hap_id, hap_label, chrom_num, chrom_pos, context,
            risk_allele_frequency, mapped_gene, so_ontology):
        tax_id = 'NCBITaxon:9606'

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        geno = Genotype(g)
        model = Model(g)
        # add the feature to the graph
        hap_description = None
        if risk_allele_frequency != '' and \
                risk_allele_frequency != 'NR':
            hap_description = \
                str(risk_allele_frequency) + \
                ' [risk allele frequency]'

        model.addIndividualToGraph(hap_id, hap_label.strip(),
                                   Feature.types['haplotype'], hap_description)
        geno.addTaxon(tax_id, hap_id)

        snp_labels = re.split(r';\s?', hap_label)
        chrom_nums = re.split(r';\s?', chrom_num)
        chrom_positions = re.split(r';\s?', chrom_pos)
        context_list = re.split(r';\s?', context)
        mapped_genes = re.split(r';\s?', mapped_gene)
        snp_curies = list()

        for index, snp in enumerate(snp_labels):
            snp_curie, snp_type = self._get_curie_and_type_from_id(snp)
            if snp_type is None:
                # make blank node
                snp_curie = self.make_id(snp, "_")

            g.addTriple(hap_id, geno.object_properties['has_variant_part'],
                        snp_curie)
            snp_curies.append(snp_curie)

        # courtesy http://stackoverflow.com/a/16720915
        length = len(snp_labels)
        if not all(len(lst) == length
                   for lst in [chrom_nums, chrom_positions, context_list]):
            logger.warn(
                "Unexpected data field for haplotype {} \n "
                "will not add snp details".format(hap_label))
            return

        variant_in_gene_count = 0
        for index, snp_curie in enumerate(snp_curies):
            self._add_snp_to_graph(
                snp_curie, snp_labels[index], chrom_nums[index],
                chrom_positions[index], context_list[index])

            if len(mapped_genes) == len(snp_labels):

                so_class = self._map_variant_type(context_list[index])

                if so_class is None:
                    raise ValueError("Unknown SO class {} in haplotype {}"
                                     .format(context_list[index], hap_label))
                so_query = """
                    SELECT ?variant_label
                    WHERE {{
                        {0} rdfs:subClassOf+ SO:0001564 ;
                            rdfs:label ?variant_label .
                    }}
                """.format(so_class)

                query_result = so_ontology.query(so_query)
                if len(list(query_result)) > 0:
                    gene_id = DipperUtil.get_ncbi_id_from_symbol(
                        mapped_genes[index])
                    if gene_id is not None:
                        geno.addAffectedLocus(snp_curie, gene_id)
                        geno.addAffectedLocus(hap_id, gene_id)
                        variant_in_gene_count += 1

                if context_list[index] == 'upstream_gene_variant':
                    gene_id = DipperUtil.get_ncbi_id_from_symbol(
                        mapped_genes[index])
                    if gene_id is not None:
                        g.addTriple(
                            snp_curie,
                            Feature.object_properties[
                                'upstream_of_sequence_of'],
                            gene_id)
                elif context_list[index] == 'downstream_gene_variant':
                    gene_id = DipperUtil.get_ncbi_id_from_symbol(
                        mapped_genes[index])
                    if gene_id is not None:
                        g.addTriple(
                            snp_curie,
                            Feature.object_properties[
                                'downstream_of_sequence_of'],
                            gene_id)
            else:
                logger.warn("More mapped genes than snps, "
                            "cannot disambiguate for {}".format(hap_label))

        # Seperate in case we want to apply a different relation
        # If not this is redundant with triples added above
        if len(mapped_genes) == variant_in_gene_count \
                and len(set(mapped_genes)) == 1:
            gene_id = DipperUtil.get_ncbi_id_from_symbol(mapped_genes[0])
            geno.addAffectedLocus(hap_id, gene_id)

        return

    def _add_snp_to_graph(
            self, snp_id, snp_label, chrom_num, chrom_pos, context,
            risk_allele_frequency=None):
        # constants
        tax_id = 'NCBITaxon:9606'
        genome_version = 'GRCh38'

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)

        if chrom_num != '' and chrom_pos != '':
            location = self._make_location_curie(chrom_num, chrom_pos)
            if location not in self.id_location_map:
                self.id_location_map[location] = set()
        else:
            location = None

        alteration = re.search(r'-(.*)$', snp_id)
        if alteration is not None \
                and re.match(r'[ATGC]', alteration.group(1)):
            # add variation to snp
            pass  # TODO

        if location is not None:
            self.id_location_map[location].add(snp_id)

        # create the chromosome
        chrom_id = makeChromID(chrom_num, genome_version, 'CHR')

        # add the feature to the graph
        snp_description = None
        if risk_allele_frequency is not None\
                and risk_allele_frequency != ''\
                and risk_allele_frequency != 'NR':
            snp_description = \
                str(risk_allele_frequency) + \
                ' [risk allele frequency]'

        f = Feature(
            g, snp_id, snp_label.strip(),
            Feature.types['SNP'], snp_description)
        if chrom_num != '' and chrom_pos != '':
            f.addFeatureStartLocation(chrom_pos, chrom_id)
            f.addFeatureEndLocation(chrom_pos, chrom_id)
        f.addFeatureToGraph()
        f.addTaxonToFeature(tax_id)
        # TODO consider adding allele frequency as property;
        # but would need background info to do that

        # also want to add other descriptive info about
        # the variant from the context
        for c in re.split(r';', context):
            cid = self._map_variant_type(c.strip())
            if cid is not None:
                model.addType(snp_id, cid)

        return

    def _add_deprecated_snp(self, snp_id, snp_id_current, merged,
                            chrom_num, chrom_pos):
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        location = self._make_location_curie(chrom_num, chrom_pos)
        # add deprecation information
        if merged == '1' and str(snp_id_current.strip()) != '':
            # get the current rs_id
            current_rs_id = 'dbSNP:'
            if not re.match(r'rs', snp_id_current):
                current_rs_id += 'rs'
            current_rs_id += str(snp_id_current)
            if location is not None:
                if location not in self.id_location_map:
                    self.id_location_map[location] = set(current_rs_id)
                else:
                    self.id_location_map[location].add(current_rs_id)
            model.addDeprecatedIndividual(snp_id, current_rs_id)
            # TODO check on this
            # should we add the annotations to the current
            # or orig?
            model.makeLeader(current_rs_id)
        else:
            model.makeLeader(snp_id)

    def _add_snp_gene_relation(self, snp_id, snp_gene_nums,
                               upstream_gene_num, downstream_gene_num):
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        geno = Genotype(g)
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
                    gene_id = 'NCBIGene:' + s
                    geno.addAffectedLocus(snp_id, gene_id)

        # add the up and downstream genes if they are available
        if upstream_gene_num != '':
            downstream_gene_id = 'NCBIGene:' + downstream_gene_num
            g.addTriple(
                snp_id,
                Feature.object_properties[
                    r'upstream_of_sequence_of'],
                downstream_gene_id)
        if downstream_gene_num != '':
            upstream_gene_id = 'NCBIGene:' + upstream_gene_num
            g.addTriple(
                snp_id,
                Feature.object_properties[
                    'downstream_of_sequence_of'],
                upstream_gene_id)

    def _add_variant_trait_association(self, variant_id, mapped_trait_uri,
                                       efo_ontology, pubmed_id,
                                       description=None):
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        # make associations to the EFO terms; there can be >1
        if mapped_trait_uri.strip() != '':
            for trait in re.split(r',', mapped_trait_uri):
                trait = trait.strip()

                cu = CurieUtil(curie_map.get())
                trait_id = cu.get_curie(trait)

                dis_query = """
                    SELECT ?trait
                    WHERE {{
                        {0} rdfs:subClassOf+ EFO:0000408 .
                        {0} rdfs:label ?trait .
                    }}
                """.format(trait_id)

                query_result = efo_ontology.query(dis_query)
                if len(list(query_result)) > 0:
                    if re.match(r'^EFO', trait_id):
                        model.addClassToGraph(trait_id, list(
                            query_result)[0][0], 'DOID:4')

                phenotype_query = """
                    SELECT ?trait
                    WHERE {{
                        {0} rdfs:subClassOf+ EFO:0000651 .
                        {0} rdfs:label ?trait .
                    }}
                """.format(trait_id)

                query_result = efo_ontology.query(phenotype_query)
                if len(list(query_result)) > 0:
                    if re.match(r'^EFO', trait_id):
                        model.addClassToGraph(
                            trait_id,
                            list(query_result)[0][0],
                            'UPHENO:0001001')

                pubmed_curie = 'PMID:' + pubmed_id

                ref = Reference(
                    g, pubmed_curie, Reference.ref_types['journal_article'])
                ref.addRefToGraph()

                assoc = G2PAssoc(
                    g, self.name, variant_id, trait_id,
                    model.object_properties['contributes_to'])
                assoc.add_source(pubmed_curie)
                # combinatorial evidence
                # used in automatic assertion
                eco_id = 'ECO:0000213'
                assoc.add_evidence(eco_id)

                if description is not None:
                    assoc.set_description(description)

                # FIXME score should get added to provenance/study
                # assoc.set_score(pvalue)
                assoc.add_association_to_graph()

    @staticmethod
    def _map_variant_type(sample_type):
        ctype = None
        type_map = {
            'stop_gained': 'SO:0001587',              # stop-gain variant
            'intron_variant': 'SO:0001627',           # intron variant
            '3_prime_UTR_variant': 'SO:0001624',      # 3'utr variant
            '5_prime_UTR_variant': 'SO:0001623',      # 5'UTR variant
            'synonymous_variant': 'SO:0001819',       # synonymous variant
            'frameshift_variant': 'SO:0001589',       # frameshift
            'intergenic_variant': 'SO:0001628',       # intergenic_variant
            'non_coding_transcript_exon_variant': 'SO:0001619', # noncoding transcript variant
            'splice_acceptor_variant': 'SO:0001574',  # splice acceptor variant
            'splice_donor_variant': 'SO:0001575',     # splice donor variant
            'missense_variant': 'SO:0001583',         # missense variant
            'downstream_gene_variant': 'SO:0001634',  # 500B_downstream_variant
            'upstream_gene_variant': 'SO:0001636',    # 2KB_upstream_variant
            'coding_sequence_variant': 'SO:0001580',  # coding_sequence_variant
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

    @staticmethod
    def _make_location_curie(chrom_num, chrom_pos):
        return 'chr' + str(chrom_num) + ':' + str(chrom_pos)

    @staticmethod
    def _make_description(disease_or_trait, initial_sample_description,
                          replicate_sample_description,
                          platform_with_snps_passing_qc, pvalue):
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
        description = ' '.join((description, '(p=' + pvalue + ')'))
        return description

    @staticmethod
    def _get_curie_and_type_from_id(variant_id):
        """
        Given a variant id, our best guess at its curie
        and type (snp, haplotype, etc)
        None will be used for both curie and type
        for IDs that we can't process
        :param variant_id:
        :return:
        """
        curie = None
        variant_type = None

        # remove space before hyphens
        variant_id = re.sub(r' -', '-', variant_id)
        if re.search(r' x ', variant_id) \
                or re.search(r',', variant_id):
            # TODO deal with rs1234 x rs234... (haplotypes?)
            logger.warning(
                "Cannot parse variant groups of this format: %s",
                variant_id)
        elif re.search(r';', variant_id):
            curie = ':haplotype_' + Source.hash_id(variant_id)
            variant_type = "haplotype"
        elif re.match(r'rs', variant_id):
            curie = 'dbSNP:' + variant_id.strip()
            curie = re.sub(r'-.*$', '', curie).strip()
            variant_type = "snp"
            # remove the alteration
        elif re.match(r'kgp', variant_id):
            # http://www.1000genomes.org/faq/what-are-kgp-identifiers
            curie = ':kgp-' + variant_id.strip()
            variant_type = "snp"
        elif re.match(r'chr', variant_id):
            # like: chr10:106180121-G
            #
            variant_id = re.sub(r'-?', '-N', variant_id)
            variant_id = re.sub(r' ', '', variant_id)
            curie = ':gwas-' + re.sub(
                r':', '-', variant_id.strip())
            variant_type = "snp"
        elif variant_id.strip() == '':
            pass
        else:
            logger.warning(
                "There's a snp id i can't manage: %s",
                variant_id)

        return curie, variant_type

    def getTestSuite(self):
        import unittest
        from tests.test_gwascatalog import GWASCatalogTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(GWASCatalogTestCase)

        return test_suite
