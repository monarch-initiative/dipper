import logging
import csv
import re
import json

from dipper.sources.Source import Source
from dipper.utils.DipperUtil import DipperUtil
from dipper.models.Model import Model
from dipper.models.Genotype import Genotype
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Reference import Reference
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv
from dipper.graph.RDFGraph import RDFGraph

logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)


class GWASCatalog(Source):
    """
    The NHGRI-EBI Catalog of published genome-wide association studies.

    We link the variants recorded here to the curated EFO-classes using a
    "contributes to" linkage because the only thing we know is that the SNPs
    are associated with the trait/disease,
    but we don't know if it is actually causative.

    Description of the GWAS catalog is here:
    http://www.ebi.ac.uk/gwas/docs/fileheaders#_file_headers_for_catalog_version_1_0_1

    GWAS also pulishes Owl files described here
    http://www.ebi.ac.uk/gwas/docs/ontology


    Status:  IN PROGRESS

    """

    GWASFTP = 'ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/'
    GWASFILE = 'gwas-catalog-associations_ontology-annotated.tsv'
    files = {
        'catalog': {
            'file': GWASFILE,
            'url': GWASFTP + GWASFILE,
            'columns': [  # expected
                # head -1 gwas-catalog-associations_ontology-annotated.tsv |
                # tr '\t' '\n'  | sed "s|\(.*\)|'\1',|g"
                'DATE ADDED TO CATALOG',
                'PUBMEDID',
                'FIRST AUTHOR',
                'DATE',
                'JOURNAL',
                'LINK',
                'STUDY',
                'DISEASE/TRAIT',
                'INITIAL SAMPLE SIZE',
                'REPLICATION SAMPLE SIZE',
                'REGION',
                'CHR_ID',
                'CHR_POS',
                'REPORTED GENE(S)',
                'MAPPED_GENE',
                'UPSTREAM_GENE_ID',
                'DOWNSTREAM_GENE_ID',
                'SNP_GENE_IDS',
                'UPSTREAM_GENE_DISTANCE',
                'DOWNSTREAM_GENE_DISTANCE',
                'STRONGEST SNP-RISK ALLELE',
                'SNPS',
                'MERGED',
                'SNP_ID_CURRENT',
                'CONTEXT',
                'INTERGENIC',
                'RISK ALLELE FREQUENCY',
                'P-VALUE',
                'PVALUE_MLOG',
                'P-VALUE (TEXT)',
                'OR or BETA',
                '95% CI (TEXT)',
                'PLATFORM [SNPS PASSING QC]',
                'CNV',
                'MAPPED_TRAIT',
                'MAPPED_TRAIT_URI',
                'STUDY ACCESSION',
                'GENOTYPING TECHNOLOGY',
            ],
        },
        'so': {
            'file': 'so.owl',
            'url': 'http://purl.obolibrary.org/obo/so.owl'},
        'mondo': {
            'file': 'mondo.json',
            'url': 'https://github.com/monarch-initiative/mondo/releases/'
                   'download/2019-04-06/mondo-minimal.json'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'gwascatalog',
            ingest_title='NHGRI-EBI Catalog of ' +
            'published genome-wide association studies',
            ingest_url='http://www.ebi.ac.uk/gwas/',
            license_url=None,
            data_rights='http://www.ebi.ac.uk/gwas/docs/about'
            # file_handle=None
        )

        if graph_type != 'rdf_graph':
            raise ValueError("GWAS Catalog requires a rdf_graph")

        if 'gene' not in self.all_test_ids:
            LOG.warning("not configured with test ids.")
        else:
            self.test_ids = self.all_test_ids

        # build a dictionary of genomic location to identifiers,
        # to try to get the equivalences
        self.id_location_map = dict()

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced:
        :return:
        """
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self.process_catalog(limit)

        LOG.info("Finished parsing.")

    def process_catalog(self, limit=None):
        """
        :param limit:
        :return:

        """
        src_key = 'catalog'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing Data from %s", raw)

        so_ontology = RDFGraph(False, "SO")
        LOG.info("Loading SO ontology in separate rdf graph")
        so_ontology.parse(self.files['so']['url'], format='xml')
        so_ontology.bind_all_namespaces()
        LOG.info("Finished loading SO ontology")

        mondo_file = '/'.join((self.rawdir, self.files['mondo']['file']))
        with open(mondo_file, 'r') as mondo_fh:
            mondo_data = json.load(mondo_fh)

        col = self.files[src_key]['columns']

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            row = next(reader)
            if self.check_fileheader(col, row):
                pass
            for row in reader:
                if len(col) != len(row):
                    LOG.error('BadRow: %i has %i columns', reader.line_num, row)
                    continue
                # head -1 gwas-catalog-associations_ontology-annotated.tsv |
                # tr '\t' '\n'  | sed "s|\(.*\)|# = row[col.index('\1')]|g"

                # = row[col.index('DATE ADDED TO CATALOG')]
                pubmed_num = row[col.index('PUBMEDID')].strip()
                # = row[col.index('FIRST AUTHOR')]
                # = row[col.index('DATE')]
                # = row[col.index('JOURNAL')]
                # = row[col.index('LINK')]
                # = row[col.index('STUDY')]
                disease_or_trait = row[col.index('DISEASE/TRAIT')].strip()
                initial_sample_description = row[
                    col.index('INITIAL SAMPLE SIZE')].strip()
                replicate_sample_description = row[
                    col.index('REPLICATION SAMPLE SIZE')].strip()
                # = row[col.index('REGION')]
                chrom_num = row[col.index('CHR_ID')].strip()
                chrom_pos = row[col.index('CHR_POS')].strip()
                # = row[col.index('REPORTED GENE(S)')]
                mapped_gene = row[col.index('MAPPED_GENE')].strip()
                upstream_gene_num = row[col.index('UPSTREAM_GENE_ID')].strip()
                downstream_gene_num = row[col.index('DOWNSTREAM_GENE_ID')].strip()
                snp_gene_nums = row[col.index('SNP_GENE_IDS')].strip()
                # = row[col.index('UPSTREAM_GENE_DISTANCE')]
                # = row[col.index('DOWNSTREAM_GENE_DISTANCE')]
                strongest_snp_risk_allele = row[
                    col.index('STRONGEST SNP-RISK ALLELE')].strip()
                # = row[col.index('SNPS')]
                merged = row[col.index('MERGED')].strip()
                snp_id_current = row[col.index('SNP_ID_CURRENT')].strip()
                context = row[col.index('CONTEXT')].strip()
                # = row[col.index('INTERGENIC')]
                risk_allele_frequency = row[col.index('RISK ALLELE FREQUENCY')].strip()
                pvalue = row[col.index('P-VALUE')].strip()
                # = row[col.index('PVALUE_MLOG')]
                # = row[col.index('P-VALUE (TEXT)')]
                # = row[col.index('OR or BETA')]
                # = row[col.index('95% CI (TEXT)')]
                platform_with_snps_passing_qc = row[
                    col.index('PLATFORM [SNPS PASSING QC]')].strip()
                # = row[col.index('CNV')]
                mapped_trait = row[col.index('MAPPED_TRAIT')].strip()
                mapped_trait_uri = row[col.index('MAPPED_TRAIT_URI')].strip()
                # = row[col.index('STUDY ACCESSION')]
                # = row[col.index('GENOTYPING TECHNOLOGY')]

                if self.test_mode:
                    continue

# 06-May-2015	25917933
#   Zai CC	20-Nov-2014	J Psychiatr Res	http://europepmc.org/abstract/MED/25917933
# A genome-wide association study of suicide severity scores in bipolar disorder.
# Suicide in bipolar disorder
# 959 European ancestry individuals	NA
# 10p11.22	10	32704340	C10orf68, CCDC7, ITGB1	CCDC7
# rs7079041-A	rs7079041	0	7079041	intron	0		2E-6	5.698970

                snp_id_current = snp_id_current.split(' ')[0]

                # note: that these will no longer pattern match other instances

                variant_curie, variant_type = self._get_curie_and_type_from_id(
                    strongest_snp_risk_allele)

                if strongest_snp_risk_allele == '':
                    LOG.debug(
                        "No strongest SNP risk allele for %s:\n%s",
                        pubmed_num, str(row))
                    # still consider adding in the EFO terms
                    # for what the study measured?
                    continue
                if variant_curie is not None and variant_curie[0] == '_' and \
                        strongest_snp_risk_allele is not None:
                    self.graph.addTriple(
                        variant_curie, self.globaltt['label'],
                        strongest_snp_risk_allele, object_is_literal=True,
                        subject_category=blv.terms.SequenceVariant)

                if variant_type == 'snp':
                    self._add_snp_to_graph(
                        variant_curie, strongest_snp_risk_allele, chrom_num,
                        chrom_pos, context, risk_allele_frequency)
                    self._add_deprecated_snp(
                        variant_curie, snp_id_current, merged, chrom_num, chrom_pos)

                    self._add_snp_gene_relation(
                        variant_curie, snp_gene_nums, upstream_gene_num,
                        downstream_gene_num)
                elif variant_type == 'haplotype':
                    self._process_haplotype(
                        variant_curie, strongest_snp_risk_allele, chrom_num,
                        chrom_pos, context, risk_allele_frequency, mapped_gene,
                        so_ontology)
                elif variant_type is None and snp_id_current != '':
                    LOG.warning(
                        "There's a snp id we can't manage: %s",
                        strongest_snp_risk_allele)
                    continue

                description = self._make_description(
                    disease_or_trait, initial_sample_description,
                    replicate_sample_description,
                    platform_with_snps_passing_qc, pvalue)

                self._add_variant_trait_association(
                    variant_curie, mapped_trait_uri, mapped_trait, mondo_data,
                    pubmed_num, description)

                if not self.test_mode and (
                        limit is not None and reader.line_num > limit):
                    break

        # TODO loop through the location hash,
        # and make all snps at that location equivalent
        for loc in self.id_location_map:
            snp_ids = self.id_location_map[loc]
            if len(snp_ids) > 1:
                LOG.info("%s has >1 snp id: %s", loc, str(snp_ids))

    def _process_haplotype(
            self, hap_id, hap_label, chrom_num, chrom_pos, context,
            risk_allele_frequency, mapped_gene, so_ontology):

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        model = Model(graph)
        # add the feature to the graph
        hap_description = None
        if risk_allele_frequency not in ['', 'NR']:
            hap_description = str(risk_allele_frequency) + ' [risk allele frequency]'

        model.addIndividualToGraph(
            hap_id, hap_label.strip(), self.globaltt['haplotype'], hap_description)
        geno.addTaxon(self.globaltt["Homo sapiens"], hap_id)

        snp_labels = re.split(r';\s?', hap_label)
        chrom_nums = re.split(r';\s?', chrom_num)
        chrom_positions = re.split(r';\s?', chrom_pos)
        context_list = re.split(r';\s?', context)
        mapped_genes = re.split(r';\s?', mapped_gene)

        # Not having four "PAX5" as a list might be better, but it breaks unit tests
        # mapped_genes = list(set(mapped_genes)) # make uniq
        # snp_labels = list(set(snp_labels)) # make uniq

        snp_curies = list()

        for snp in snp_labels:
            snp_curie, snp_type = self._get_curie_and_type_from_id(snp)
            if snp_type is None:
                LOG.info('cant find type for SNP in %s', snp)
                # make blank node
                snp_curie = self.make_id(snp, "_")
                model.addLabel(snp_curie, snp)
            elif snp_curie[0] == '_':   # arrived an unlabeled blanknode
                model.addLabel(snp_curie, snp)

            graph.addTriple(hap_id, self.globaltt['has_variant_part'], snp_curie)
            snp_curies.append(snp_curie)

        # courtesy http://stackoverflow.com/a/16720915
        # check lengths of mutiple lists
        length = len(snp_curies)
        if not all(len(lst) == length
                   for lst in [snp_labels, chrom_nums, chrom_positions, context_list]):
            LOG.warning(
                "Incongruous data field(s) for haplotype %s \n "
                "will not add snp details", hap_label)
        else:

            variant_in_gene_count = 0
            for index, snp_curie in enumerate(snp_curies):
                self._add_snp_to_graph(
                    snp_curie, snp_labels[index], chrom_nums[index],
                    chrom_positions[index], context_list[index])

                if mapped_genes and len(mapped_genes) != len(snp_labels):
                    LOG.warning(
                        "More mapped genes than snps,"
                        " cannot disambiguate for\n%s\n%s",
                        mapped_genes, snp_labels)  # hap_label)
                else:
                    so_class = self.resolve(context_list[index])
                    so_query = """
        SELECT ?variant_label
        WHERE {{
            {0} rdfs:subClassOf+ {1} ;
            rdfs:label ?variant_label .
        }}
                    """.format(so_class, self.globaltt['gene_variant'])

                    query_result = so_ontology.query(so_query)

                    gene_id = DipperUtil.get_hgnc_id_from_symbol(mapped_genes[index])

                    if gene_id is not None and len(list(query_result)) == 1:
                        if context_list[index] in ['upstream_gene_variant',
                                                   'downstream_gene_variant']:
                            graph.addTriple(
                                snp_curie, self.resolve(context_list[index]), gene_id)
                        else:
                            geno.addAffectedLocus(snp_curie, gene_id)
                            variant_in_gene_count += 1

            # Seperate in case we want to apply a different relation
            # If not this is redundant with triples added above
            if len(mapped_genes) == variant_in_gene_count and \
                    len(set(mapped_genes)) == 1:
                gene_id = DipperUtil.get_hgnc_id_from_symbol(mapped_genes[0])
                geno.addAffectedLocus(hap_id, gene_id)

    def _add_snp_to_graph(
            self, snp_id, snp_label, chrom_num, chrom_pos, context,
            risk_allele_frequency=None):

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)

        if chrom_num != '' and chrom_pos != '':
            location = self._make_location_curie(chrom_num, chrom_pos)
            if location not in self.id_location_map:
                self.id_location_map[location] = set()
        else:
            location = None

        alteration = re.search(r'-(.*)$', snp_id)
        if alteration is not None and re.match(r'[ATGC]', alteration.group(1)):
            # add variation to snp
            pass  # TODO

        if location is not None:
            self.id_location_map[location].add(snp_id)

        # create the chromosome
        chrom_id = makeChromID(chrom_num, self.localtt['reference assembly'], 'CHR')

        # add the feature to the graph
        snp_description = None
        if risk_allele_frequency is not None\
                and risk_allele_frequency != ''\
                and risk_allele_frequency != 'NR':
            snp_description = str(risk_allele_frequency) + ' [risk allele frequency]'

        feat = Feature(
            graph, snp_id, snp_label.strip(), self.globaltt['SNP'], snp_description,
            feature_category=blv.terms.SequenceVariant)
        if chrom_num != '' and chrom_pos != '':
            feat.addFeatureStartLocation(chrom_pos, chrom_id)
            feat.addFeatureEndLocation(chrom_pos, chrom_id)
        feat.addFeatureToGraph()
        feat.addTaxonToFeature(self.globaltt['Homo sapiens'])
        # TODO consider adding allele frequency as property;
        # but would need background info to do that

        # also want to add other descriptive info about
        # the variant from the context
        for ctx in context.split(';'):
            ctx = ctx.strip()
            if ctx:
                cid = self.resolve(ctx, False)
                if cid != ctx:
                    model.addType(snp_id, cid,
                                  subject_category=blv.terms.SequenceVariant)

    def _add_deprecated_snp(
            self, snp_id, snp_id_current, merged, chrom_num, chrom_pos):
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        location = self._make_location_curie(chrom_num, chrom_pos)
        # add deprecation information
        if merged == '1' and snp_id_current != '':
            current_rs_id = 'dbSNP:rs' + snp_id_current
            if location is not None:
                if location not in self.id_location_map:
                    self.id_location_map[location] = set(current_rs_id)
                else:
                    self.id_location_map[location].add(current_rs_id)
            model.addDeprecatedIndividual(snp_id, current_rs_id,
                                          blv.terms.SequenceVariant,
                                          blv.terms.SequenceVariant)
            # TODO check on this
            # should we add the annotations to the current
            # or orig?
            model.makeLeader(current_rs_id,
                             node_category=blv.terms.SequenceVariant)
        else:
            model.makeLeader(snp_id,
                             node_category=blv.terms.SequenceVariant)

    def _add_snp_gene_relation(
            self, snp_id, snp_gene_nums, upstream_gene_num, downstream_gene_num):
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        # add the feature as a sequence alteration
        # affecting various genes
        # note that intronic variations don't necessarily list
        # the genes such as for rs10448080  FIXME
        if snp_gene_nums != '':
            for geneid in re.split(r',', snp_gene_nums):
                geneid = geneid.strip()
                # still have to test for this,
                # because sometimes there's a leading comma
                if geneid != '':
                    geno.addAffectedLocus(snp_id, 'ENSEMBL:' + geneid)

        # add the up and downstream genes if they are available
        if upstream_gene_num != '':
            downstream_gene_id = 'ENSEMBL:' + downstream_gene_num
            graph.addTriple(
                snp_id, self.globaltt['is upstream of sequence of'], downstream_gene_id,
                subject_category=blv.terms.SequenceVariant,
                object_category=blv.terms.Gene)
        if downstream_gene_num != '':
            upstream_gene_id = 'ENSEMBL:' + upstream_gene_num
            graph.addTriple(
                snp_id, self.globaltt['is downstream of sequence of'], upstream_gene_id,
                subject_category=blv.terms.SequenceVariant,
                object_category=blv.terms.Gene)

    def _add_variant_trait_association(
            self, variant_id, mapped_trait_uri, mapped_trait, mondo_data, pubmed_id,
            description=None):
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)

        # Use mapped traits for labels, hope that labels do not contain commas

        mapped_traits = [trait.strip() for trait in mapped_trait.split(',')]
        mapped_trait_uris = [iri.strip() for iri in mapped_trait_uri.split(',')]
        if mapped_trait_uris:
            for index, trait in enumerate(mapped_trait_uris):
                trait_curie = trait.replace("http://www.ebi.ac.uk/efo/EFO_", "EFO:")

                if not DipperUtil.is_id_in_mondo(trait_curie, mondo_data):
                    if re.match(r'^EFO', trait_curie):
                        model.addClassToGraph(
                            trait_curie, mapped_traits[index],
                            self.globaltt['phenotype'],
                            class_category=blv.terms.PhenotypicFeature)
                    LOG.debug("{} not in mondo".format(trait_curie))
                else:
                    LOG.debug("{} in mondo".format(trait_curie))

                pubmed_curie = 'PMID:' + pubmed_id

                ref = Reference(
                    graph, pubmed_curie, self.globaltt['journal article'])
                ref.addRefToGraph()

                assoc = G2PAssoc(
                    graph, self.name, variant_id, trait_curie,
                    model.globaltt['contributes to condition'])
                assoc.add_source(pubmed_curie)

                assoc.add_evidence(
                    self.globaltt['combinatorial evidence used in automatic assertion'])

                if description is not None:
                    assoc.set_description(description)

                # FIXME score should get added to provenance/study
                # assoc.set_score(pvalue)
                if trait_curie is not None and variant_id is not None:
                    assoc.add_association_to_graph()

    @staticmethod
    def _make_location_curie(chrom_num, chrom_pos):
        return 'chr' + str(chrom_num) + ':' + str(chrom_pos)

    @staticmethod
    def _make_description(
            disease_or_trait, initial_sample_description, replicate_sample_description,
            platform_with_snps_passing_qc, pvalue):
        description = 'A study of '+disease_or_trait+' in '+initial_sample_description
        if replicate_sample_description != '':
            description = ' '.join(
                (description, 'with', replicate_sample_description))
        if platform_with_snps_passing_qc != '':
            description = ' '.join(
                (description, 'on platform', platform_with_snps_passing_qc))
        description = ' '.join((description, '(p=' + pvalue + ')'))
        return description

    @staticmethod
    def _get_curie_and_type_from_id(variant_id):
        """
        Given a variant id, our best guess at its curie and type (snp, haplotype, etc)
        'None' will be used for both curie and type  for IDs that we can't process

        # 2019-May three snp-id have  ' e' or ' a'  appended. note space.
        # examples: 'rs2440154 e-A'  and 'rs2440154 e'
        # including the suffix in the url is a web noop but breaks rdflib

        :param variant_id:
        :return:
        """
        curie = None
        variant_type = None

        # remove space before hyphens
        variant_id = re.sub(r' -', '-', variant_id).strip()
        if re.search(r' x ', variant_id) or re.search(r',', variant_id):
            # TODO deal with rs1234 x rs234... (haplotypes?)
            LOG.warning("Cannot parse variant groups of this format: %s", variant_id)
        elif re.search(r';', variant_id):
            curie = ':haplotype_' + Source.hash_id(variant_id)   # deliberate 404
            variant_type = "haplotype"
        elif variant_id[:2] == 'rs':
            # remove whitespace from errant id, rs6194 5053-?
            curie = 'dbSNP:' + variant_id.split('-')[0].replace(' ', '')
            # curie = re.sub(r'-.*$', '', curie).strip()
            variant_type = "snp"
            # remove the alteration
        elif variant_id[:3] == 'kgp':
            # http://www.1000genomes.org/faq/what-are-kgp-identifiers
            curie = 'GWAS:' + variant_id.split('-')[0]
            variant_type = "snp"
        elif variant_id[:3] == 'chr':
            # like: chr10:106180121-G
            variant_id = re.sub(r'-?', '-N', variant_id)
            variant_id = re.sub(r' ', '', variant_id)
            # going to hate myself but ...
            # moving this from a broken base node to yet another blank node
            # It had produced this monstrocity with the embedded quote
            # :gwas--Nc-Nh-Nr-N1-N1-N--N1-N0-N2-N7-N5-N1-N1-N0-N2-N"-N?-N
            curie = Source.make_id('gwas-' + re.sub(r':', '-', variant_id), '_')
            variant_type = "snp"
        elif variant_id.strip() == '':
            pass
        else:
            LOG.warning("There's a snp id i can't manage: %s", variant_id)

        return curie, variant_type
