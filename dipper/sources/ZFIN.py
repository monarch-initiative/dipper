import csv
import re
import logging
import os
import yaml

from intermine.webservice import Service
from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Genotype import Genotype
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Environment import Environment
from dipper.models.GenomicFeature import makeChromID
from dipper.models.GenomicFeature import Feature
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)
ZFDL = 'http://zfin.org/downloads'
ZPCURATION = 'http://purl.obolibrary.org/obo/zp/'


class ZFIN(Source):
    """
    This is the parser for the
    [Zebrafish Model Organism Database (ZFIN)](http://www.zfin.org),
    from which we process genotype and phenotype data for laboratory zebrafish.

    We generate the zfin graph to include the following information:
    * genes
    * sequence alterations
    (includes SNPs/del/ins/indel and large chromosomal rearrangements)
    * transgenic constructs
    * morpholinos, talens, crisprs as expression-affecting reagents
    * genotypes, and their components
    * fish (as comprised of intrinsic and extrinsic genotypes)
    * publications (and their mapping to PMIDs, if available)
    * genotype-to-phenotype associations
    (including environments and stages at which they are assayed)
    * environmental components
    * orthology to human genes
    * genetic positional information for genes and sequence alterations
    * fish-to-disease model associations

    Genotypes leverage the GENO genotype model and include
    both intrinsic and extrinsic genotypes.
    Where necessary, we create anonymous nodes of the genotype partonomy
    (such as for variant single locus complements,
    genomic variation complements, variant loci, extrinsic genotypes,
    and extrinsic genotype parts).

    Genotype labels are output as ZFIN genotype name + "[background]".
    We also process the genotype components to build labels
    in a monarch-style, and these are added as synonyms. The monarch-style
    genotype label includes:
    *  all genes targeted by reagents (morphants, crisprs, etc),
    in addition to the ones that the reagent was designed against.
    *  all affected genes within deficiencies
    *  complex hets being listed as gene<mutation1>/gene<mutation2>
    rather than gene<mutation1>/+; gene<mutation2>/+


    see: resources/zfin/README for column extraction from downloads page

    """

    files = {
        'geno': {
            'file': 'genotype_features.txt',
            'url': ZFDL + '/genotype_features.txt',
            'columns': [
                'Genotype ID',
                'Genotype Name',
                'Genotye Unique Name',
                'Allele ID',
                'Allele Name',
                'Allele Abbreviation',
                'Allele Type',
                'Allele Display Type',
                'Gene or Construct Symbol',
                'Corresponding ZFIN Gene ID/Construct ID',
                'Allele Zygosity',
                'Construct Name',
                'Construct ZdbId',
            ]},

        'pheno': {
            'file': 'phenotype_fish.txt',
            'url': ZFDL + '/phenotype_fish.txt',
            'columns': [
                'Fish ID',
                'Fish Name',
                'Start Stage ID',
                'Start Stage Name',
                'End Stage ID',
                'End Stage Name',
                'Affected Structure or Process 1 subterm ID',
                'Affected Structure or Process 1 subterm Name',
                'Post-composed Relationship ID',
                'Post-composed Relationship Name',
                'Affected Structure or Process 1 superterm ID',
                'Affected Structure or Process 1 superterm Name',
                'Phenotype Keyword ID',
                'Phenotype Keyword Name',
                'Phenotype Tag',
                'Affected Structure or Process 2 subterm ID',
                'Affected Structure or Process 2 subterm name',
                'Post-composed Relationship (rel) ID',
                'Post-composed Relationship (rel) Name',
                'Affected Structure or Process 2 superterm ID',
                'Affected Structure or Process 2 superterm name',
                'Publication ID',
                'Environment ID',
            ]},
        'pubs': {
            'file': 'zfinpubs.txt',
            'url': ZFDL + '/zfinpubs.txt',
            'columns': [
                'Publication ID',
                'pubMed ID (none or blank when not available)',
                'Authors',
                'Title',
                'Journal',
                'Year',
                'Volume',
                'Pages',
            ]},
        'zpmap': {  # this is a zp mapping file generated by Nico
            'file': 'id_map_zfin.tsv',
            'url': ZPCURATION + '/id_map_zfin.tsv',
            'columns': [
                'iri',
                'id',
            ]},
        'morph': {
            'file': 'Morpholinos.txt',
            'url': ZFDL + '/Morpholinos.txt',
            'columns': [
                'Gene ID',
                'Gene SO ID',
                'Gene Symbol',
                'Morpholino ID',
                'Morpholino SO ID',
                'Morpholino Symbol',
                'Morpholino Sequence',
                'Publication(s)',
                'Note',
            ]},
        # 'enviro': {'file': 'pheno_environment.txt',
        #            'url': ZFDL + '/pheno_environment.txt','columns': []},
        'enviro': {
            'file': 'pheno_environment_fish.txt',
            'url': ZFDL + '/pheno_environment_fish.txt',
            'columns': [
                'Environment ID',
                'ZECO Term Name',
                'ZECO Term ID (ZECO:ID)',
                'Chebi Term Name',
                'Chebi Term ID (Chebi:ID)',
                'ZFA Term Name',
                'ZFA Term ID (ZFA:ID)',
                'Affected Structure Subterm Name',
                'Affected Structure Subterm ID (GO-CC:ID)',
                'NCBI Taxon Name',
                'NCBI Taxon ID (NCBI Taxon:ID)',
            ]},
        'stage': {
            'file': 'stage_ontology.txt',
            'url': 'http://zfin.org/Downloads/stage_ontology.txt',
            'columns': [
                'Stage ID',
                'Stage OBO ID',
                'Stage Name',
                'Begin Hours',
                'End Hours',
            ]},
        # 'wild_expression': {
        #   'file': 'wildtype-expression.txt',
        #   'url': ZFDL + '/wildtype-expression.txt','columns': []},
        'mappings': {
            'file': 'mappings.txt',
            'url': ZFDL + '/mappings.txt',
            'columns': [
                'ZFIN ID',
                'Symbol',
                'SO_id',
                'Panel Symbol',
                'Chromosome',
                'Location',
                'Metric',
            ]},
        'backgrounds': {
            'file': 'genotype_backgrounds.txt',
            'url': ZFDL + '/genotype_backgrounds.txt',
            'columns': [
                'Genotype ID',
                'Genotype Name',
                'Background',
                'Background Name',
            ]},
        'genbank': {
            'file': 'genbank.txt',
            'url': ZFDL + '/genbank.txt',
            'columns': [
                'ZFIN ID',
                'SO ID',
                'Name',
                'GenBank ID',
            ]},
        'uniprot': {
            'file': 'uniprot.txt',
            'url': ZFDL + '/uniprot.txt',
            'columns': [
                'ZFIN ID',
                'SO ID',
                'Symbol',
                'UniProt ID',
            ]},
        'gene': {
            'file': 'gene.txt',
            'url': 'http://zfin.org/downloads/gene.txt',
            'columns': [
                'ZFIN ID',
                'SO ID',
                'Symbol',
                'NCBI Gene ID',
            ]},

        'wild': {
            'file': 'wildtypes.txt',
            'url': ZFDL + '/wildtypes_fish.txt',
            'columns': [
                'Fish ID',
                'Fish Name',
                'Fish Abbreviation',
                'Genotype ID',
            ]},
        'human_orthos': {
            'file': 'human_orthos.txt',
            'url': ZFDL + '/human_orthos.txt',
            'columns': [
                'ZFIN ID',
                'ZFIN Symbol',
                'ZFIN Name',
                'Human Symbol',
                'Human Name',
                'OMIM ID',
                'Gene ID',
                'HGNC ID',
                'Evidence',
                'Pub ID',
            ]},
        'features': {
            'file': 'features.txt',
            'url': ZFDL + '/features.txt',
            'columns': [
                'Genomic Feature ID',
                'Feature SO ID',
                'Genomic Feature Abbreviation',
                'Genomic Feature Name',
                'Genomic Feature Type',
                'Mutagen',
                'Mutagee',
                'Construct ID',
                'Construct name',
                'Construct SO ID',
                'TALEN/CRISPR ID',
                'TALEN/CRISPR Name',
            ]},
        'feature_affected_gene': {
            'file': 'features-affected-genes.txt',
            'url': ZFDL + '/features-affected-genes.txt',
            'columns': [
                'Genomic Feature ID',
                'Feature SO ID',
                'Genomic Feature Abbreviation',
                'Gene Symbol',
                'Gene ID',
                'Gene SO ID',
                'Genomic Feature - Marker Relationship',
                'Feature Type',
                'DNA/cDNA Change SO ID',
                'Reference Nucleotide',
                'Mutant Nucleotide',
                'Base Pairs Added',
                'Base Pairs Removed',
                'DNA/cDNA Change Position Start',
                'DNA/cDNA Change Position End',
                'DNA/cDNA Reference Sequence',
                'DNA/cDNA Change Localization',
                'DNA/cDNA Change Localization SO ID',
                'DNA/cDNA Change Localization Exon',
                'DNA/cDNA Change Localization Intron',
                'Transcript Consequence',
                'Transcript Consequence SO ID',
                'Transcript Consequence Exon',
                'Transcript Consequence Intron',
                'Protein Consequence',
                'Protein Consequence SO ID',
                'Reference Amino Acid',
                'Mutant Amino Acid',
                'Amino Acids Added',
                'Amino Acids Removed',
                'Protein Consequence Position Start',
                'Protein Consequence Position End',
                'Protein Reference Sequence',
            ]},
        'gene_marker_rel': {
            'file': 'gene_marker_relationship.txt',
            'url': ZFDL + '/gene_marker_relationship.txt',
            'columns': [
                'Gene ID',
                'Gene SO ID',
                'Gene Symbol',
                'Marker ID',
                'Marker SO ID',
                'Marker Symbol',
                'Relationship',
            ]},
        'crispr': {
            'file': 'CRISPR.txt',
            'url': ZFDL + '/CRISPR.txt',
            'columns': [
                'Gene ID',
                'Gene SO ID',
                'Gene Symbol',
                'CRISPR ID',
                'CRISPR SO ID',
                'CRISPR Symbol',
                'CRISPR Target Sequence',
                'Publication(s)',
                'Note',
            ]},
        'talen': {
            'file': 'TALEN.txt',
            'url': ZFDL + '/TALEN.txt',
            'columns': [
                'Gene ID',
                'Gene SO ID',
                'Gene Symbol',
                'TALEN ID',
                'TALEN SO ID',
                'TALEN Symbol',
                'TALEN Target Sequence 1',
                'TALEN Target Sequence 2',
                'Publication(s)',
                'Note',
            ]},
        'pub2pubmed': {
            'file': 'pub_to_pubmed_id_translation.txt',
            'url': ZFDL + '/pub_to_pubmed_id_translation.txt',
            'columns': [
                'Publication ZFIN ID',
                'PubMed ID (none or blank when not available)',
            ]},
        'gene_coordinates': {
            'file': 'E_zfin_gene_alias.gff3',
            'url': ZFDL + '/E_zfin_gene_alias.gff3',
            'columns': [  # basic gff3?
                'Chromosome',
                'Source',
                'Type',
                'Start',
                'End',
                'Score',
                'Strand',
                'Phase',
                'Attributes',
            ]},
        'fish_disease_models': {
            'file': 'fish_model_disease.txt',
            'url': ZFDL + '/fish_model_disease.txt',
            'columns': [
                'Fish ZDB ID',
                'Environment ZDB ID',
                'is_a_model',
                'DO Term ID',
                'DO Term Name',
                'Publication ZDB ID',
                'PubMed ID',
                'Evidence Code',
            ]},
        'fish_components': {
            'file': 'fish_components_fish.txt',
            'url': ZFDL + '/fish_components_fish.txt',
            'columns': [
                'Fish ID',
                'Fish Name',
                'Gene ID',
                'Gene Symbol',
                'Affector ID',
                'Affector Symbol',
                'Construct ID',
                'Construct Symbol',
                'Background ID',
                'Background Name',
                'Genotype ID',
                'Genotype Name',
            ]},
        'zmine_ortho_evidence': {
            'file': 'zmine_ortho_evidence.txt',
            # yep it is odd; see: get_orthology_sources_from_zebrafishmine()
            'url': 'http://0.0.0.0',
            'columns': [
                'zfin_gene_num',
                'zfin_gene_symbol',
                'ortholog_gene_symbol',
                'ortholog_ncbigene_num',
                'evidence_code',
                'zfin_pub_num',
                'pubmed_num',
            ]},
    }

    # load test_ids
    with open(
            os.path.join(os.path.dirname(__file__),
                         '../../tests/resources/zfin/zfin_test_ids.yaml')) as fhandle:
        test_ids = yaml.safe_load(fhandle)

    def __init__(self, graph_type, are_bnodes_skolemized, data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='zfin',
            ingest_title='Zebra Fish Information Network',
            ingest_url='https://zfin.org',
            ingest_logo="source-zfin.png",
            license_url=None,
            data_rights='http://zfin.org/warranty.html'
            # file_handle=None
        )

        self.dataset.set_citation(
            'https://wiki.zfin.org/display/general/ZFIN+db+information')

        self.fish_parts = {}
        self.geno_alleles = {}
        # to hold any label for a given id
        self.id_label_map = {}
        # to hold the mappings between genotype and background
        self.genotype_backgrounds = {}
        self.extrinsic_id_to_enviro_id_hash = {}
        # to hold the parts that are introduced from a construct
        self.transgenic_parts = {}
        # to hold the genes variant due to a seq alt
        self.variant_loci_genes = {}
        # to hold the parts of an environment
        self.environment_hash = {}
        self.wildtype_genotypes = []
        self.zp_map = {}

        if 'disease' not in self.all_test_ids:
            LOG.warning("not configured with gene test ids.")
            self.test_ids['disease'] = []
        else:
            self.test_ids['disease'] = self.all_test_ids['disease']
        self.default_taxon_id = self.globaltt['Danio rerio']
        self.mapped_zpids = list()

    def fetch(self, is_dl_forced=False):
        # fetch all the files; versions are set by the date of download.
        self.get_files(is_dl_forced)
        self.get_orthology_sources_from_zebrafishmine()

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)
        LOG.info("Parsing files...")

        self.zp_map = self._load_zp_mappings('zpmap')
        if self.test_only:
            self.test_mode = True

        # basic information on classes and instances
        self._process_genes(limit)
        self._process_stages(limit)
        self._process_pubinfo(limit)
        self._process_pub2pubmed(limit)

        # The knockdown reagents
        for kdr in ['morph', 'crispr', 'talen']:
            self._process_targeting_reagents(kdr, limit)

        self._process_gene_marker_relationships(limit)
        self._process_features(limit)
        self._process_feature_affected_genes(limit)
        # only adds features on chromosomes, not positions
        self._process_mappings(limit)

        # These must be processed before G2P and expression
        self._process_wildtypes(limit)
        self._process_genotype_backgrounds(limit)
        # REVIEWED - NEED TO REVIEW LABELS ON Deficiencies
        self._process_genotype_features(limit)

        self.process_fish(limit)
        # Must be processed after morpholinos/talens/crisprs id/label
        # self._process_pheno_enviro(limit)  # TODO waiting on issue #385

        # once the genotypes and environments are processed,
        # we can associate these with the phenotypes
        self._process_g2p(limit)
        self.process_fish_disease_models(limit)

        # zfin-curated orthology calls to human genes
        self._process_human_orthos(limit)
        self.process_orthology_evidence(limit)

        # coordinates of all genes - from ensembl
        self._process_gene_coordinates(limit)

        # FOR THE FUTURE - needs verification
        # self._process_wildtype_expression(limit)
        # self._process_uniprot_ids(limit)

        LOG.info("Finished parsing all files.")

    def process_fish(self, limit=None):
        """
        Fish give identifiers to the "effective genotypes" that we create.
        We can match these by:
        Fish = (intrinsic) genotype + set of morpholinos

        We assume here that the intrinsic genotypes and their parts
        will be processed separately, prior to calling this function.

        :param limit:
        :return:

        """
        src_key = 'fish_components'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing Fish Parts from file: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        taxon_id = self.default_taxon_id

        allele_to_construct_hash = {}
        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                fish_num = row[col.index('Fish ID')]
                fish_name = row[col.index('Fish Name')]
                # gene_num = row[col.index('Gene ID')]
                # gene_symbol = row[col.index('Gene Symbol')]
                affector_num = row[col.index('Affector ID')]
                # affector_symbol = row[col.index('Affector Symbol')]
                construct_num = row[col.index('Construct ID')]
                # construct_symbol = row[col.index('Construct Symbol')]
                # background_num = row[col.index('Background ID')]
                # background_symbol = row[col.index('Background Name')]
                genotype_num = row[col.index('Genotype ID')]
                # genotype_name = row[col.index('Genotype Name')]

                # fish have the following components:
                #  *  genotype, which is the intrinsic genotype;
                #     this may be a genetic background (WT)
                #  *  an optional background for the intrinsic genotype
                #  *  affectors == alleles or morphants
                #  *  constructs which may give rise to the affectors
                #  *  affected genes

                if fish_num not in self.fish_parts:
                    self.fish_parts[fish_num] = {}
                    self.fish_parts[fish_num] = {
                        'intrinsic_genotype': genotype_num,
                        'affectors': set(),
                        'fish_label': fish_name
                    }

                # HACK - bad allele id - replace it with the new one  FIXME
                if affector_num == 'ZDB-ALT-090504-1':
                    affector_num = 'ZDB-ALT-040723-4'

                self.fish_parts[fish_num]['affectors'].add(affector_num)

                # add the constructs that the allele comes from
                if construct_num != '':
                    if affector_num not in allele_to_construct_hash:
                        allele_to_construct_hash[affector_num] = set()
                    allele_to_construct_hash[affector_num].add(construct_num)

            # ### finish looping through fish file

        # given the components of a fish,
        # subtract out the intrinsic parts to just leave the extrinsic
        # to create the extrinsic genotypes.

        for fish_num in self.fish_parts:
            if self.test_mode and fish_num not in self.test_ids['fish']:
                continue

            fish_id = 'ZFIN:' + fish_num
            fish = self.fish_parts[fish_num]

            # get the intrinsic parts
            intrinsic_genotype_num = fish['intrinsic_genotype']
            intrinsic_genotype_id = 'ZFIN:' + intrinsic_genotype_num
            intrinsic_genotype_label = self.id_label_map.get(
                intrinsic_genotype_id)
            if intrinsic_genotype_num not in self.geno_alleles:
                intrinsic_parts = set()
            else:
                intrinsic_parts = self.geno_alleles[intrinsic_genotype_num]

            # subtract out the intrinsic parts, to get the extrinsic parts
            extrinsic_parts = fish['affectors'] - intrinsic_parts
            extrinsic_list = list(sorted(extrinsic_parts))

            # build up the extrinsic genotype from it's parts.
            # these will be reagents/morphants.
            if len(extrinsic_list) > 0:
                list_of_targeted_genes = []
                gene_to_reagent_hash = {}
                for eid in extrinsic_list:
                    # link the morpholino to the genes that it affects
                    eid = 'ZFIN:' + eid

                    # just in case, skip over the ALTs
                    if re.search(r'ALT', eid):
                        continue
                    ag = self.variant_loci_genes.get(eid)

                    # LOG.debug("%s affected genes %s", eid, str(ag))
                    if ag is None:
                        pass
                        # LOG.warn("No affected genes for %s", eid)
                    else:
                        # turn the gene-targeting-reagents inside out,
                        # such that instead of morph -> set(genes)
                        # we make a gene -> set(morphs)

                        for gid in ag:
                            if gid not in gene_to_reagent_hash:
                                gene_to_reagent_hash[gid] = set()
                            gene_to_reagent_hash[gid].add(eid)
                    # end loop through each extrinsic component

                for gid in gene_to_reagent_hash:
                    reagent_list = sorted(list(gene_to_reagent_hash.get(gid)))
                    # create variant gene(s) that have been targeted
                    # by the reagent

                    if gid not in self.id_label_map:
                        # should not happen, except maybe in testing
                        LOG.error("%s not in id-label-hash", gid)
                        glabel = gid
                    else:
                        glabel = self.id_label_map[gid]

                    eid = '-'.join(reagent_list)

                    targeted_gene_id = self.make_targeted_gene_id(
                        gid, eid)
                    # get the reagent labels
                    elabel = ', '.join(
                        self.id_label_map.get(l) for l in reagent_list)
                    if elabel is None:
                        elabel = eid  # should not happen, but just in case
                    targeted_gene_label = glabel + '<' + elabel + '>'

                    for r in reagent_list:
                        geno.addReagentTargetedGene(r, gid, targeted_gene_id,
                                                    targeted_gene_label)
                    self.id_label_map[targeted_gene_id] = targeted_gene_label
                    list_of_targeted_genes += [targeted_gene_id]
                    # end loop through each gene that is targeted
                list_of_targeted_genes = sorted(list_of_targeted_genes)
                extrinsic_id = self.make_id('-'.join(list_of_targeted_genes), '_')
                extrinsic_label = '; '.join([
                    str(self.id_label_map.get(l))
                    for l in list_of_targeted_genes])
                self.id_label_map[extrinsic_id] = extrinsic_label

                # add the parts
                for tg in list_of_targeted_genes:
                    if tg != extrinsic_id:
                        geno.addParts(
                            tg, extrinsic_id, self.globaltt['has_variant_part'])

            else:
                extrinsic_id = None
                extrinsic_label = None
            if extrinsic_id is not None:
                geno.addGenotype(
                    extrinsic_id, extrinsic_label, self.globaltt['extrinsic_genotype'])
                geno.addParts(
                    extrinsic_id, fish_id, self.globaltt['has_variant_part']
                )

            # check if the intrinsic is in the wildtype genotypes,
            # then it's a genomic background
            if intrinsic_genotype_id in self.wildtype_genotypes:
                intrinsic_rel = self.globaltt['has_reference_part']
                intrinsic_type = self.globaltt['genomic_background']
            else:
                intrinsic_rel = self.globaltt['has_variant_part']
                intrinsic_type = self.globaltt['intrinsic genotype']
            geno.addGenotype(
                intrinsic_genotype_id, intrinsic_genotype_label, intrinsic_type)

            # add the intrinsic to the fish
            geno.addParts(intrinsic_genotype_id, fish_id, intrinsic_rel)

            # build the fish label
            if extrinsic_id is None:
                fish_label = intrinsic_genotype_label
            else:
                fish_label = '; '.join((
                    str(intrinsic_genotype_label), extrinsic_label))

            fish_type = self.globaltt['effective_genotype']

            geno.addGenotype(fish_id, intrinsic_genotype_label, fish_type)
            geno.addTaxon(taxon_id, fish_id)

            # since we re-create a label,
            # add the zfin fish label as the synonym
            if fish['fish_label'] != '':
                model.addSynonym(fish_id, fish['fish_label'])
            self.id_label_map[fish_id] = fish_label

            if not self.test_mode and limit is not None and reader.line_num > limit:
                break

            # ###finish iterating over fish

        # iterate of the alleles and attache the constructs to them
        LOG.info("Adding Allele/Construct relationships")
        for a in allele_to_construct_hash:
            if self.test_mode and a not in self.test_ids['allele']:
                continue
            allele_id = 'ZFIN:' + a
            constructs = allele_to_construct_hash.get(a)
            if len(constructs) > 0:
                for c in constructs:
                    cid = 'ZFIN:' + c
                    geno.addSequenceDerivesFrom(allele_id, cid)
                    # LOG.info("constructs for %s: %s", allele_id,
                    #             str(constructs))
                    # migrate the transgenic features to be alternate parts
                    # of the transgene insertion/alteration
                    if cid in self.transgenic_parts:
                        tg_parts = self.transgenic_parts.get(cid)
                        if tg_parts is not None:
                            for p in tg_parts:
                                # HACK - if it's a promoter part,
                                # then make it a simple has_part
                                if re.search(r'promoter', p):
                                    r = self.globaltt['has_part']
                                else:
                                    r = self.globaltt['has_variant_part']
                                geno.addParts(p, allele_id, r)

    def _process_genotype_features(self, limit=None):
        """
        Here we process the genotype_features file, which lists genotypes
        together with any intrinsic sequence alterations, their zygosity,
        and affected gene.
        Because we don't necessarily get allele pair (VSLC) ids
        in a single row, we iterate through the file and build up a hash
        that contains all of a genotype's partonomy.
        We then assemble a genotype based on that partonomy.
        This table does not list the background genotype/strain:
        that is listed elsewhere.

        ZFIN "ALT" objects are mapped to sequence alterations in our system.

        By the end of this method, we have built up the intrinsic genotype,
        with Monarch-style labels.
        Original ZFIN labels are added (including the "sup" html tags), and
        Monarch-style labels are added as synonyms.

        We make assumptions here that any variants that affect the same locus
        are in trans.
        All of the genotype parts are created as BNodes at the moment,
        to avoid minting new Monarch ids, which means for anyone consuming this
        data they are inherently unstable.  This may change in the future.

        :param limit:
        :return:

        """
        src_key = 'geno'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing Genotypes from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        taxon_id = self.default_taxon_id
        geno_hash = {}  # This is used to store the genotype partonomy
        gvc_hash = {}
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                genotype_num = row[col.index('Genotype ID')].strip()
                genotype_name = row[col.index('Genotype Name')]
                genotype_unique_name = row[col.index('Genotye Unique Name')]
                allele_num = row[col.index('Allele ID')].strip()
                allele_name = row[col.index('Allele Name')]
                allele_ab = row[col.index('Allele Abbreviation')]
                allele_type = row[col.index('Allele Type')]
                # allele_disp_type = row[col.index('Allele Display Type')]
                gene_symbol = row[col.index('Gene or Construct Symbol')]
                gene_num = row[col.index(
                    'Corresponding ZFIN Gene ID/Construct ID')].strip()
                zygosity = row[col.index('Allele Zygosity')]
                construct_name = row[col.index('Construct Name')]
                construct_num = row[col.index('Construct ZdbId')].strip()

                if self.test_mode and genotype_num not in self.test_ids['genotype']:
                    continue

                # add the genotype to the graph
                # not adding the genotype label here,
                # since it doesn't include the background
                # that will be done in another method
                genotype_curie = 'ZFIN:' + genotype_num
                geno.addGenotype(genotype_curie, None)

                # add the given name from ZFIN and uniquename as synonyms
                if genotype_name != '':
                    model.addSynonym(genotype_curie, genotype_name)
                if genotype_unique_name != '':
                    model.addSynonym(genotype_curie, genotype_unique_name)

                # store the alleles of the genotype,
                # in order to use when processing fish
                if genotype_num not in self.geno_alleles:
                    self.geno_alleles[genotype_num] = set()
                self.geno_alleles[genotype_num].add(allele_num)

                if genotype_curie not in geno_hash:
                    geno_hash[genotype_curie] = {}

                genoparts = geno_hash[genotype_curie]

                # reassign the allele_type to a proper GENO or SO class
                # allele_type = self._map_allele_type_to_geno(allele_type)
                allele_type_id = self.resolve(allele_type, False)
                if allele_type_id == allele_type:
                    allele_type_id = self.globaltt['unspecified']  # is geno: not zfa:

                allele_curie = ':'.join(('ZFIN', allele_num))

                if allele_num != '':
                    self.id_label_map[allele_curie] = allele_name
                else:
                    allele_curie = None

                # alleles in zfin are really sequence alterations in our system
                geno.addSequenceAlteration(allele_curie, allele_name, allele_type_id)
                if allele_ab != '':
                    model.addSynonym(allele_curie, allele_ab)

                # here, we assemble the items into a genotype hash
                # we need to do this because each row only holds one allele
                # of a gene but a genotype may have many alleles and therefore
                # many rows so we loop through the file once to build a hash of
                # genotype components
                if gene_num is not None and gene_num != '':
                    # add the gene to the graph, along with it's symbol
                    # as the primary label
                    gene_curie = ':'.join(('ZFIN', gene_num))
                    geno.addGene(gene_curie, gene_symbol)
                    self.id_label_map[gene_curie] = gene_symbol

                    # if it's a transgenic construct,
                    # then we'll have to get the other bits
                    if construct_num is not None and construct_num != '':
                        construct_curie = ':'.join(('ZFIN', construct_num))
                        geno.addSequenceDerivesFrom(allele_curie, construct_curie)
                        self.id_label_map[construct_curie] = construct_name

                    # allele to gene
                    if allele_curie not in self.variant_loci_genes:
                        self.variant_loci_genes[allele_curie] = [gene_curie]
                    elif gene_curie not in self.variant_loci_genes[allele_curie]:
                        self.variant_loci_genes[allele_curie] += [gene_curie]

                    if gene_curie not in genoparts:
                        genoparts[gene_curie] = [allele_curie]
                    else:
                        genoparts[gene_curie] += [allele_curie]

                    other_allele = self._get_other_allele_by_zygosity(
                        allele_curie, zygosity)
                    if other_allele is not None:
                        genoparts[gene_curie] += [other_allele]

                else:
                    # if the gene is not known,
                    # still need to add the allele to the genotype hash
                    # these will be added as sequence alterations.
                    genoparts[allele_curie] = [allele_curie]
                    other_allele = self._get_other_allele_by_zygosity(
                        allele_curie, zygosity)
                    if other_allele is not None:
                        genoparts[allele_curie] += [other_allele]

                geno_hash[genotype_curie] = genoparts

                # fetch the other affected genes,
                # and make sure they are in the geno hash
                # we have to do this because some of the affected genes
                # are not listed in this file
                genes_from_hash = None
                if allele_curie in self.variant_loci_genes:
                    genes_from_hash = self.variant_loci_genes[allele_curie]
                else:
                    pass
                    # LOG.info('no gene found for %s', allele_curie)

                if genes_from_hash is not None \
                        and genes_from_hash != [gene_curie] \
                        and gene_curie not in genes_from_hash:
                    LOG.info(
                        "***Found genes not in genotype_features for %s: %s",
                        allele_curie, genes_from_hash)
                    for gh in genes_from_hash:
                        if gh not in genoparts:
                            genoparts[gh] = [allele_curie]
                        else:
                            genoparts[gh] += [allele_curie]

                        other_allele = self._get_other_allele_by_zygosity(
                            allele_curie, zygosity)
                        if other_allele is not None:
                            genoparts[gh].append(other_allele)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Finished parsing file: %s", raw)

        # ############## BUILD THE INTRINSIC GENOTYPES ###############
        # using the geno_hash, build the genotype parts,
        # and add them to the graph
        # the hash is organized like:
        # genotype_id : {
        #   gene_curie : [list, of, alleles], # for located things
        #   allele_curie : [list, of, alleles] # for unlocated things
        #   }
        # now loop through the geno_hash, and build the vslcs

        LOG.info("Building intrinsic genotypes from partonomy")
        for gt in geno_hash:
            if self.test_mode and re.sub(r'ZFIN:', '', gt) \
                    not in self.test_ids['genotype']:
                print('skipping ', gt)
                continue

            if gt not in gvc_hash:
                gvc_hash[gt] = []
            gvcparts = gvc_hash[gt]

            for locus_id in geno_hash[gt]:
                # LOG.info("locus id %s",locus_id)
                locus_label = self.id_label_map[locus_id]
                variant_locus_parts = geno_hash.get(gt).get(locus_id)
                # LOG.info(
                #   'vl parts: %s',pprint.pformat(variant_locus_parts))
                # if the locus == part, then it isn't a gene,
                # rather a variant not in a specific gene
                if locus_id in variant_locus_parts:
                    # set the gene_curie to none
                    gene_curie = None
                else:
                    gene_curie = locus_id

                allele1_id = variant_locus_parts[0]
                if allele1_id not in self.id_label_map:
                    allele1_label = allele1_id
                    LOG.error('allele1 %s not in hash', allele1_id)
                else:
                    allele1_label = self.id_label_map[allele1_id]
                allele2_id = None
                allele2_label = None
                zygosity_id = None

                if len(variant_locus_parts) > 2:
                    LOG.error(
                        "There may be a problem. >2 parts for this locus (%s): %s",
                        locus_id, variant_locus_parts)
                elif len(variant_locus_parts) > 1:
                    allele2_id = variant_locus_parts[1]
                    if allele2_id not in ['0', '?']:
                        allele2_label = self.id_label_map[allele2_id]
                    else:
                        allele2_label = allele2_id
                if allele2_id is not None:
                    if allele2_id == '?':
                        zygosity_id = self.globaltt['indeterminate']
                        allele2_id = 'UN'
                    elif allele2_id == '0':
                        zygosity_id = self.globaltt['hemizygous']
                    elif allele1_id != allele2_id:
                        zygosity_id = self.globaltt['compound heterozygous']
                    elif allele1_id == allele2_id:
                        zygosity_id = self.globaltt['homozygous']
                else:
                    zygosity_id = self.globaltt['simple heterozygous']
                    allele2_label = '+'
                    allele2_id = 'WT'

                # make variant_loci
                vloci2 = vloci2_label = None
                if gene_curie is not None:
                    vloci1 = self._make_variant_locus_id(gene_curie, allele1_id)
                    vloci1_label = geno.make_variant_locus_label(
                        locus_label, allele1_label)
                    geno.addSequenceAlterationToVariantLocus(
                        allele1_id, vloci1)
                    geno.addAlleleOfGene(vloci1, gene_curie)
                    model.addIndividualToGraph(
                        vloci1, vloci1_label, self.globaltt['variant_locus'])
                    if allele2_id is not None and allele2_id not in ['WT', '0', 'UN']:
                        vloci2 = self._make_variant_locus_id(
                            gene_curie, allele2_id)
                        vloci2_label = geno.make_variant_locus_label(
                            locus_label, allele2_label)
                        geno.addSequenceAlterationToVariantLocus(
                            allele2_id, vloci2)
                        model.addIndividualToGraph(
                            vloci2, vloci2_label, self.globaltt['variant_locus']
                        )
                        geno.addAlleleOfGene(vloci2, gene_curie)
                else:
                    vloci1 = allele1_id
                    vloci1_label = allele1_label
                    vloci2 = None
                    if allele2_id not in ['WT', '0', 'UN']:
                        vloci2 = allele2_id
                    vloci2_label = allele2_label

                # create the vslc
                gene_label = ''
                if gene_curie is None:
                    gn = 'UN'
                else:
                    gn = gene_curie
                    gene_label = self.id_label_map[gene_curie]

                # TODO also consider adding this to Genotype.py
                vslc_id = '-'.join((gn, allele1_id, allele2_id))
                vslc_id = self.make_id(re.sub(r'(ZFIN)?:', '', vslc_id), '_')

                vslc_label = geno.make_vslc_label(
                    gene_label, allele1_label, allele2_label)

                # add to global hash
                self.id_label_map[vslc_id] = vslc_label

                model.addIndividualToGraph(
                    vslc_id,
                    vslc_label,
                    self.globaltt['variant single locus complement']
                )
                geno.addPartsToVSLC(
                    vslc_id, vloci1, vloci2, zygosity_id,
                    self.globaltt['has_variant_part'],
                    self.globaltt['has_variant_part']
                )

                gvcparts += [vslc_id]

            gvc_hash[gt] = gvcparts

        # end loop through geno_hash

        LOG.info('Finished finding all the intrinsic genotype parts')
        LOG.info('Build pretty genotype labels')
        # now loop through the gvc_hash, and build the gvc
        for gt in gvc_hash:
            if self.test_mode and re.sub(r'ZFIN:', '', gt) \
                    not in self.test_ids['genotype']:
                continue

            gvc_parts = gvc_hash[gt]

            # only need to make a gvc specifically if there's >1 vslc
            if len(gvc_parts) > 1:
                gvc_labels = []
                # put these in order so they will always make the same id
                gvc_parts.sort()

                gvc_id = '-'.join(gvc_parts)
                gvc_id = re.sub(r'(ZFIN)?:', '', gvc_id)
                gvc_id = self.make_id(re.sub(r'^_*', '', gvc_id), '_')

                for vslc_id in gvc_parts:
                    # add the vslc to the gvc
                    geno.addVSLCtoParent(vslc_id, gvc_id)

                    # build the gvc label
                    vslc_label = self.id_label_map[vslc_id]
                    if vslc_label is not None:
                        gvc_labels += [vslc_label]
                    else:
                        gvc_labels += [vslc_id]

                gvc_labels.sort()
                gvc_label = '; '.join(gvc_labels)

                # add the gvc to the id-label hash
                self.id_label_map[gvc_id] = gvc_label

                # add the gvc
                model.addIndividualToGraph(
                    gvc_id, gvc_label, self.globaltt['genomic_variation_complement']
                )
            elif len(gvc_parts) == 1:
                # assign the vslc to be also a gvc
                vslc_id = gvc_parts[0]
                gvc_id = vslc_id
                gvc_label = self.id_label_map[vslc_id]
                model.addType(vslc_id, self.globaltt['genomic_variation_complement'])
            else:
                gvc_id = None
                gvc_label = ''
                LOG.error("No GVC parts for %s", gt)

            if gt in self.genotype_backgrounds:
                background_id = self.genotype_backgrounds[gt]
                if background_id in self.id_label_map:
                    background_label = self.id_label_map[background_id]
                else:
                    background_label = background_id
                    LOG.error("We don't have the label for %s stored", background_id)

            else:
                background_num = re.sub(r'ZFIN:', '', gt)
                background_id = self.make_id('bkgd-' + background_num, '_')
                background_label = 'unspecified background'
                background_label_and_num = \
                    'unspecified background (' + background_num + ')'
                background_desc = 'This genomic background is unknown. ' +\
                    'This is a placeholder background for ' + gt + '.'
                # there is no background for this genotype;
                # need to add the taxon to this one!
                # make an anonymous background for this genotype
                geno.addGenomicBackground(
                    background_id, background_label_and_num, None, background_desc)
                geno.addGenomicBackgroundToGenotype(background_id, gt)

            geno.addTaxon(taxon_id, background_id)

            # add genotype_name + background as label for this gt
            genotype_name_background = genotype_name + ' (' + background_label + ')'
            geno.addGenotype(gt, genotype_name_background)

            # add monarch_genotype_name as synonym
            monarch_genotype_name = gvc_label + ' [' + background_label + ']'
            if monarch_genotype_name != '':
                model.addSynonym(gt, monarch_genotype_name)
                self.id_label_map[gt] = monarch_genotype_name

            # Add the GVC to the genotype
            geno.addParts(gvc_id, gt, self.globaltt['has_variant_part'])

            # end of gvc loop

        # end of genotype loop

        # TODO this is almost complete;
        # deficiencies with >1 locus deleted are still not right
        LOG.info("Finished building genotype labels")
        LOG.info("Done with genotypes")

    def _process_genotype_backgrounds(self, limit=None):
        """
        This table provides a mapping of genotypes to background genotypes
        Note that the background_id is also a genotype_id.

        Makes these triples:
        <ZFIN:genotype_id> GENO:has_reference_part <ZFIN:background_id>
        <ZFIN:background_id> a GENO:genomic_background
        <ZFIN:background_id> in_taxon <taxon_id>
        <taxon_id> a class
        :param limit:
        :return:

        """
        src_key = 'backgrounds'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing genotype backgrounds from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)

        # Add the taxon as a class
        taxon_id = self.default_taxon_id
        model.addClassToGraph(taxon_id, None)

        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                genotype_id = row[col.index('Genotype ID')].strip()
                # genotype_name = row[col.index('Genotype Name')]
                background_id = row[col.index('Background')].strip()
                # background_name = row[col.index('Background Name')]

                if self.test_mode and genotype_id not in self.test_ids['genotype']:
                    continue

                genotype_curie = ':'.join(('ZFIN', genotype_id))
                bg_curie = 'ZFIN:' + background_id

                # store this in the hash for later lookup
                # when building fish genotypes
                self.genotype_backgrounds[genotype_curie] = bg_curie

                # add the background into the graph,
                # in case we haven't seen it before
                geno.addGenomicBackground(bg_curie, None)

                # hang the taxon from the background
                geno.addTaxon(taxon_id, bg_curie)

                # add the intrinsic genotype to the graph
                # we DO NOT ADD THE LABEL here
                # as it doesn't include the background
                geno.addGenotype(
                    genotype_curie, None, self.globaltt['intrinsic genotype'])

                # Add background to the intrinsic genotype
                geno.addGenomicBackgroundToGenotype(bg_curie, genotype_curie)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with genotype backgrounds")

    def _process_wildtypes(self, limit=None):
        """
        This table provides the genotype IDs, name,
        and abbreviation of the wildtype genotypes.
        These are the typical genomic backgrounds...there's about 20 of them.
        http://zfin.org/downloads/wildtypes_fish.txt

        Triples created:
        <genotype id> a GENO:wildtype
        <genotype id> rdfs:label genotype_abbreviation
        <genotype id> dc:description genotype_name

        :param limit:
        :return:

        """
        src_key = 'wild'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing wildtype genotypes from file: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        # model = Model(graph)  # unused
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                fish_num = row[col.index('Fish ID')]
                fish_name = row[col.index('Fish Name')]
                fish_abbreviation = row[col.index('Fish Abbreviation')]
                genotype_num = row[col.index('Genotype ID')].strip()

                # ZDB-FISH-150901-10750 INDO    INDO    ZDB-GENO-980210-32
                fish_id = 'ZFIN:' + fish_num
                genotype_id = 'ZFIN:' + genotype_num
                background_type = self.globaltt['genomic_background']

                # Add genotype to graph with label and description,
                # as a genomic_background genotype
                unspecified_background = 'ZDB-GENO-030619-2'
                if re.match(genotype_num.strip(), unspecified_background):
                    background_type = self.globaltt['unspecified_genomic_background']
                geno.addGenomicBackground(
                    genotype_id, fish_abbreviation, background_type, fish_name)

                graph.addTriple(fish_id, self.globaltt['has_genotype'], genotype_id)

                # Build the hash for the wild type genotypes.
                self.id_label_map[genotype_id] = fish_abbreviation

                # store these in a special hash to look up later
                self.wildtype_genotypes += [genotype_id]

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with wildtype genotypes")

    def _process_stages(self, limit=None):
        """
        This table provides mappings between ZFIN stage IDs and ZFS terms,
        and includes the starting and ending hours for the developmental stage.
        Currently only processing the mapping from the ZFIN stage ID
        to the ZFS ID.

        :param limit:
        :return:

        """
        src_key = 'stage'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing stages from file: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)

        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)
                stage_id = row[col.index('Stage ID')].strip()
                stage_obo_id = row[col.index('Stage OBO ID')]
                stage_name = row[col.index('Stage Name')]
                # begin_hours = row[col.index('Begin Hours')]
                # end_hours = row[col.index('End Hours')]

                # Add the stage as a class, and it's obo equivalent
                stage_curie = ':'.join(('ZFIN', stage_id))
                model.addClassToGraph(stage_curie, stage_name,
                                      class_category=blv.terms['LifeStage'])
                model.addEquivalentClass(stage_curie, stage_obo_id,
                                         subject_category=blv.terms['LifeStage'],
                                         object_category=blv.terms['LifeStage'])

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with stages")

    def _process_g2p(self, limit=None):
        """
        Here, we process the fish-to-phenotype associations,
        which also include environmental perturbations.
        The phenotypes may also be recorded as observed at specific stages.
        We create association objects with as much of the information
        as possible.

        A full association object may include:

        assoc hasSubject effective genotype
        assoc hasObject phenotype
        assoc hasSource publication (PMID or ZFIN pub id)
        assoc hasEnvironment environment
        assoc hasStartStage start_stage_id
        assoc hasEndStage end_stage_id

        :param limit:
        :return:

        """

        src_key = 'pheno'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing G2P from file: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        mapped_zpids = list()
        model = Model(graph)
        eco_id = self.globaltt['experimental phenotypic evidence']
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                fish_num = row[col.index('Fish ID')].strip()
                # fish_name = row[col.index('Fish Name')]
                start_stage_id = row[col.index('Start Stage ID')].strip()
                # start_stage_name = row[col.index('Start Stage Name')]
                end_stage_id = row[col.index('End Stage ID')].strip()
                # end_stage_name = row[col.index('End Stage Name')]

                subterm1_id = row[col.index(
                    'Affected Structure or Process 1 subterm ID')]
                subterm1_name = row[col.index(
                    'Affected Structure or Process 1 subterm Name')]
                postcomp1_rel_id = row[col.index('Post-composed Relationship ID')]
                # postcomp1_rel_name = row[col.index('Post-composed Relationship Name')]
                superterm1_id = row[col.index(
                    'Affected Structure or Process 1 superterm ID')]
                superterm1_name = row[col.index(
                    'Affected Structure or Process 1 superterm Name')]
                quality_id = row[col.index('Phenotype Keyword ID')]
                quality_name = row[col.index('Phenotype Keyword Name')]
                modifier = row[col.index('Phenotype Tag')]
                subterm2_id = row[col.index(
                    'Affected Structure or Process 2 subterm ID')]
                subterm2_name = row[col.index(
                    'Affected Structure or Process 2 subterm name')]
                postcomp2_rel_id = row[col.index('Post-composed Relationship (rel) ID')]
                # postcomp2_rel_name = row[col.index(
                #    'Post-composed Relationship (rel) Name')]
                superterm2_id = row[col.index(
                    'Affected Structure or Process 2 superterm ID')]
                superterm2_name = row[col.index(
                    'Affected Structure or Process 2 superterm name')]
                pub_id = row[col.index('Publication ID')].strip()
                env_id = row[col.index('Environment ID')].strip()

                if self.test_mode and (
                        fish_num not in self.test_ids['fish'] or
                        env_id not in self.test_ids['environment']):
                    continue

                fish_curie = ':'.join(('ZFIN', fish_num))
                env_curie = ':'.join(('ZFIN', env_id))

                # ########### PHENOTYPES ##########
                phenotype_id = self._map_octuple_to_phenotype(
                    subterm1_id, postcomp1_rel_id, superterm1_id, quality_id,
                    subterm2_id, postcomp2_rel_id, superterm2_id, modifier)

                if phenotype_id is not None:
                    mapped_zpids.append([
                        subterm1_id, postcomp1_rel_id, superterm1_id, quality_id,
                        subterm2_id, postcomp2_rel_id, superterm2_id, modifier])

                if pub_id != '':
                    pub_id = ':'.join(('ZFIN', pub_id))
                    ref = Reference(graph, pub_id)
                    ref.addRefToGraph()

                if modifier[:6] != 'normal':
                    if phenotype_id is None:
                        continue
                    if start_stage_id != '':
                        start_stage_id = ':'.join(('ZFIN', start_stage_id))
                    if end_stage_id != '':
                        end_stage_id = ':'.join(('ZFIN', end_stage_id))

                    # add association
                    assoc = G2PAssoc(graph, self.name, fish_curie, phenotype_id)

                    # only add the environment if there's components to it
                    if env_curie in self.environment_hash and len(
                            self.environment_hash.get(env_curie)) > 0:
                        assoc.set_environment(env_curie)
                    assoc.set_stage(start_stage_id, end_stage_id)
                    assoc.add_evidence(eco_id)
                    assoc.add_source(pub_id)
                    assoc.add_association_to_graph()
                    assoc_id = assoc.get_association_id()

                    if env_curie not in self.environment_hash or len(
                            self.environment_hash.get(env_curie)) > 0:
                        model.addComment(assoc_id, 'Legacy environment id ' + env_curie)
                else:
                    # TODO add normal phenotypes as associations #134 when
                    # https://github.com/sba1/bio-ontology-zp/issues/9
                    # is finished, we can use these add normal phenotypes
                    # as a comment on the genotype for now
                    clist = []
                    for x in [superterm1_name, subterm1_name, quality_name,
                              superterm2_name, subterm2_name, modifier]:
                        if x != '':
                            clist += [x]

                    c = '+'.join(clist)
                    c = ' '.join(("Normal phenotype observed:", c, "(" + pub_id + ")"))
                    if pub_id != '':
                        graph.addTriple(pub_id, self.globaltt['mentions'], fish_curie)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        self.mapped_zpids = mapped_zpids
        myset = set([','.join(x) for x in mapped_zpids])
        LOG.info("Phenotype-octuples: %i mapped", len(myset))

    def _process_genes(self, limit=None):
        """
        This table provides the ZFIN gene id, the SO type of the gene,
        the gene symbol, and the NCBI Gene ID.

        Triples created:
        <gene id> a class
        <gene id> rdfs:label gene_symbol
        <gene id> equivalent class <ncbi_gene_id>
        <gene id> in_taxon <ncbitaxon_curie>
        :param limit:
        :return:

        """
        src_key = 'gene'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing genes from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        taxon_id = self.default_taxon_id
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                gene_id = row[col.index('ZFIN ID')]
                # gene_so_id = row[col.index('SO ID')]
                gene_symbol = row[col.index('Symbol')]
                ncbi_gene_id = row[col.index('NCBI Gene ID')]

                if self.test_mode and gene_id not in self.test_ids['gene']:
                    continue

                gene_id = 'ZFIN:' + gene_id.strip()
                ncbi_gene_id = 'NCBIGene:' + ncbi_gene_id.strip()

                self.id_label_map[gene_id] = gene_symbol

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

                geno.addGene(gene_id, gene_symbol)
                model.addEquivalentClass(gene_id, ncbi_gene_id)
                geno.addTaxon(taxon_id, gene_id)

        LOG.info("Done with genes")

    def _process_features(self, limit=None):
        """
        This module provides information for the intrinsic
        and extrinsic genotype features of zebrafish.
        All items here are 'alterations', and are therefore instances.

        sequence alteration ID, SO type, abbreviation, and relationship to
        the affected gene, with the gene's ID, symbol,
        and SO type (gene/pseudogene).

        Triples created:
        <gene id> a class:
        :param limit:
        :return:

        """
        src_key = 'features'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing features from filr: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                genomic_feature_id = row[col.index('Genomic Feature ID')].strip()
                feature_so_id = row[col.index('Feature SO ID')]
                genomic_feature_abbreviation = row[col.index(
                    'Genomic Feature Abbreviation')]
                genomic_feature_name = row[col.index('Genomic Feature Name')]
                # genomic_feature_type = row[col.index('Genomic Feature Type')]
                # mutagen = row[col.index('Mutagen')]
                # mutagee = row[col.index('Mutagee')]
                construct_id = row[col.index('Construct ID')].strip()
                construct_name = row[col.index('Construct name')]
                construct_so_id = row[col.index('Construct SO ID')]
                # talen_crispr_id = row[col.index('TALEN/CRISPR ID')]
                # talen_crispr_name = row[col.index('TALEN/CRISPR Name')]

                if self.test_mode and (
                        genomic_feature_id not in self.test_ids['allele']):
                    continue

                genomfeat_curie = ':'.join(('ZFIN', genomic_feature_id))
                model.addIndividualToGraph(
                    genomfeat_curie, genomic_feature_name, feature_so_id
                )
                if genomic_feature_abbreviation != '':
                    model.addSynonym(genomfeat_curie, genomic_feature_abbreviation)

                if construct_id is not None and construct_id != '':
                    construct_curie = ':'.join(('ZFIN', construct_id))
                    geno.addConstruct(construct_curie, construct_name, construct_so_id)
                    geno.addSequenceDerivesFrom(genomfeat_curie, construct_curie)
                    self.id_label_map[construct_curie] = construct_name

                # Note, we don't really care about how the variant was derived.
                # so we skip that.
                # add to the id-label map
                self.id_label_map[genomfeat_curie] = genomic_feature_abbreviation

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with features")

    def _process_feature_affected_genes(self, limit=None):
        """
        This table lists (intrinsic) genomic sequence alterations
        and their affected gene(s).
        It provides the sequence alteration ID, SO type, abbreviation,
        and relationship to the affected gene, with the gene's ID, symbol,
        and SO type (gene/pseudogene).

        Triples created:
        <gene id> a class:
        <gene id> rdfs:label gene_symbol
        <gene id> subclass of gene/pseudogene

        <variant locus id> is a GENO:allele
        <variant locus id> rdfs:label <variant_locus_label>
        <variant locus id> is an allele of <gene id>
        <variant locus id> has alternate part <sequence alteration id>

        <sequence alteration id> is an allele of <gene id>
        <sequence alteration id> rdf:type <sequence alteration type>

        :param limit:
        :return:

        """

        # can use this to process and build the variant locus.
        # but will need to process through some kind of helper hash,
        # just like we do with the genotype file.
        # that's because each gene is on a separate line
        # for example, ZDB-ALT-021021-2 is a deficiency that affects 4 genes
        # that case is when the relationship is != 'is allele of'

        src_key = 'feature_affected_gene'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing feature affected genes from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                genomic_feature_id = row[col.index('Genomic Feature ID')].strip()
                feature_so_id = row[col.index('Feature SO ID')]
                genomic_feature_abbreviation = row[col.index(
                    'Genomic Feature Abbreviation')].strip()
                gene_symbol = row[col.index('Gene Symbol')]
                gene_id = row[col.index('Gene ID')].strip()
                gene_so_id = row[col.index('Gene SO ID')]
                # note: there are ~27 more columns in this file

                # alteration Seq types present in column 2 of the file (2020 Jun)
                # cut -f2 raw/zfin/features-affected-genes.txt|sort|uniq -c|sort -nr
                #  38087 SO:1000008  - point_mutation
                #   6505 SO:0001218  - transgenic insertion
                #   3433 SO:0001060  - sequence_variant
                #   2364 SO:0000159  - sequence_alteration
                #    833 SO:1000032  - indel
                #    424 SO:0000667  - insertion
                #    191 SO:1000005  - complex_substitution
                #    118 SO:1000029  - chromosomal_deletion
                #     95 SO:0001059  - deletion
                #     43 SO:0000199  - translocation
                # gene Seq types present in column 6 of the file:
                # cut -f6 raw/zfin/features-affected-genes.txt|sort|uniq -c|sort -nr
                #  51942 SO:0001217
                #     65 SO:0001265
                #     58 SO:0001641
                #     22 SO:0000336
                #      5 SO:0000165
                #      1 SO:0005836

                if self.test_mode and (
                        gene_id not in self.test_ids['gene'] and
                        genomic_feature_id not in self.test_ids['allele']):
                    continue

                genomfeat_curie = ':'.join(('ZFIN', genomic_feature_id))
                gene_curie = ':'.join(('ZFIN', gene_id))

                self.id_label_map[genomfeat_curie] = genomic_feature_abbreviation
                self.id_label_map[gene_curie] = gene_symbol

                geno.addGene(gene_curie, gene_symbol, gene_so_id)

                # add the gene to the list of things altered by this thing
                if genomfeat_curie not in self.variant_loci_genes:
                    self.variant_loci_genes[genomfeat_curie] = [gene_curie]
                elif gene_curie not in self.variant_loci_genes[genomfeat_curie]:
                    self.variant_loci_genes[genomfeat_curie] += [gene_curie]

                # Add the sequence alteration id, label, and type(-as-curie)
                geno.addSequenceAlteration(
                    genomfeat_curie, genomic_feature_abbreviation, feature_so_id)

                if self.globaltcid[feature_so_id] == 'is allele of':
                    vl_label = geno.make_variant_locus_label(
                        gene_symbol, genomic_feature_abbreviation)

                    vl_id = self._make_variant_locus_id(gene_curie, genomfeat_curie)

                    self.id_label_map[vl_id] = vl_label

                    # create the variant locus,
                    # add it's parts and relationship to the gene
                    geno.addSequenceAlterationToVariantLocus(genomic_feature_id, vl_id)
                    model.addIndividualToGraph(
                        vl_id, vl_label, self.globaltt['variant_locus']
                    )
                    geno.addAlleleOfGene(vl_id, gene_curie)

                    # note that deficiencies or translocations
                    # that affect only one gene are considered alleles here
                    # by zfin, which is appropriate.
                    # I don't yet see duplications
                else:
                    # don't make the variant loci for the other things
                    # which include deficiencies, translocations, transgenes
                    # TODO review this
                    pass

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with feature affected genes")

    def _process_gene_marker_relationships(self, limit=None):
        """
        Gene-marker relationships include:
            clone contains gene,
            coding sequence of,
            contains polymorphism,
            gene contains small segment,
            gene encodes small segment,
            gene has artifact,
            gene hybridized by small segment,
            gene produces transcript,
            gene product recognized by antibody,
            knockdown reagent targets gene,
            promoter of,
            transcript targets gene

        Here, we only process the following:
            knockdown reagent targets gene,
            coding sequence of,
            promoter of,
            transcript targets gene

        We only take a fraction of these here...
        we are interested in the knockdown reagents, promoters, and
        the transgenic constructs with coding bits.

        :param limit:
        :return:

        """
        src_key = 'gene_marker_rel'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing gene marker relationships from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)
                gene_id = row[col.index('Gene ID')].strip()
                gene_so_id = row[col.index('Gene SO ID')]
                gene_symbol = row[col.index('Gene Symbol')]
                marker_id = row[col.index('Marker ID')].strip()
                marker_so_id = row[col.index('Marker SO ID')]
                marker_symbol = row[col.index('Marker Symbol')]
                relationship = row[col.index('Relationship')]

                if self.test_mode and not (
                        gene_id in self.test_ids['gene'] or
                        marker_id in self.test_ids['allele'] or
                        marker_id in self.test_ids['morpholino']):
                    continue

                # there are many relationships, but we only take a few for now
                if relationship not in [
                        'knockdown reagent targets gene',
                        'coding sequence of',
                        'gene product recognized by antibody',
                        'promoter of',
                        'transcript targets gene']:
                    continue

                gene_curie = ':'.join(('ZFIN', gene_id))
                geno.addGene(gene_curie, gene_symbol, gene_so_id)

                marker_curie = 'ZFIN:' + marker_id
                if relationship == 'knockdown reagent targets gene':
                    geno.addGeneTargetingReagent(
                        marker_curie, marker_symbol, marker_so_id, gene_curie)
                    # waiting to add the reagent_targeted_gene
                    # until processing environments
                elif relationship == 'coding sequence of':
                    # we add the partonomy
                    # directly to the allele in the process_fish method
                    geno.addConstruct(
                        marker_curie, marker_symbol, marker_so_id)
                    transgene_part_id = self._make_transgene_part_id(
                        marker_curie, gene_curie, relationship)
                    transgene_part_label = 'Tg(' + relationship + ' ' +\
                        gene_symbol + ')'
                    model.addIndividualToGraph(
                        transgene_part_id,
                        transgene_part_label,
                        self.globaltt['coding_transgene_feature']
                    )
                    geno.addSequenceDerivesFrom(transgene_part_id, gene_curie)

                    # save the transgenic parts in a hashmap for later
                    if marker_curie not in self.transgenic_parts:
                        self.transgenic_parts[marker_curie] = set()
                    self.transgenic_parts[marker_curie].add(transgene_part_id)
                    self.id_label_map[transgene_part_id] = transgene_part_label

                elif relationship == 'gene product recognized by antibody':
                    # TODO for ticket #32
                    pass
                elif relationship == 'promoter of':
                    # transgenic constructs with promoters regions
                    # we add the partonomy
                    # directly to the allele in the process_fish method
                    geno.addConstruct(marker_curie, marker_symbol, marker_so_id)
                    transgene_part_id = self._make_transgene_part_id(
                        marker_curie, gene_curie, relationship)
                    transgene_part_label = 'Tg(' + relationship + ' ' +\
                        gene_symbol + ')'
                    model.addIndividualToGraph(
                        transgene_part_id, transgene_part_label,
                        self.globaltt['regulatory_transgene_feature'])
                    geno.addSequenceDerivesFrom(transgene_part_id, gene_curie)

                    # save the transgenic parts in a hashmap for later
                    if marker_curie not in self.transgenic_parts:
                        self.transgenic_parts[marker_curie] = set()
                    self.transgenic_parts[marker_curie].add(transgene_part_id)

                elif relationship == 'transcript targets gene':  # miRNAs
                    # TODO should this be an interaction
                    # instead of this special relationship?
                    model.addIndividualToGraph(
                        marker_curie, marker_symbol, marker_so_id
                    )
                    graph.addTriple(
                        marker_curie, self.globaltt['targets_gene'], gene_curie
                    )

                self.id_label_map[marker_curie] = marker_symbol
                # just in case we haven't seen it before
                self.id_label_map[gene_curie] = gene_symbol

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with gene marker relationships")

    @staticmethod
    def _make_transgene_part_id(construct_id, part_id, relationship):
        transgene_part_id = '-'.join((
            construct_id, part_id, re.sub(r'\W+', '-', relationship)))
        transgene_part_id = re.sub(r'ZFIN:', '', transgene_part_id)
        transgene_part_id = Source.make_id(transgene_part_id, '_')

        return transgene_part_id

    def _process_pubinfo(self, limit=None):
        """
        This will pull the zfin internal publication information,
        and map them to their equivalent pmid, and make labels.

        Triples created:
        <pub_curie> is an individual
        <pub_curie> rdfs:label <pub_label>
        <pubmed_id> is an individual
        <pubmed_id> rdfs:label <pub_label>

        <pub_curie> sameIndividual <pubmed_id>
        :param limit:
        :return:

        """
        src_key = 'pubs'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing Pubs info from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)

        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="latin-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                pub_id = row[col.index('Publication ID')].strip()
                pubmed_id = row[col.index(
                    'pubMed ID (none or blank when not available)')]
                authors = row[col.index('Authors')]
                title = row[col.index('Title')]
                journal = row[col.index('Journal')]
                year = row[col.index('Year')]
                vol = row[col.index('Volume')]
                pages = row[col.index('Pages')]

                pub_curie = ':'.join(('ZFIN', pub_id))
                if self.test_mode and (
                        pub_curie not in self.test_ids['pub'] and
                        'PMID:' + pubmed_id not in self.test_ids['pub']):
                    continue

                # trim the author list for ease of reading
                alist = re.split(r',', authors)
                if len(alist) > 1:
                    astring = ' '.join((alist[0].strip(), 'et al'))
                else:
                    astring = authors

                pub_label = '; '.join((astring, title, journal, year, vol, pages))
                ref = Reference(graph, pub_curie)
                ref.setShortCitation(pub_label)
                ref.setYear(year)
                ref.setTitle(title)

                if pubmed_id is not None and pubmed_id != '':
                    # let's make an assumption that if there's a pubmed id,
                    # that it is a journal article
                    ref.setType(self.globaltt['journal article'])

                    pubmed_id = 'PMID:' + pubmed_id.strip()
                    rpm = Reference(graph, pubmed_id, self.globaltt['journal article'])
                    rpm.addRefToGraph()

                    model.addSameIndividual(pub_curie, pubmed_id)
                    model.makeLeader(pubmed_id)

                ref.addRefToGraph()

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

    def _process_pub2pubmed(self, limit=None):
        """
        This will pull the zfin internal publication to pubmed mappings.
        Somewhat redundant with the process_pubinfo method,
        but this includes additional mappings.

        <pub_id> is an individual
        <pub_id> rdfs:label <pub_label>
        <pubmed_id> is an individual
        <pubmed_id> rdfs:label <pub_label>

        <pub_id> sameIndividual <pubmed_id>
        :param limit:
        :return:
        """
        src_key = 'pub2pubmed'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing Pub to PunMed from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="latin-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                pub_id = row[col.index('Publication ZFIN ID')].strip()
                pubmed_id = row[col.index(
                    'PubMed ID (none or blank when not available)')].strip()

                if pubmed_id in ['', None]:
                    continue

                pub_curie = ':'.join(('ZFIN', pub_id))
                pubmed_curie = ':'.join(('PMID', pubmed_id))

                if self.test_mode and (
                        pub_curie not in self.test_ids['pub'] and
                        pubmed_curie not in self.test_ids['pub']):
                    continue

                rtype = self.globaltt['journal article']
                rpm = Reference(graph, pubmed_curie, rtype)
                rpm.addRefToGraph()
                model.addSameIndividual(pub_curie, pubmed_curie)
                ref = Reference(graph, pub_curie, rtype)
                ref.addRefToGraph()
                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

    def _process_targeting_reagents(self, src_key, limit=None):
        """
        This method processes the gene targeting knockdown reagents,
        such as morpholinos, talens, and crisprs.
        We create triples for the reagents and pass the data into a hash map
        for use in the pheno_enviro method.

        Morpholinos work similar to RNAi.
        TALENs are artificial restriction enzymes
            that can be used for genome editing in situ.
        CRISPRs are knockdown reagents, working similar to RNAi
            but at the transcriptional level instead of mRNA level.

        You can read more about TALEN and CRISPR techniques in review
        [Gaj et al]
        http://www.cell.com/trends/biotechnology/abstract/S0167-7799%2813%2900087-5

        TODO add sequences

        Triples created:
        <reagent_id> is a gene_targeting_reagent
        <reagent_id> rdfs:label <reagent_symbol>
        <reagent_id> has type <reagent_so_id>
        <reagent_id> has comment <note>

        <publication_id> is an individual
        <publication_id> mentions <morpholino_id>
        :param src_key: should be one of: morph, talen, crispr
        :param limit:
        :return:
        """
        # src_key passed in
        LOG.info("Processing Gene Targeting Reagents of type %s", src_key,)
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info('from file %s', raw)
        reagent_types = ['morph', 'talen', 'crispr']
        if src_key not in reagent_types:
            LOG.error(
                '%s is not an expected reagent type,\nKnown types are : %s ',
                src_key, reagent_types)
            return

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        kol = [
            'gene_num',
            'gene_so_id',
            'gene_symbol',
            # here we are using common 'reagent_' names for columns in different files
            'reagent_num',
            'reagent_so_id',
            'reagent_symbol',
            'reagent_sequence',
            'publication',
            'note']

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)
                    if row[0][0:9] != 'ZDB-GENE-':  # too messed up (row overflow)
                        LOG.error('Cannot be valid: %s', row)
                        continue

                publication = note = ''
                gene_num = row[col.index('Gene ID')].strip()
                # gene_so_id = row[col.index('Gene SO ID')]
                # gene_symbol = row[col.index('Gene Symbol')]
                reagent_num = row[kol.index('reagent_num')].strip()
                reagent_so_id = row[kol.index('reagent_so_id')]
                reagent_symbol = row[kol.index('reagent_symbol')]
                # reagent_sequence = row[kol.index('reagent_sequence')]
                # if src_key == 'talen':
                #    reagent_sequence2 = row[col.index('TALEN Target Sequence 2')]
                try:
                    publication = row[col.index('Publication(s)')].strip()
                    note = row[col.index('Note')].strip()
                except IndexError:
                    LOG.error(
                        'Row: %i is missing values\nExp:\t%s\nGot:\t%s',
                        reader.line_num, col, row)

                reagent_id = 'ZFIN:' + reagent_num
                gene_id = 'ZFIN:' + gene_num

                self.id_label_map[reagent_id] = reagent_symbol

                if self.test_mode and (
                        reagent_num not in self.test_ids['morpholino'] and
                        gene_num not in self.test_ids['gene']):
                    continue

                geno.addGeneTargetingReagent(
                    reagent_id, reagent_symbol, reagent_so_id, gene_id)
                # The reagent targeted gene is added
                # in the pheno_environment processing function.

                # Add publication
                # note that the publications can be comma-delimited,
                # like: ZDB-PUB-100719-4,ZDB-PUB-130703-22
                if publication != '':
                    pubs = re.split(r',', publication)
                    for pub in pubs:
                        pub = pub.strip()
                        pub_id = 'ZFIN:' + pub
                        ref = Reference(graph, pub_id)
                        ref.addRefToGraph()
                        graph.addTriple(pub_id, self.globaltt['mentions'], reagent_id)

                # Add comment?
                if note != '':
                    model.addComment(reagent_id, note)

                # use the variant hash for reagents to list the affected genes
                if reagent_id not in self.variant_loci_genes:
                    self.variant_loci_genes[reagent_id] = [gene_id]
                else:
                    if gene_id not in self.variant_loci_genes[reagent_id]:
                        self.variant_loci_genes[reagent_id] += [gene_id]

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with Reagent type %s", src_key)

    def _process_pheno_enviro(self, limit=None):
        """
        The pheno_environment.txt (became pheno_environment_fish.txt?)
        file ties experimental conditions
        to an environment ID.
        An environment ID may have one or more associated conditions.
        Condition groups present:
        * chemical, physical, physiological,
        * salinity, temperature, _Generic-control

        First, we build a nice human-readable label
        of all of the components of the environment.
        This is added to our global id-label hash.
        TODO
        Eventually, we will add each component of the environment
        to an environmental object, but needs to be modeled first

        Triples created:
        <env_curie> is an Individual
        <env_curie> rdfs:label <environment_label>
        <env_curie> has type GENO:environment

        :param limit:
        :return:

        """
        src_key = 'enviro'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing environments from file: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        env_hash = {}
        envo = Environment(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                environment_id = row[col.index('Environment ID')].strip()
                # zeco_term_name = row[col.index('ZECO Term Name')]
                zeco_term_id = row[col.index('ZECO Term ID (ZECO:ID)')].strip()
                # chebi_term_name = row[col.index('Chebi Term Name')]
                # chebi_term_id = row[col.index('Chebi Term ID (Chebi:ID)')]
                # zfa_term_name = row[col.index('ZFA Term Name')]
                # fa_term_id = row[col.index('ZFA Term ID (ZFA:ID)')]
                # altered_structure_name = row[col.index(
                #    'Affected Structure Subterm Name')]
                # altered_structure_id = row[col.index(
                #   'Affected Structure Subterm ID (GO-CC:ID)')]
                # ncbi_taxon_name = row[col.index('NCBI Taxon Name')]
                # ncbi_taxon_id = row[col.index('NCBI Taxon ID (NCBI Taxon:ID)')]

                env_curie = ':'.join(('ZFIN', environment_id))
                if self.test_mode and env_curie not in self.test_ids['environment']:
                    continue

                # We can start to build the extrinsic genotype using this file.
                # A single environment can have >1 row in the file,
                # so we build a hash to store all the components of
                # the environment first.
                # Then we can build a label containing all the parts.
                # Using a strategy similar to what is used for genotypes
                # to get the VSLCs and GVCs.

                if env_curie not in env_hash:
                    env_hash[env_curie] = []

                # create environmental components, and add to the hash
                # cleanup the "condition" to remove non-id-friendly chars
                cond_id = zeco_term_id
                cond_id = re.sub(r'\W+', '-', cond_id)

                # TODO  Matt re model
                # description is gone
                # condition is gone
                # condition_group is gone
                # subcond_id = description.strip()
                # subcond_id = re.sub(r'\W+', '-', subcond_id)
                # env_component_id = '-'.join((condition_group.strip(),
                #                             cond_id.strip()))
                # if subcond_id != '':
                #   env_component_id = '-'.join((env_component_id, subcond_id))
                # make them blank nodes
                # env_component_id = self.make_id(env_component_id, '_')  # bnode
                # env_condition = condition.strip()
                # env_component_label = condition_group + '[' + condition + ']'
                # if description != '':
                #    env_component_label += ': ' + description

                # self.id_label_map[env_component_id] = env_component_label
                # env_hash[env_curie] += [env_component_id]

                # if env_curie not in enviro_label_hash:
                #   enviro_label_hash[env_curie] = [env_component_id]
                # else:
                #    enviro_label_hash[env_curie].append(env_component_id)

                # add each component to the environment as a part
                # envo.addEnvironmentalCondition(
                #    env_component_id, env_component_label)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

                # End of loop through pheno_env file

        LOG.info("Building complex environments from components")

        self.environment_hash = env_hash

        # iterate through the env hash to build the full environment label
        for env_id in env_hash:
            environment_labels = []
            env_hash[env_id].sort()
            env_component_list = env_hash[env_id]
            for env_comp_id in env_component_list:
                env_comp_label = self.id_label_map[env_comp_id]
                environment_labels += [env_comp_label]
                envo.addComponentToEnvironment(env_id, env_comp_id)
            environment_labels.sort()
            env_label = 'Environment that includes: ' + '; '.join(environment_labels)
            envo.addEnvironment(env_id, env_label)
            self.id_label_map[env_id] = env_label

        LOG.info("Done with environments")

    def _process_mappings(self, limit=None):
        """
        This function imports linkage mappings of various entities
        to genetic locations in cM or cR.
        Entities include sequence variants, BAC ends, cDNA, ESTs, genes,
        PAC ends, RAPDs, SNPs, SSLPs, and  STSs.
        Status: NEEDS REVIEW
        :param limit:
        :return:

        """
        src_key = 'mappings'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing chromosome mappings from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        taxon_id = self.default_taxon_id
        taxon_label = self.globaltcid[taxon_id]

        # genome_id = geno.makeGenomeID(taxon_id)
        geno.addGenome(taxon_id, taxon_label)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                zfin_num = row[col.index('ZFIN ID')].strip()
                # symbol = row[col.index('Symbol')]
                # so_id = row[col.index('SO_id')]  # !!! Does not exist 2020 Jun!
                panel_symbol = row[col.index('Panel Symbol')].strip()
                chromosome = row[col.index('Chromosome')].strip()
                # location = row[col.index('Location')]
                # metric = row[col.index('Metric')]

                if self.test_mode and zfin_num \
                        not in self.test_ids['gene'] + self.test_ids['allele']:
                    continue

                zfin_curie = ':'.join(('ZFIN', zfin_num))
                if zfin_num[:9] == 'ZDB-GENE-':
                    # assume type and label get added elsewhere
                    model.addClassToGraph(zfin_curie, None)
                    geno.addTaxon(taxon_id, zfin_curie)
                elif zfin_num[:8] == 'ZDB-ALT-':
                    # assume type and label get added elsewhere
                    model.addIndividualToGraph(zfin_curie, None)
                    geno.addTaxon(taxon_id, zfin_curie)
                else:
                    continue
                    # skip any of the others
                # ZFIN don't catalog non-fish things, thankfully
                model.makeLeader(zfin_curie)
                # make the chromosome class
                chr_id = makeChromID(chromosome, taxon_id, 'CHR')
                # chr_label = makeChromLabel(chromosome, taxon_label)
                geno.addChromosomeClass(chromosome, taxon_id, taxon_label)

                pinfo = self._get_mapping_panel_info(panel_symbol)
                panel_label = ' '.join((panel_symbol, pinfo['type'], 'map'))
                if pinfo is not None:
                    # add the panel as a genome build
                    panel_curie = ':'.join(('ZFIN', pinfo['id']))
                    geno.addReferenceGenome(panel_curie, panel_label, taxon_id)

                    if panel_symbol != '':
                        model.addSynonym(panel_curie, panel_symbol)
                    if pinfo['name'] != '':
                        model.addDescription(panel_curie, pinfo['name'])

                    # add the mapping-panel chromosome
                    chr_inst_id = makeChromID(chromosome, panel_curie, 'MONARCH')
                    geno.addChromosomeInstance(
                        chromosome, panel_curie, panel_label, chr_id)
                    # add the feature to the mapping-panel chromosome
                    feat = Feature(graph, zfin_curie, None, None)
                    feat.addSubsequenceOfFeature(chr_inst_id)
                    # TODO add the coordinates see:
                    # https://github.com/JervenBolleman/FALDO/issues/24
                else:
                    LOG.error(
                        "There's a panel (%s) we don't have info for", panel_symbol)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with chromosome mappings")

    def _process_uniprot_ids(self, limit=None):
        """
        This method processes the mappings from ZFIN gene IDs to UniProtKB IDs.

        Triples created:
        <zfin_gene_id> a class
        <zfin_gene_id> rdfs:label gene_symbol

        <uniprot_id> is an Individual
        <uniprot_id> has type <polypeptide>

        <zfin_gene_id> has_gene_product <uniprot_id>
        :param limit:
        :return:

        """
        src_key = 'uniprot'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing UniProt IDs from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                gene_id = row[col.index('ZFIN ID')].strip()
                # gene_so_id = row[col.index('SO ID')]
                gene_symbol = row[col.index('Symbol')]
                uniprot_id = row[col.index('UniProt ID')].strip()

                if self.test_mode and gene_id not in self.test_ids['gene']:
                    continue

                gene_curie = ':'.join(('ZFIN', gene_id))
                prot_curie = ':'.join(('UniProtKB', uniprot_id))

                geno.addGene(gene_curie, gene_symbol)
                # TODO: Abstract to one of the model utilities
                model.addIndividualToGraph(
                    prot_curie, None, self.globaltt['polypeptide']
                )
                graph.addTriple(
                    gene_curie, self.globaltt['has gene product'], prot_curie
                )

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with UniProt IDs")

    def _process_human_orthos(self, limit=None):
        """
        This table provides ortholog mappings between zebrafish and humans.
        ZFIN has their own process of creating orthology mappings,
        that we take in addition to other orthology-calling sources
        (like PANTHER). We ignore the omim ids, and only use the gene_id.

        Triples created:
        <zfin gene id> a class
        <zfin gene id> rdfs:label gene_symbol
        <zfin gene id> dc:description gene_name

        <human gene id> a class
        <human gene id> rdfs:label gene_symbol
        <human gene id> dc:description gene_name
        <human gene id> equivalent class <omim id>

        <zfin gene id> orthology association <human gene id>
        :param limit:
        :return:

        """
        src_key = 'human_orthos'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing human orthos from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        # model = Model(graph)  # unused
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                zfin_id = row[col.index('ZFIN ID')].strip()
                zfin_symbol = row[col.index('ZFIN Symbol')]
                zfin_name = row[col.index('ZFIN Name')]
                human_symbol = row[col.index('Human Symbol')]
                human_name = row[col.index('Human Name')]
                # omim_id = row[col.index('OMIM ID')]
                gene_id = row[col.index('Gene ID')].strip()
                # hgnc_id = row[col.index('HGNC ID')]
                evidence_code = row[col.index('Evidence')]
                pub_id = row[col.index('Pub ID')]

                if self.test_mode and zfin_id not in self.test_ids['gene']:
                    continue

                # Add the zebrafish gene.
                zfin_id = 'ZFIN:' + zfin_id
                geno.addGene(zfin_id, zfin_symbol, None, zfin_name)

                # Add the human gene.
                ncbigene_curie = ':'.join(('NCBIGene', gene_id))
                geno.addGene(ncbigene_curie, human_symbol, None, human_name)

                # make the association
                assoc = OrthologyAssoc(graph, self.name, zfin_id, ncbigene_curie)
                # we don't know anything about the orthology type,
                # so we just use the default

                if pub_id[:8] == 'ZDB-PUB-':
                    assoc.add_source('ZFIN:' + pub_id)

                eco_id = self.get_orthology_evidence_code(evidence_code)
                if eco_id is not None:
                    assoc.add_evidence(eco_id)

                assoc.add_association_to_graph()

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with human orthos")

    def _process_gene_coordinates(self, limit=None):

        src_key = 'gene_coordinates'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing human orthos from file: %s", raw)
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        model = Model(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row = next(reader)
            if row[0] != '##gff-version 3':
                LOG.warning('expected "##gff-version 3"; got %s', row)
            row = next(reader)  # blank line expected
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)
                chrom = row[col.index('Chromosome')]
                # source = row[col.index('Source')]
                # ftype = row[col.index('Type')]
                start = row[col.index('Start')]
                end = row[col.index('End')]
                # score = row[col.index('Score')]
                strand = row[col.index('Strand')]
                # phase = row[col.index('Phase')]
                attributes = row[col.index('Attributes')]

                gene_id = None
                if attributes == '':
                    continue

                attribute_dict = dict(
                    item.split("=")
                    for item in re.sub(r'"', '', attributes).split(";"))
                gene_id = attribute_dict.get('gene_id')

                if self.test_mode and gene_id not in self.test_ids['gene']:
                    continue
                gene_curie = ':'.join(('ZFIN', gene_id))

                # make chrom
                chrom_id = makeChromID(chrom, self.globaltt['Danio rerio'], 'CHR')
                # assume it gets added elsewhere
                model.addClassToGraph(chrom_id, None)
                # FIXME - remove this hardcoding,
                # difficult. zfin is not forwarding the build metadata
                build_label = 'danRer10'
                build_curie = ':'.join(('UCSC', build_label))
                chrom_in_build = makeChromID(chrom, build_curie, 'MONARCH')
                geno.addChromosomeInstance(
                    chrom, build_curie, build_label, chrom_id)
                feat = Feature(graph, gene_curie, None, None)
                feat.addFeatureStartLocation(start, chrom_in_build, strand)
                feat.addFeatureEndLocation(end, chrom_in_build, strand)
                feat.addFeatureToGraph(True, None, True)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with gene coordinates")

    def process_fish_disease_models(self, limit=None):
        '''

        '''
        src_key = 'fish_disease_models'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing fish models from file %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        model = Model(graph)
        col = self.files[src_key]['columns']
        collen = len(col)
        fish_taxon = self.globaltt['Danio rerio']
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row = next(reader)

            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                if re.match(r'^(\s|#|$)', ''.join(row)):
                    continue  # skip header

                # ZDB-FISH-150901-9014
                # ZDB-EXP-041102-1
                # is_a_model
                # DOID:5603
                # acute T cell leukemia
                # ZDB-PUB-110523-12
                # 21552289
                # TAS, ECO:0000033

                fish = row[col.index('Fish ZDB ID')].strip()
                environment = row[col.index('Environment ZDB ID')]
                # rel = row[col.index('is_a_model')]
                disease_id = row[col.index('DO Term ID')]
                disease_label = row[col.index('DO Term Name')]
                zfin_pub_id = row[col.index('Publication ZDB ID')]
                pubmed_id = row[col.index('PubMed ID')]
                # evidence_code = row[col.index('Evidence Code')]

                # if no fish is listed, then don't add it.
                if fish == '':
                    continue

                if self.test_mode and (
                        fish not in self.test_ids['fish'] and
                        disease_id not in self.test_ids['disease']):
                    continue

                # make effective genotype = fish+environment
                fish_id = 'ZFIN:' + fish
                fish_label = self.id_label_map.get(fish_id)
                if fish_label is None:
                    fish_label = fish_id
                environment_id = 'ZFIN:' + environment
                environment_label = self.id_label_map.get(environment_id)
                if environment_label is None:
                    environment_label = environment_id
                geno.make_experimental_model_with_genotype(
                    fish_id, fish_label, fish_taxon, 'zebrafish')

                assoc = Assoc(graph, self.name)
                assoc.set_subject(fish_id)
                assoc.set_object(disease_id)
                assoc.set_relationship(self.globaltt['is model of'])
                desc = ' '.join((
                    'A fish with genotype', fish_label,
                    'is a model for disease', disease_label,
                    'under the condition of', environment_label))
                assoc.set_description(desc)

                if zfin_pub_id != '':
                    assoc.add_source('ZFIN:' + zfin_pub_id)

                # make the pubmed id, if it exists
                if pubmed_id != '':
                    pubmed_id = 'PMID:' + pubmed_id
                    model.addSameIndividual('ZFIN:' + zfin_pub_id, pubmed_id)
                    model.makeLeader(pubmed_id)
                assoc.add_association_to_graph()

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

    def _map_octuple_to_phenotype(
            self,
            subterm1_id, post_composed_relationship_id_1, superterm1_id, quality_id,
            subterm2_id, post_composed_relationship_id_2, superterm2_id, modifier):
        """
        This will take the 8-part EQ-style annotation
        used in zp-mapping.txt and return the ZP id.

        :param subterm1_id,
        :param post_composed_relationship_id_1,
        :param superterm1_id,
        :param quality_id,
        :param subterm2_id,
        :param post_composed_relationship_id_2,
        :param superterm2_id,
        :param modifier
        """

        zp_id = None

        # zfin uses free-text modifiers,
        # but we need to convert them to proper PATO classes for the mapping
        mod_id = self.resolve(modifier, False)

        if modifier == mod_id:
            LOG.warning("no mapping for pato modifier %s", modifier)

        key = self._make_zpkey(
            [subterm1_id, post_composed_relationship_id_1, superterm1_id, quality_id,
             subterm2_id, post_composed_relationship_id_2, superterm2_id, mod_id])

        mapping = self.zp_map.get(key)

        if mapping is None:
            if modifier == 'normal':
                pass
            else:
                LOG.warning(
                    "Couldn't map ZP id to %s with modifier %s", "_"
                    .join((
                        subterm1_id,
                        post_composed_relationship_id_1,
                        superterm1_id,
                        quality_id,
                        subterm2_id,
                        post_composed_relationship_id_2,
                        superterm2_id,
                        mod_id)), modifier)
        else:
            zp_id = mapping['zp_id']

        return zp_id

    def _load_zp_mappings(self, src_key):
        """
        Given a file that defines the mapping between
        ZFIN-specific EQ definitions and the automatically derived ZP ids,
        create a mapping here.
        :return:
        dict zp_map with
        """
        #  'zpmap'
        zp_file = '/'.join((self.rawdir, self.files[src_key]['file']))
        zp_map = {}
        col = self.files[src_key]['columns']
        kol = [
            'zp_id',
            'subterm1_id',
            'post_composed_relationship_id_1',
            'superterm1_id',
            'quality_id',
            'subterm2_id'
            'post_composed_relationship_id_2',
            'superterm2_id',
            'modifier'
        ]
        # id_map_zfin.tsv only contains data for abnormal phenotypes
        modifier = self.globaltt['abnormal']
        LOG.info("Loading ZP-to-EQ mappings from %s", zp_file)
        with open(zp_file, 'r', encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row = next(reader)
            if not self.check_fileheader(col, row):
                pass
            for row in reader:
                zp_id = row[col.index('iri')]
                # the 'id' column is a hyphen separated list of 7 items.
                val = row[col.index('id')].split('-')
                val.append(modifier)
                key = self._make_zpkey(val)
                val.insert(0, zp_id)
                zp_map[key] = dict(zip(kol, val))

        LOG.info("Loaded %s zp terms", zp_map.__len__())

        return zp_map

    def _make_zpkey(self, args_for_key, replace_empty_args_with_zeros=True):
        # empty args -> zero is the convention in Nico's file and elsewhere
        if replace_empty_args_with_zeros:
            args_for_key = [element or '0' for element in args_for_key]

        key = self.make_id('_'.join((args_for_key)))
        return key

    @staticmethod
    def _get_other_allele_by_zygosity(allele_id, zygosity):
        """
        A helper function to switch on the zygosity,
        and return the appropriate allele id, or symbol.
        :param allele_id:
        :param zygosity:
        :return:
        """
        other_allele = None
        if zygosity == 'homozygous':
            other_allele = allele_id
        elif zygosity == 'hemizygous':
            other_allele = '0'
        elif zygosity == 'unknown':  # we'll use this as a convention
            other_allele = '?'
        elif zygosity == 'complex':  # transgenics
            other_allele = '0'
        elif zygosity == 'heterozygous':
            # passing on hets until a different fxn
            pass
        else:
            LOG.warning("Unconfigured zygosity: %s", zygosity)

        return other_allele

    @staticmethod
    def _get_mapping_panel_info(panel):
        panel_hash = {
            'HS': {
                'id': 'ZDB-REFCROSS-000320-1',
                'name': 'Heat Shock',
                'type': 'meiotic',
                'num_meioses': 42},
            'GAT': {
                'id': 'ZDB-REFCROSS-990308-7',
                'name': 'Gates et al',
                'type': 'meiotic',
                'num_meioses': 96},
            'LN54': {
                'id': 'ZDB-REFCROSS-990426-6',
                'name': 'Loeb/NIH/5000/4000',
                'dose': '4000 rads',
                'type': 'Radiation Hybrid'},
            'MGH': {
                'id': 'ZDB-REFCROSS-980521-11',
                'name': 'Boston MGH Cross',
                'type': 'meiotic',
                'num_meioses': 40},
            'MOP': {
                'id': 'ZDB-REFCROSS-980526-5',
                'name': 'Mother of Pearl',
                'type': 'meiotic',
                'num_meioses': 96},
            'T51': {
                'id': 'ZDB-REFCROSS-990707-1',
                'name': 'Goodfellow T51',
                'dose': '3000 rads',
                'type': 'Radiation Hybrid'},
        }

        return panel_hash.get(panel)

    @staticmethod
    def _make_variant_locus_id(gene_id, allele_id):
        """
        A convenience method to uniformly create variant loci.
        If we want to materialize these in the monarch space,
        then we wrap with the self.make_id function.
        :param gene_id:
        :param allele_id:
        :return:

        """

        varloci = '-'.join((gene_id, allele_id))
        varloci = Source.make_id(re.sub(r'(ZFIN)?:', '', varloci), '_')

        return varloci

    def _make_effective_genotype_id(self, intrinsic_id, extrinsic_id):
        effective_genotype_id = self.make_id(
            '-'.join((intrinsic_id, extrinsic_id)))

        return effective_genotype_id

    def get_orthology_sources_from_zebrafishmine(self):
        """
        Fetch the zfin gene to other species orthology annotations,
        together with the evidence for the assertion.
        Write the file locally to be read in a separate function.
        :return:

        """

        # For further documentation you can visit:
        #     http://www.intermine.org/wiki/PythonClient

        # The following two lines will be needed in every python script:
        service = Service("http://zebrafishmine.org/service")

        # Get a new query on the class (table) you will be querying:
        query = service.new_query("Gene")

        # The view specifies the output columns
        query.add_view(
            "primaryIdentifier", "symbol", "homologues.homologue.symbol",
            "homologues.evidence.evidenceCode.abbreviation",
            "homologues.evidence.publications.primaryIdentifier",
            "homologues.evidence.publications.pubMedId",
            "homologues.crossReferences.identifier",
            # only "orthologue" is used
            # "homologues.crossReferences.linkType",
            # only needed if >1 source
            # "homologues.crossReferences.source.name"
        )

        # This query's custom sort order is specified below:
        query.add_sort_order("Gene.name", "ASC")

        # You can edit the constraint values below
        query.add_constraint(
            "homologues.dataSets.name", "=",
            "ZFIN Curated Human, Mouse, Fly, Yeast Orthologue Data Set",
            code="A")
        query.add_constraint(
            "homologues.homologue.organism.name", "=",
            "Homo sapiens", code="B")
        query.add_constraint(
            "homologues.crossReferences.source.name", "=",
            "Gene", code="D")  # NCBIGene
        query.add_constraint("symbol", "=", "*", code="C")

        # Uncomment and edit the code below to specify your own custom logic:
        # query.set_logic("C and A and B and D and D")

        src_key = 'zmine_ortho_evidence'
        outfile = '/'.join((self.rawdir, self.files[src_key]['file']))

        with open(outfile, 'w', encoding="utf-8", newline='\n') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            for row in query.rows():
                stuff = [
                    row["primaryIdentifier"],
                    row["symbol"],
                    row["homologues.homologue.symbol"],
                    row["homologues.crossReferences.identifier"],
                    row["homologues.evidence.evidenceCode.abbreviation"],
                    row["homologues.evidence.publications.primaryIdentifier"],
                    row["homologues.evidence.publications.pubMedId"],
                    # row["homologues.crossReferences.linkType"],
                    # row["homologues.crossReferences.source.name"]
                ]
                filewriter.writerow(stuff)

    def process_orthology_evidence(self, limit):
        # this ia an odd  one see:      get_orthology_sources_from_zebrafishmine()

        src_key = 'zmine_ortho_evidence'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Processing orthology evidence from file: %s", raw)

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        col = self.files[src_key]['columns']
        collen = len(col)
        with open(raw, 'r', encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in reader:
                if len(row) != collen:
                    LOG.warning('Row: %i has unexpected format', reader.line_num)

                zfin_gene_num = row[col.index('zfin_gene_num')].strip()
                # zfin_gene_symbol = row[col.index('zfin_gene_symbol')]
                # ortholog_gene_symbol = row[col.index('ortholog_gene_symbol')]
                ortholog_ncbigene_num = row[col.index('ortholog_ncbigene_num')]
                evidence_code = row[col.index('evidence_code')]
                zfin_pub_num = row[col.index('zfin_pub_num')]
                pubmed_num = row[col.index('pubmed_num')]

                if self.test_mode and (zfin_gene_num not in self.test_ids['gene']):
                    continue

                zfin_gene_id = 'ZFIN:' + str(zfin_gene_num)
                ortho_gene_id = 'NCBIGene:' + str(ortholog_ncbigene_num)
                zfin_pub_id = 'ZFIN:' + zfin_pub_num
                pubmed_id = 'PMID:' + str(pubmed_num)
                if zfin_gene_num != '' and ortholog_ncbigene_num != '':
                    assoc = OrthologyAssoc(
                        graph, self.name, zfin_gene_id, ortho_gene_id)
                    if zfin_pub_num != '':
                        ref = Reference(graph, zfin_pub_id)
                        ref.addRefToGraph()
                        assoc.add_source(zfin_pub_id)
                    if pubmed_num != '':
                        ref = Reference(
                            graph, pubmed_id, self.globaltt['journal article'])
                        ref.addRefToGraph()
                        assoc.add_source(pubmed_id)
                    if evidence_code != '':
                        eco_id = self.get_orthology_evidence_code(evidence_code)
                        assoc.add_evidence(eco_id)

                    assoc.add_association_to_graph()
                # FIXME need to update with proper provenance model
                # so the papers get attached with the relevant eco code

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

    @staticmethod
    def get_orthology_evidence_code(abbrev):
        '''
            move to localtt & globltt
        '''
        # AA    Amino acid sequence comparison.
        # CE    Coincident expression.
        # CL    Conserved genome location (synteny).
        # FC    Functional complementation.
        # FH    Formation of functional heteropolymers.
        # IX    Immunological cross-reaction.
        # NS    Not specified.
        # NT    Nucleotide sequence comparison.
        # SI    Similar response to inhibitors.
        # SL    Similar subcellular location.
        # SS    Similar substrate specificity.
        # SU    Similar subunit structure.
        # XH    Cross-hybridization to same molecular probe.
        # PT    Phylogenetic Tree.
        # OT    Other

        eco_abbrev_map = {
            'AA': 'ECO:0000031',  # BLAST protein sequence similarity evidence
            'CE': 'ECO:0000008',  # expression evidence
            'CL': 'ECO:0000044',  # sequence similarity FIXME
            'FC': 'ECO:0000012',  # functional complementation
            # functional complementation in a heterologous system
            'FH': 'ECO:0000064',
            'IX': 'ECO:0000040',  # immunological assay evidence
            'NS': None,
            'NT': 'ECO:0000032',  # nucleotide blast
            'SI': 'ECO:0000094',  # biological assay evidence FIXME
            'SL': 'ECO:0000122',  # protein localization evidence FIXME
            'SS': 'ECO:0000024',  # protein binding evidence  FIXME
            'SU': 'ECO:0000027',  # structural similarity evidence
            'XH': 'ECO:0000002',  # direct assay evidence  FIXME
            'PT': 'ECO:0000080',  # phylogenetic evidence
            'OT': None,
        }

        if abbrev not in eco_abbrev_map:
            LOG.warning("Evidence code for orthology (%s) not mapped", str(abbrev))

        return eco_abbrev_map.get(abbrev)

    @staticmethod
    def make_targeted_gene_id(geneid, reagentid):

        targeted_gene_id = '-'.join((geneid, reagentid))
        # these are not zfin resolvable, so make BNodes
        targeted_gene_id = re.sub(r'(ZFIN)?:', '', targeted_gene_id)
        targeted_gene_id = Source.make_id(targeted_gene_id, '_')
        return targeted_gene_id

    # #### SOME OLD (COBOL?) CODE #####

# # LINK THE MORPHOLINO TO THE GENES THAT IT AFFECTS
# AG = SELF.VARIANT_LOCI_GENES.GET(MORPH_ID)
# # LOGGER.INFO("%S AFFECTED GENES %S", MORPH_ID, PP.PFORMAT(AG))
# LIST_OF_TARGETED_GENES = []
# IF AG IS NONE:
#     LOGGER.WARN("NO AFFECTED GENES FOR %S", MORPH_ID)
# ELSE:
#     # CREATE VARIANT GENE(S) THAT HAVE BEEN TARGETED BY THE REAGENT
#     FOR GID IN AG:
#         IF GID NOT IN SELF.ID_LABEL_MAP:
#             # SHOULD NOT HAPPEN, EXCEPT MAYBE IN TESTING
#             LOGGER.ERROR("%S NOT IN ID-LABEL-HASH", GID)
#             GLABEL = GID
#         ELSE:
#             GLABEL = SELF.ID_LABEL_MAP[GID]
#
#         TARGETED_GENE_ID = '-'.JOIN((GID, APPLIED_MORPH_ID))
#         # THESE ARE NOT ZFIN RESOLVABLE, SO MAKE BNODES
#         TARGETED_GENE_ID = RE.SUB('(ZFIN)?:', '', TARGETED_GENE_ID)
#         TARGETED_GENE_ID = '_' + TARGETED_GENE_ID
#         IF SELF.NOBNODES:
#             TARGETED_GENE_ID = ':' + TARGETED_GENE_ID
#         TARGETED_GENE_LABEL = GLABEL + '<' + APPLIED_MORPH_LABEL + '>'
#
#         GENO.ADD(MORPH_ID, GID, TARGETED_GENE_ID, TARGETED_GENE_LABEL)
#         SELF.ID_LABEL_MAP[TARGETED_GENE_ID] = TARGETED_GENE_LABEL
#         LIST_OF_TARGETED_GENES += [TARGETED_GENE_ID]
