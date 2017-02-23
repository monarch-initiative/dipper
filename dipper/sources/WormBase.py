import csv
import re
import logging
import gzip
import io
from ftplib import FTP

from dipper.sources.Source import Source
from dipper.models.Genotype import Genotype
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.GenomicFeature import makeChromID, Feature
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper.models.assoc.InteractionAssoc import InteractionAssoc

logger = logging.getLogger(__name__)


class WormBase(Source):
    """
    This is the parser for the
    [C. elegans Model Organism Database (WormBase)](http://www.wormbase.org),
    from which we process genotype and phenotype data for laboratory worms
    (C.elegans and other nematodes).

    We generate the wormbase graph to include the following information:
    * genes
    * sequence alterations (includes SNPs/del/ins/indel and
        large chromosomal rearrangements)
    * RNAi as expression-affecting reagents
    * genotypes, and their components
    * strains
    * publications (and their mapping to PMIDs, if available)
    * allele-to-phenotype associations (including variants by RNAi)
    * genetic positional information for genes and sequence alterations

    Genotypes leverage the GENO genotype model and includes both
    intrinsic and extrinsic genotypes.  Where necessary, we create anonymous
    nodes of the genotype partonomy (i.e. for variant single locus complements,
    genomic variation complements, variant loci, extrinsic genotypes, and
    extrinsic genotype parts).

    TODO:  get people and gene expression
    """
    wbrel = 'ftp://ftp.wormbase.org/pub/wormbase/releases'
    wbdev = wbrel + '/current-development-release'
    wbprod = wbrel + '/current-production-release'
    species = '/species/c_elegans/PRJNA13758'
    files = {
        'gene_ids': {
            'file': 'c_elegans.PRJNA13758.geneIDs.txt.gz',
            'url': wbprod + species +
            '/annotation/c_elegans.PRJNA13758.WSNUMBER.geneIDs.txt.gz'},
        # 'gene_desc': { # TEC: missing as of 2016 Mar 03
        #    'file': 'c_elegans.PRJNA13758.functional_descriptions.txt.gz',
        #    'url': wbdev + species +
        #     '/annotation/c_elegans.PRJNA13758.WSNUMBER.functional_descriptions.txt.gz'},
        'allele_pheno': {
            'file': 'phenotype_association.wb',
            'url': wbprod + '/ONTOLOGY/phenotype_association.WSNUMBER.wb'},
        'rnai_pheno': {
            'file': 'rnai_phenotypes.wb',
            'url': wbprod + '/ONTOLOGY/rnai_phenotypes.WSNUMBER.wb'},
        'pub_xrefs': {
            'file': 'pub_xrefs.txt',
            'url': 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref'},
        'feature_loc': {
            'file': 'c_elegans.PRJNA13758.annotations.gff3.gz',
            'url': wbprod + species +
            '/c_elegans.PRJNA13758.WSNUMBER.annotations.gff3.gz'},
        'disease_assoc': {
            'file': 'disease_association.wb',
            'url': 'ftp://ftp.sanger.ac.uk/pub/wormbase/releases/WSNUMBER/ONTOLOGY/disease_association.WSNUMBER.wb'},
        # 'genes_during_development': {
        #   'file': 'development_association.wb',
        #   'url wbdev+'/ONTOLOGY/development_association.WS249.wb'},
        # 'genes_in_anatomy': {
        #   'file': 'anatomy_association.wb',
        #   'url': wbdev+'/ONTOLOGY/anatomy_association.WS249.wb'},
        # 'gene_interaction': {
        #   'file': 'c_elegans.PRJNA13758.gene_interactions.txt.gz',
        #   'url': wbdev + species +
        #   '/annotation/c_elegans.PRJNA13758.WSNUMBER.gene_interactions.txt.gz'},
        # 'orthologs': {
        #   'file': 'c_elegans.PRJNA13758.orthologs.txt.gz',
        #   'url': wbdev + species +
        #   '/annotation/c_elegans.PRJNA13758.WS249.orthologs.txt.gz'},
        'xrefs': {
            'file': 'c_elegans.PRJNA13758.xrefs.txt.gz',
            'url': wbprod + species +
            '/c_elegans.PRJNA13758.WSNUMBER.xrefs.txt.gz'},
        # 'letter': { # no longer exists 2016-11-18
        #    'file': 'letter.WSNUMBER',
        #    'url': wbprod + '/letter.WSNUMBER'},
        'checksums': {
            'file': 'CHECKSUMS',
            'url':  wbprod + '/CHECKSUMS'}
    }

    test_ids = {
        'gene': [
            'WBGene00001414', 'WBGene00004967', 'WBGene00003916',
            'WBGene00004397', 'WBGene00001531'],
        'allele': [
            'WBVar00087800', 'WBVar00087742', 'WBVar00144481', 'WBVar00248869',
            'WBVar00250630'],
        'strain': ['BA794', 'RK1', 'HE1006'],
        'pub': []  # FIXME
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'wormbase')

        # update the dataset object with details about this resource
        # NO LICENSE for this resource
        self.dataset = Dataset(
            'wormbase', 'WormBase', 'http://www.wormbase.org', None, None,
            'http://www.wormbase.org/about/policies#012')

        self.version_num = None
        return

    def fetch(self, is_dl_forced=False):

        # figure out the version number by probing the "current_release",
        # then edit the file dict accordingly
        # connect to wormbase ftp
        current_dev_release_dir = \
            'pub/wormbase/releases/current-production-release'
        ftp = FTP('ftp.wormbase.org')
        ftp.login()
        ftp.cwd(current_dev_release_dir)
        # the current release dir is a redirect to a versioned release.
        # pull that from the pwd.
        pwd = ftp.pwd()
        ftp.quit()
        wsver = re.search(r'releases\/(WS\d+)', pwd)
        if wsver is None or len(wsver.groups()) < 1:
            logger.error(
                "Couldn't figure out version number from FTP site.  Exiting.")
            exit(1)
        else:

            self.update_wsnum_in_files(wsver.group(1))

        self.dataset.set_version_by_num(self.version_num)
        # fetch all the files
        self.get_files(is_dl_forced)
        return

    def update_wsnum_in_files(self, vernum):
        """
        With the given version number ```vernum```,
        update the source's version number, and replace in the file hashmap.
        the version number is in the CHECKSUMS file.
        :param vernum:
        :return:

        """
        self.version_num = vernum
        # replace the WSNUMBER in the url paths with the real WS###
        for f in self.files:
            url = self.files[f].get('url')
            url = re.sub(r'WSNUMBER', self.version_num, url)
            self.files[f]['url'] = url
            logger.debug(
                "Replacing WSNUMBER in %s with %s", f, self.version_num)

        # also the letter file - keep this so we know the version number
        # self.files['checksums']['file'] = re.sub(
        #    r'WSNUMBER', self.version_num, self.files['checksums']['file'])
        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows of each file", limit)

        if self.version_num is None:
            logger.info("Figuring out version num for files")
            # probe the raw directory for the WSnumber incthe "CHECKSUMS" file.
            # 20f7d39c73012c9cfc8444a657af2b80  acedb/md5sum.WS255

            checksums = open(self.rawdir + '/CHECKSUMS', 'r')
            checksum = checksums.readline()
            vernum = re.search(r'\.(WS\d+)', checksum)
            self.update_wsnum_in_files(vernum.group(1))
            checksums.close()

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        # to hold any label for a given id
        self.id_label_map = {}
        # to hold the mappings between genotype and background
        self.genotype_backgrounds = {}
        self.extrinsic_id_to_enviro_id_hash = {}
        # to hold the genes variant due to a seq alt
        self.variant_loci_genes = {}
        # to hold the parts of an environment
        self.environment_hash = {}
        self.wildtype_genotypes = []
        # stores the rnai_reagent to gene targets
        self.rnai_gene_map = {}

        self.process_gene_ids(limit)
        # self.process_gene_desc(limit)   #TEC imput file is mia 2016-Mar-03
        self.process_allele_phenotype(limit)
        self.process_rnai_phenotypes(limit)
        self.process_pub_xrefs(limit)
        self.process_feature_loc(limit)
        self.process_disease_association(limit)
        # TODO add this when when complete
        # self.process_gene_interaction(limit)

        logger.info("Finished parsing.")
        return

    def process_gene_ids(self, limit):
        raw = '/'.join((self.rawdir, self.files['gene_ids']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        model = Model(g)
        logger.info("Processing Gene IDs")
        line_counter = 0
        geno = Genotype(g)
        with gzip.open(raw, 'rb') as csvfile:
            filereader = csv.reader(
                io.TextIOWrapper(csvfile, newline=""), delimiter=',',
                quotechar='\"')
            for row in filereader:
                line_counter += 1
                (taxon_num, gene_num, gene_symbol, gene_synonym, live) = row
                # 6239,WBGene00000001,aap-1,Y110A7A.10,Live

                if self.testMode and gene_num not in self.test_ids['gene']:
                    continue

                taxon_id = 'NCBITaxon:'+taxon_num
                gene_id = 'WormBase:'+gene_num
                if gene_symbol == '':
                    gene_symbol = gene_synonym
                if gene_symbol == '':
                    gene_symbol = None
                model.addClassToGraph(
                    gene_id, gene_symbol, Genotype.genoparts['gene'])
                if live == 'Dead':
                    model.addDeprecatedClass(gene_id)
                geno.addTaxon(taxon_id, gene_id)
                if gene_synonym != '' and gene_synonym is not None:
                    model.addSynonym(gene_id, gene_synonym)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def process_gene_desc(self, limit):
        raw = '/'.join((self.rawdir, self.files['gene_desc']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("Processing Gene descriptions")
        line_counter = 0
        # geno = Genotype(g)  # TODO unused
        with gzip.open(raw, 'rb') as csvfile:
            filereader = csv.reader(
                io.TextIOWrapper(csvfile, newline=""), delimiter='\t',
                quotechar='\"')
            for row in filereader:
                if re.match(r'\#', ''.join(row)):
                    continue
                line_counter += 1
                if line_counter == 1:
                    continue
                (gene_num, public_name, molecular_name, concise_description,
                 provisional_description, detailed_description,
                 automated_description, gene_class_description) = row

                if self.testMode and gene_num not in self.test_ids['gene']:
                    continue

                gene_id = 'WormBase:'+gene_num

                if concise_description != 'none available':
                    model.addDefinition(gene_id, concise_description)

                # remove the description if it's identical to the concise
                descs = {
                    'provisional': provisional_description,
                    'automated': automated_description,
                    'detailed': detailed_description,
                    'gene class': gene_class_description
                }
                for d in descs:
                    text = descs.get(d)
                    if text == concise_description \
                            or re.match(r'none', text) or text == '':
                        pass  # don't use it
                    else:
                        text = ' '.join((text, '['+d+']'))
                        descs[d] = text
                        model.addDescription(gene_id, text)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def process_allele_phenotype(self, limit=None):
        """
        This file compactly lists variant to phenotype associations,
        such that in a single row, there may be >1 variant listed
        per phenotype and paper.  This indicates that each variant is
        individually assocated with the given phenotype,
        as listed in 1+ papers.
        (Not that the combination of variants is producing the phenotype.)
        :param limit:
        :return:

        """

        raw = '/'.join((self.rawdir, self.files['allele_pheno']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        logger.info("Processing Allele phenotype associations")
        line_counter = 0
        geno = Genotype(g)
        with open(raw, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                if re.match(r'!', ''.join(row)):  # header
                    continue
                line_counter += 1
                (db, gene_num, gene_symbol, is_not, phenotype_id, ref,
                 eco_symbol, with_or_from, aspect, gene_name, gene_synonym,
                 gene_class, taxon, date, assigned_by, blank, blank2) = row

                if self.testMode and gene_num not in self.test_ids['gene']:
                    continue

                # TODO add NOT phenotypes
                if is_not == 'NOT':
                    continue

                eco_id = None
                if eco_symbol == 'IMP':
                    eco_id = 'ECO:0000015'
                elif eco_symbol.strip() != '':
                    logger.warning(
                        "Encountered an ECO code we don't have: %s",
                        eco_symbol)

                # according to the GOA spec, persons are not allowed to be
                # in the reference column, therefore they the variant and
                # persons are swapped between the reference and with column.
                # we unswitch them here.
                temp_var = temp_ref = None
                if re.search(r'WBVar|WBRNAi', ref):
                    temp_var = ref
                    # move the paper from the with column into the ref
                if re.search(r'WBPerson', with_or_from):
                    temp_ref = with_or_from
                if temp_var is not None or temp_ref is not None:
                    with_or_from = temp_var
                    ref = temp_ref

                allele_list = re.split(r'\|', with_or_from)
                if len(allele_list) == 0:
                    logger.error(
                        "Missing alleles from phenotype assoc at line %d",
                        line_counter)
                    continue
                else:
                    for a in allele_list:
                        allele_num = re.sub(r'WB:', '', a.strip())
                        allele_id = 'WormBase:'+allele_num
                        gene_id = 'WormBase:'+gene_num

                        if re.search(r'WBRNAi', allele_id):
                            # make the reagent-targeted gene,
                            # & annotate that instead of the RNAi item directly
                            rnai_num = re.sub(r'WormBase:', '', allele_id)
                            rnai_id = allele_id
                            rtg_id = self.make_reagent_targeted_gene_id(
                                gene_num, rnai_num)
                            geno.addReagentTargetedGene(
                                rnai_id, 'WormBase:'+gene_num, rtg_id)
                            geno.addGeneTargetingReagent(
                                rnai_id, None, geno.genoparts['RNAi_reagent'],
                                gene_id)
                            allele_id = rtg_id
                        elif re.search(r'WBVar', allele_id):
                            # this may become deprecated by using wormmine
                            # make the allele to gene relationship
                            # the WBVars are really sequence alterations

                            # the public name will come from elsewhere
                            geno.addSequenceAlteration(allele_id, None)
                            vl_id = '_:'+'-'.join((gene_num, allele_num))
                            geno.addSequenceAlterationToVariantLocus(
                                allele_id, vl_id)
                            geno.addAlleleOfGene(vl_id, gene_id)
                        else:
                            logger.warning(
                                "Some kind of allele I don't recognize: %s",
                                allele_num)
                            continue
                        assoc = G2PAssoc(g, self.name, allele_id, phenotype_id)

                        if eco_id is not None:
                            assoc.add_evidence(eco_id)

                        if ref is not None and ref != '':
                            ref = re.sub(r'(WB:|WB_REF:)', 'WormBase:', ref)
                            reference = Reference(g, ref)
                            if re.search(r'Person', ref):
                                reference.setType(reference.ref_types['person'])
                                # also add
                                # inferred from background scientific knowledge
                                assoc.add_evidence('ECO:0000001')
                            reference.addRefToGraph()
                            assoc.add_source(ref)

                        assoc.add_association_to_graph()

                        # finish looping through all alleles

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def process_rnai_phenotypes(self, limit=None):

        raw = '/'.join((self.rawdir, self.files['rnai_pheno']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)

        logger.info("Processing RNAi phenotype associations")
        line_counter = 0
        geno = Genotype(g)
        with open(raw, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_num, gene_alt_symbol, phenotype_label, phenotype_id,
                 rnai_and_refs) = row
# WBGene00001908	F17E9.9	locomotion variant	WBPhenotype:0000643	WBRNAi00025129|WBPaper00006395 WBRNAi00025631|WBPaper00006395
# WBGene00001908	F17E9.9	avoids bacterial lawn	WBPhenotype:0000402	WBRNAi00095640|WBPaper00040984
# WBGene00001908	F17E9.9	RAB-11 recycling endosome localization variant	WBPhenotype:0002107	WBRNAi00090830|WBPaper00041129

                if self.testMode and gene_num not in self.test_ids['gene']:
                    continue

                gene_id = 'WormBase:'+gene_num
                # refs = list()  # TODO unused

                # the rnai_and_refs has this so that
                # WBRNAi00008687|WBPaper00005654 WBRNAi00025197|WBPaper00006395 WBRNAi00045381|WBPaper00025054
                # space delimited between RNAi sets;
                # then each RNAi should have a paper

                rnai_sets = re.split(r' ', rnai_and_refs)

                for s in rnai_sets:

                    # get the rnai_id
                    (rnai_num, ref_num) = re.split(r'\|', s)
                    if len(re.split(r'\|', s)) > 2:
                        logger.warning(
                            "There's an unexpected number of items in %s", s)
                    if rnai_num not in self.rnai_gene_map:
                        self.rnai_gene_map[rnai_num] = set()

                    # to use for looking up later
                    self.rnai_gene_map[rnai_num].add(gene_num)

                    rnai_id = 'WormBase:'+rnai_num
                    geno.addGeneTargetingReagent(
                        rnai_id, None, geno.genoparts['RNAi_reagent'], gene_id)

                    # make the "allele" of the gene
                    # that is targeted by the reagent
                    allele_id = self.make_reagent_targeted_gene_id(
                        gene_num, rnai_num)
                    allele_label = gene_alt_symbol+'<'+rnai_num+'>'
                    geno.addReagentTargetedGene(
                        rnai_id, gene_id, allele_id, allele_label)

                    assoc = G2PAssoc(g, self.name, allele_id, phenotype_id)
                    assoc.add_source('WormBase:'+ref_num)
                    # eco_id = 'ECO:0000019'  # RNAi evidence  # TODO unused
                    assoc.add_association_to_graph()

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def process_pub_xrefs(self, limit=None):

        raw = '/'.join((self.rawdir, self.files['pub_xrefs']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        model = Model(g)
        logger.info("Processing publication xrefs")
        line_counter = 0
        with open(raw, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (wb_ref, xref) = row
                # WBPaper00000009 pmid8805<BR>
                # WBPaper00000011 doi10.1139/z78-244<BR>
                # WBPaper00000012 cgc12<BR>

                if self.testMode and wb_ref not in self.test_ids['pub']:
                    continue

                ref_id = 'WormBase:'+wb_ref
                xref_id = None
                r = None
                xref = re.sub(r'<BR>', '', xref)
                xref = xref.strip()
                if re.match(r'pmid', xref):
                    xref_id = 'PMID:'+re.sub(r'pmid\s*', '', xref)
                    reference = Reference(
                        g, xref_id, Reference.ref_types['journal_article'])
                elif re.search(r'[\(\)\<\>\[\]\s]', xref):
                    continue
                elif re.match(r'doi', xref):
                    xref_id = 'DOI:'+re.sub(r'doi', '', xref.strip())
                    reference = Reference(g, xref_id)
                elif re.match(r'cgc', xref):
                    # TODO not sure what to do here with cgc xrefs
                    continue
                else:
                    # logger.debug("Other xrefs like %s", xref)
                    continue

                if xref_id is not None:
                    reference.addRefToGraph()
                    model.addSameIndividual(ref_id, xref_id)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def process_feature_loc(self, limit):

        raw = '/'.join((self.rawdir, self.files['feature_loc']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("Processing Feature location and attributes")
        line_counter = 0
        geno = Genotype(g)
        strain_to_variant_map = {}
        build_num = self.version_num
        build_id = 'WormBase:'+build_num
        with gzip.open(raw, 'rb') as csvfile:
            filereader = csv.reader(
                io.TextIOWrapper(csvfile, newline=""), delimiter='\t',
                quotechar='\"')
            for row in filereader:
                if re.match(r'\#', ''.join(row)):
                    continue
                (chrom, db, feature_type_label, start, end, score, strand,
                 phase, attributes) = row

# I	interpolated_pmap_position	gene	1	559768	.	.	.	ID=gmap:spe-13;gmap=spe-13;status=uncloned;Note=-21.3602 cM (+/- 1.84 cM)
# I	WormBase	gene	3747	3909	.	-	.	ID=Gene:WBGene00023193;Name=WBGene00023193;interpolated_map_position=-21.9064;sequence_name=Y74C9A.6;biotype=snoRNA;Alias=Y74C9A.6
# I	absolute_pmap_position	gene	4119	10230	.	.	.	ID=gmap:homt-1;gmap=homt-1;status=cloned;Note=-21.8252 cM (+/- 0.00 cM)

                # dbs = re.split(
                #   r' ', 'assembly_component expressed_sequence_match Coding_transcript Genomic_canonical Non_coding_transcript Orfeome Promoterome Pseudogene RNAi_primary RNAi_secondary Reference Transposon Transposon_CDS cDNA_for_RNAi miRanda ncRNA operon polyA_signal_sequence polyA_site snlRNA')
                #
                # if db not in dbs:
                #     continue

                if feature_type_label not in [
                        'gene', 'point_mutation', 'deletion', 'RNAi_reagent',
                        'duplication', 'enhancer', 'binding_site',
                        'biological_region', 'complex_substitution',
                        'substitution', 'insertion', 'inverted_repeat']:
                    # note biological_regions include balancers
                    # other options here: promoter, regulatory_region, reagent
                    continue
                line_counter += 1

                attribute_dict = {}
                if attributes != '':
                    attribute_dict = dict(
                        item.split("=")for item in
                        re.sub(r'"', '', attributes).split(";"))

                fid = flabel = desc = None
                if 'ID' in attribute_dict:
                    fid = attribute_dict.get('ID')
                    if re.search(r'WB(Gene|Var|sf)', fid):
                        fid = re.sub(r'^\w+:WB', 'WormBase:WB', fid)
                    elif re.match(r'(gmap|landmark)', fid):
                        continue
                    else:
                        logger.info('other identifier %s', fid)
                        fid = None
                elif 'variation' in attribute_dict:
                    fid = 'WormBase:'+attribute_dict.get('variation')
                    flabel = attribute_dict.get('public_name')
                    sub = attribute_dict.get('substitution')
                    ins = attribute_dict.get('insertion')
                    # if it's a variation:
                    # variation=WBVar00604246;public_name=gk320600;strain=VC20384;substitution=C/T
                    desc = ''
                    if sub is not None:
                        desc = 'substitution='+sub
                    if ins is not None:
                        desc = 'insertion='+ins

                    # keep track of the strains with this variation,
                    # for later processing
                    strain_list = attribute_dict.get('strain')
                    if strain_list is not None:
                        for s in re.split(r',', strain_list):
                            if s.strip() not in strain_to_variant_map:
                                strain_to_variant_map[s.strip()] = set()
                            strain_to_variant_map[s.strip()].add(fid)

                # if feature_type_label == 'RNAi_reagent':
                    # Target=WBRNAi00096030 1 4942
                    # this will tell us where the RNAi is actually binding
                    # target = attribute_dict.get('Target') # TODO unused
                    # rnai_num = re.split(r' ', target)[0]  # TODO unused
                    # it will be the reagent-targeted-gene that has a position,
                    # (i think)
                    # TODO finish the RNAi binding location

                name = attribute_dict.get('Name')
                polymorphism = attribute_dict.get('polymorphism')

                if fid is None:
                    if name is not None and re.match(r'WBsf', name):
                        fid = 'WormBase:'+name
                        name = None
                    else:
                        continue

                if self.testMode \
                        and re.sub(r'WormBase:', '', fid) \
                        not in self.test_ids['gene']+self.test_ids['allele']:
                    continue

                # these really aren't that interesting
                if polymorphism is not None:
                    continue

                if name is not None and not re.search(name, fid):
                    if flabel is None:
                        flabel = name
                    else:
                        model.addSynonym(fid, name)

                if desc is not None:
                    model.addDescription(fid, desc)

                alias = attribute_dict.get('Alias')

                biotype = attribute_dict.get('biotype')
                note = attribute_dict.get('Note')
                other_name = attribute_dict.get('other_name')
                for n in [alias, other_name]:
                    if n is not None:
                        model.addSynonym(fid, other_name)

                ftype = self.get_feature_type_by_class_and_biotype(
                    feature_type_label, biotype)

                chr_id = makeChromID(chrom, build_id, 'CHR')
                geno.addChromosomeInstance(chrom, build_id, build_num)

                feature = Feature(g, fid, flabel, ftype)
                feature.addFeatureStartLocation(start, chr_id, strand)
                feature.addFeatureEndLocation(start, chr_id, strand)

                feature_is_class = False
                if feature_type_label == 'gene':
                    feature_is_class = True

                feature.addFeatureToGraph(True, None, feature_is_class)

                if note is not None:
                    model.addDescription(fid, note)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

                # RNAi reagents:
# I	RNAi_primary	RNAi_reagent	4184	10232	.	+	.	Target=WBRNAi00001601 1 6049 +;laboratory=YK;history_name=SA:yk326e10
# I	RNAi_primary	RNAi_reagent	4223	10147	.	+	.	Target=WBRNAi00033465 1 5925 +;laboratory=SV;history_name=MV_SV:mv_G_YK5052
# I	RNAi_primary	RNAi_reagent	5693	9391	.	+	.	Target=WBRNAi00066135 1 3699 +;laboratory=CH

                # TODO TF bindiing sites and network:
# I	TF_binding_site_region	TF_binding_site	1861	2048	.	+	.	Name=WBsf292777;tf_id=WBTranscriptionFactor000025;tf_name=DAF-16
# I	TF_binding_site_region	TF_binding_site	3403	4072	.	+	.	Name=WBsf331847;tf_id=WBTranscriptionFactor000703;tf_name=DPL-1

        return

    def process_disease_association(self, limit):

        raw = '/'.join((self.rawdir, self.files['disease_assoc']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        model = Model(g)
        logger.info("Processing disease models")
        geno = Genotype(g)
        line_counter = 0
        worm_taxon = 'NCBITaxon:6239'
        with open(raw, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                if re.match(r'!', ''.join(row)):  # header
                    continue
                line_counter += 1
                (db, gene_num, gene_symbol, is_not, disease_id, ref,
                 eco_symbol, with_or_from, aspect, gene_name, gene_synonym,
                 gene_class, taxon, date, assigned_by, blank, blank2) = row

                if self.testMode and gene_num not in self.test_ids['gene']:
                    continue

                # TODO add NOT phenotypes
                if is_not == 'NOT':
                    continue

                # WB	WBGene00000001	aap-1		DOID:2583	PMID:19029536	IEA	ENSEMBL:ENSG00000145675|OMIM:615214	D		Y110A7A.10	gene	taxon:6239	20150612	WB
                gene_id = 'WormBase:'+gene_num

                # make a variant of the gene
                vl = '_:'+'-'.join((gene_num, 'unspecified'))
                vl_label = 'some variant of '+gene_symbol
                geno.addAffectedLocus(vl, gene_id)
                model.addBlankNodeAnnotation(vl)
                animal_id = geno.make_experimental_model_with_genotype(
                    vl, vl_label, worm_taxon, 'worm')

                assoc = G2PAssoc(
                    g, self.name, animal_id,
                    disease_id, model.object_properties['model_of'])
                ref = re.sub(r'WB_REF:', 'WormBase:', ref)
                if ref != '':
                    assoc.add_source(ref)
                eco_id = None
                if eco_symbol == 'IEA':
                    eco_id = 'ECO:0000501'  # IEA is this now
                if eco_id is not None:
                    assoc.add_evidence(eco_id)

                assoc.add_association_to_graph()

        return

    def process_gene_interaction(self, limit):
        """
        The gene interaction file includes identified interactions,
        that are between two or more gene (products).
        In the case of interactions with >2 genes, this requires creating
        groups of genes that are involved in the interaction.
        From the wormbase help list: In the example WBInteraction000007779
        it would likely be misleading to suggest that lin-12 interacts with
        (suppresses in this case) smo-1 ALONE or that lin-12 suppresses let-60
        ALONE; the observation in the paper; see Table V in paper PMID:15990876
        was that a lin-12 allele (heterozygous lin-12(n941/+)) could suppress
        the "multivulva" phenotype induced synthetically by simultaneous
        perturbation of BOTH smo-1 (by RNAi) AND let-60 (by the n2021 allele).
        So this is necessarily a three-gene interaction.

        Therefore, we can create groups of genes based on their "status" of
        Effector | Effected.

        Status:  IN PROGRESS

        :param limit:
        :return:

        """

        raw = '/'.join((self.rawdir, self.files['gene_interaction']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("Processing gene interaction associations")
        line_counter = 0

        with gzip.open(raw, 'rb') as csvfile:
            filereader = csv.reader(
                io.TextIOWrapper(csvfile, newline=""), delimiter='\t',
                quotechar="'")

            for row in filereader:
                line_counter += 1
                if re.match(r'#', ''.join(row)):
                    continue

                (interaction_num, interaction_type, interaction_subtype,
                 summary, citation) = row[0:5]
                # print(row)
                interaction_id = 'WormBase:'+interaction_num

                # TODO deal with subtypes
                interaction_type_id = None
                if interaction_type == 'Genetic':
                    interaction_type_id = \
                        InteractionAssoc.interaction_object_properties[
                            'genetically_interacts_with']
                elif interaction_type == 'Physical':
                    interaction_type_id = \
                        InteractionAssoc.interaction_object_properties[
                            'molecularly_interacts_with']
                elif interaction_type == 'Regulatory':
                    interaction_type_id = \
                        InteractionAssoc.interaction_object_properties[
                            'regulates']
                else:
                    logger.info(
                        "An interaction type I don't understand %s",
                        interaction_type)

                num_interactors = (len(row) - 5) / 3
                if num_interactors != 2:
                    logger.info(
                        "Skipping interactions with !=2 participants:\n %s",
                        str(row))
                    continue

                gene_a_id = 'WormBase:'+row[5]
                gene_b_id = 'WormBase:'+row[8]

                if self.testMode \
                        and gene_a_id not in self.test_ids['gene'] \
                        and gene_b_id not in self.test_ids['gene']:
                    continue

                assoc = InteractionAssoc(
                    g, self.name, gene_a_id, gene_b_id, interaction_type_id)
                assoc.set_association_id(interaction_id)
                assoc.add_association_to_graph()
                assoc_id = assoc.get_association_id()
                # citation is not a pmid or WBref - get this some other way
                model.addDescription(assoc_id, summary)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def get_feature_type_by_class_and_biotype(self, ftype, biotype):
        ftype_id = None
        biotype_map = {
            'lincRNA': 'SO:0001641',
            'miRNA': 'SO:0001265',
            'ncRNA': 'SO:0001263',
            'piRNA': 'SO:0001638',
            'rRNA': 'SO:0001637',
            'scRNA': 'SO:0001266',
            'snRNA': 'SO:0001268',
            'snoRNA': 'SO:0001267',
            'tRNA': 'SO:0001272',
            # transposable element gene
            'transposon_protein_coding': 'SO:0000111',
            'transposon_pseudogene': 'SO:0001897',
            'pseudogene': 'SO:0000336',
            'protein_coding': 'SO:0001217',
            'asRNA': 'SO:0001263',  # using ncRNA gene  TODO make term request
        }

        ftype_map = {
            'point_mutation': 'SO:1000008',
            'deletion': 'SO:0000159',
            'RNAi_reagent': 'SO:0000337',
            'duplication': 'SO:1000035',
            'enhancer': 'SO:0000165',
            'binding_site': 'SO:0000409',
            'biological_region': 'SO:0001411',
            'complex_substitution': 'SO:1000005'
        }
        if ftype == 'gene':
            if biotype in biotype_map:
                ftype_id = biotype_map.get(biotype)
        else:
            ftype_id = ftype_map.get(ftype)

        return ftype_id

    def make_reagent_targeted_gene_id(
            self, gene_id, reagent_id):

        rtg_id = '_:'+'-'.join((gene_id, reagent_id))
        # TODO targeted_gene_id unused
        # targeted_gene_id = re.sub(r'W(orm)?B(ase)?:', '', rtg_id)

        return rtg_id

    def getTestSuite(self):
        import unittest
        from tests.test_wormbase import WormBaseTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(WormBaseTestCase)

        return test_suite
