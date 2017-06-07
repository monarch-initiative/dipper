import csv
import logging
import re
import gzip

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Model import Model
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.GenomicFeature import Feature, makeChromID


logger = logging.getLogger(__name__)
AQDL = 'http://www.animalgenome.org/QTLdb'


class AnimalQTLdb(Source):
    """
    The Animal Quantitative Trait Loci (QTL) database (Animal QTLdb)
    is designed to house publicly all available QTL and
    single-nucleotide polymorphism/gene association data on livestock
    animal species.  This includes:
        * chicken
        * horse
        * cow
        * sheep
        * rainbow trout
        * pig
    While most of the phenotypes here are related to animal husbandry,
    production, and rearing, integration of these phenotypes with other species
    may lead to insight for human disease.

    Here, we use the QTL genetic maps and their computed genomic locations to
    create associations between the QTLs and their traits.  The traits come in
    their internal Animal Trait ontology vocabulary, which they further map to
    [Vertebrate Trait](http://bioportal.bioontology.org/ontologies/VT),
    Product Trait, and Clinical Measurement Ontology vocabularies.

    Since these are only associations to broad locations,
    we link the traits via "is_marker_for", since there is no specific
    causative nature in the association.  p-values for the associations
    are attached to the Association objects.  We default to the UCSC build for
    the genomic coordinates, and make equivalences.

    Any genetic position ranges that are <0, we do not include here.

    """

    files = {
        # defaulting to this
        'cattle_bp': {
            'file': 'QTL_Btau_4.6.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_Btau_4.6.gff.txt.gz'},
        # disabling this for now
        # 'cattle_umd_bp': {
        #   'file': 'QTL_UMD_3.1.gff.txt.gz',
        #   'url': AQDL + '/tmp/QTL_UMD_3.1.gff.txt.gz'},
        'cattle_cm': {
            'file': 'cattle_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/cattle_QTLdata.txt'},
        'chicken_bp': {
            'file': 'QTL_GG_4.0.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_GG_4.0.gff.txt.gz'},
        'chicken_cm': {
            'file': 'chicken_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/chicken_QTLdata.txt'},
        'pig_bp': {
            'file': 'QTL_SS_10.2.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_SS_10.2.gff.txt.gz'},
        'pig_cm': {
            'file': 'pig_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/pig_QTLdata.txt'},
        'sheep_bp': {
            'file': 'QTL_OAR_3.1.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_OAR_3.1.gff.txt.gz'},
        'sheep_cm': {
            'file': 'sheep_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/sheep_QTLdata.txt'},
        'horse_bp': {
            'file': 'QTL_EquCab2.0.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_EquCab2.0.gff.txt.gz'},
        'horse_cm': {
            'file': 'horse_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/horse_QTLdata.txt'},
        'rainbow_trout_cm': {
            'file': 'rainbow_trout_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/rainbow_trout_QTLdata.txt'},
        # TODO add rainbow_trout_bp when available
        'trait_mappings': {
            'file': 'trait_mappings',
            'url': AQDL + '/export/trait_mappings.csv'}
    }

    # QTL ids
    test_ids = {
        28483, 29016, 29018, 8945, 29385, 12532, 31023, 14234, 17138, 1795,
        1798, 32133
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'animalqtldb')

        # update the dataset object with details about this resource
        self.dataset = Dataset(
            'animalqtldb', 'Animal QTL db',
            'http://www.animalgenome.org/cgi-bin/QTLdb/index', None, None,
            AQDL + '/faq#23', graph_type=graph_type)

        # source-specific warnings.  will be cleared when resolved.
        logger.warning(
            "No licences or rights exist for the raw data from this resource.")

        return

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        """

        :param limit:
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %s rows fo each file", str(limit))

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True
            g = self.testgraph
        else:
            g = self.graph

        tmap = '/'.join((self.rawdir, self.files['trait_mappings']['file']))
        self._process_trait_mappings(tmap, limit)

        geno = Genotype(g)
        # organisms  = ['chicken']
        organisms = [
            'chicken', 'pig', 'horse', 'rainbow_trout', 'sheep', 'cattle']

        for o in organisms:
            tax_id = self._get_tax_by_common_name(o)
            geno.addGenome(tax_id, o)
            build_id = None
            build = None

            k = o + '_bp'
            if k in self.files:
                file = self.files[k]['file']
                m = re.search(r'QTL_([\w\.]+)\.gff.txt.gz', file)
                if m is None:
                    logger.error("Can't match a gff build")
                else:
                    build = m.group(1)
                    build_id = self._map_build_by_abbrev(build)
                    logger.info("Build = %s", build_id)
                    geno.addReferenceGenome(build_id, build, tax_id)
                if build_id is not None:
                    self._process_QTLs_genomic_location(
                        '/'.join((self.rawdir, file)), tax_id, build_id, build,
                        limit)

            k = o+'_cm'
            if k in self.files:
                file = self.files[k]['file']
                self._process_QTLs_genetic_location(
                    '/'.join((self.rawdir, file)), tax_id, o, limit)

        logger.info("Finished parsing")
        return

    def _process_QTLs_genetic_location(
            self, raw, taxon_id, common_name, limit=None):
        """
        This function processes

        Triples created:

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        model = Model(g)
        eco_id = "ECO:0000061"  # Quantitative Trait Analysis Evidence

        logger.info(
            "Processing genetic location for %s from %s", taxon_id, raw)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (qtl_id,
                 qtl_symbol,
                 trait_name,
                 assotype,
                 empty,
                 chromosome,
                 position_cm,
                 range_cm,
                 flankmark_a2,
                 flankmark_a1,
                 peak_mark,
                 flankmark_b1,
                 flankmark_b2,
                 exp_id,
                 model_id,
                 test_base,
                 sig_level,
                 lod_score,
                 ls_mean,
                 p_values,
                 f_statistics,
                 variance,
                 bayes_value,
                 likelihood_ratio,
                 trait_id, dom_effect,
                 add_effect,
                 pubmed_id,
                 gene_id,
                 gene_id_src,
                 gene_id_type,
                 empty2) = row

                if self.testMode and int(qtl_id) not in self.test_ids:
                    continue

                qtl_id = 'AQTL:'+qtl_id.strip()
                trait_id = 'AQTLTrait:'+trait_id.strip()

                # Add QTL to graph
                feature = Feature(g, qtl_id, qtl_symbol, geno.genoparts['QTL'])
                feature.addTaxonToFeature(taxon_id)

                # deal with the chromosome
                chrom_id = makeChromID(chromosome, taxon_id, 'CHR')

                # add a version of the chromosome which is defined as
                # the genetic map
                build_id = 'MONARCH:'+common_name.strip()+'-linkage'
                build_label = common_name+' genetic map'
                geno.addReferenceGenome(build_id, build_label, taxon_id)
                chrom_in_build_id = makeChromID(
                    chromosome, build_id, 'MONARCH')
                geno.addChromosomeInstance(
                    chromosome, build_id, build_label, chrom_id)
                start = stop = None
                # range_cm sometimes ends in "(Mb)"  (i.e pig 2016 Nov)
                range_mb = re.split(r'\(', range_cm)
                if range_mb is not None:
                    range_cm = range_mb[0]

                if re.search(r'[0-9].*-.*[0-9]', range_cm):
                    range_parts = re.split(r'-', range_cm)

                    # check for poorly formed ranges
                    if len(range_parts) == 2 and\
                            range_parts[0] != '' and range_parts[1] != '':
                        (start, stop) = \
                            [int(float(x.strip()))
                             for x in re.split(r'-', range_cm)]
                    else:
                        logger.info(
                            "A cM range we can't handle for QTL %s: %s",
                            qtl_id, range_cm)
                elif position_cm != '':
                    match = re.match(r'([0-9]*\.[0-9]*)', position_cm)
                    if match is not None:
                        position_cm = match.group()
                        start = stop = int(float(position_cm))

                # FIXME remove converion to int for start/stop
                # when schema can handle floats add in the genetic location
                # based on the range
                feature.addFeatureStartLocation(
                    start, chrom_in_build_id, None,
                    [Feature.types['FuzzyPosition']])
                feature.addFeatureEndLocation(
                    stop, chrom_in_build_id, None,
                    [Feature.types['FuzzyPosition']])
                feature.addFeatureToGraph()

                # sometimes there's a peak marker, like a rsid.
                # we want to add that as a variant of the gene,
                # and xref it to the qtl.
                dbsnp_id = None
                if peak_mark != '' and peak_mark != '.' and \
                        re.match(r'rs', peak_mark.strip()):
                    dbsnp_id = 'dbSNP:'+peak_mark.strip()

                    model.addIndividualToGraph(
                        dbsnp_id, None,
                        geno.genoparts['sequence_alteration'])
                    model.addXref(qtl_id, dbsnp_id)

                gene_id = gene_id.replace('uncharacterized ', '')
                if gene_id is not None and gene_id != '' and gene_id != '.'\
                        and re.fullmatch(r'[^ ]*', gene_id) is not None:

                    # we assume if no src is provided
                    # and gene_id is an integer, it's NCBI
                    if (gene_id_src == 'NCBIgene' or gene_id_src == '') and \
                            gene_id.strip().isdigit() :
                        gene_id = 'NCBIGene:' + gene_id.strip()
                        # we will expect that these labels provided elsewhere
                        geno.addGene(gene_id, None)
                        # FIXME what is the right relationship here?
                        geno.addAffectedLocus(qtl_id, gene_id)

                        if dbsnp_id is not None:
                            # add the rsid as a seq alt of the gene_id
                            vl_id = \
                                '_:' + re.sub(
                                    r':', '', gene_id) + '-' + peak_mark.strip()
                            geno.addSequenceAlterationToVariantLocus(
                                dbsnp_id, vl_id)
                            geno.addAffectedLocus(vl_id, gene_id)

                # add the trait
                model.addClassToGraph(trait_id, trait_name)

                # Add publication
                reference = None
                if re.match(r'ISU.*', pubmed_id):
                    pub_id = 'AQTLPub:'+pubmed_id.strip()
                    reference = Reference(g, pub_id)
                elif pubmed_id != '':
                    pub_id = 'PMID:'+pubmed_id.strip()
                    reference = Reference(g,
                        pub_id, Reference.ref_types['journal_article'])

                if reference is not None:
                    reference.addRefToGraph()

                # make the association to the QTL
                assoc = G2PAssoc(
                    g, self.name, qtl_id, trait_id,
                    model.object_properties['is_marker_for'])
                assoc.add_evidence(eco_id)
                assoc.add_source(pub_id)

                # create a description from the contents of the file
                # desc = ''

                # assoc.addDescription(g, assoc_id, desc)

                # TODO add exp_id as evidence
                # if exp_id != '':
                #     exp_id = 'AQTLExp:'+exp_id
                #     gu.addIndividualToGraph(g, exp_id, None, eco_id)

                if p_values != '':
                    s = re.sub(r'<', '', p_values)
                    s = re.sub(r',', '.', s)  # international notation
                    if s.isnumeric():
                        score = float(s)
                        assoc.set_score(score)  # todo add score type
                # TODO add LOD score?
                assoc.add_association_to_graph()

                # make the association to the dbsnp_id, if found
                if dbsnp_id is not None:
                    # make the association to the dbsnp_id
                    assoc = G2PAssoc(
                        g, self.name, dbsnp_id, trait_id,
                        model.object_properties['is_marker_for'])
                    assoc.add_evidence(eco_id)
                    assoc.add_source(pub_id)

                    # create a description from the contents of the file
                    # desc = ''
                    # assoc.addDescription(g, assoc_id, desc)

                    # TODO add exp_id
                    # if exp_id != '':
                    #     exp_id = 'AQTLExp:'+exp_id
                    #     gu.addIndividualToGraph(g, exp_id, None, eco_id)

                    if p_values != '':
                        s = re.sub(r'<', '', p_values)
                        s = re.sub(r',', '.', s)
                        if s.isnumeric():
                            score = float(s)
                            assoc.set_score(score)  # todo add score type
                    # TODO add LOD score?

                    assoc.add_association_to_graph()

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        logger.info("Done with QTL genetic info")
        return

    def _process_QTLs_genomic_location(
            self, raw, taxon_id, build_id, build_label, limit=None):
        """
        This method

        Triples created:

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        geno = Genotype(g)
        # assume that chrs get added to the genome elsewhere
        # genome_id = geno.makeGenomeID(taxon_id)  # TODO unused

        eco_id = "ECO:0000061"  # Quantitative Trait Analysis Evidence
        logger.info("Processing QTL locations for %s", taxon_id)
        with gzip.open(raw, 'rt', encoding='ISO-8859-1') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            # bad_attr_flag = False  # TODO unused
            for row in reader:
                line_counter += 1
                if re.match(r'^#', ' '.join(row)):
                    continue

                (chromosome, qtl_source, qtl_type, start_bp, stop_bp, frame,
                 strand, score, attr) = row

                # Chr.Z   Animal QTLdb    Production_QTL  33954873      34023581        .       .       .
                # QTL_ID=2242;Name="Spleen percentage";Abbrev="SPLP";PUBMED_ID=17012160;trait_ID=2234;
                # trait="Spleen percentage";breed="leghorn";"FlankMarkers=ADL0022";VTO_name="spleen mass";
                # CMO_name="spleen weight to body weight ratio";Map_Type="Linkage";Model="Mendelian";
                # Test_Base="Chromosome-wise";Significance="Significant";P-value="<0.05";F-Stat="5.52";
                # Variance="2.94";Dominance_Effect="-0.002";Additive_Effect="0.01"

                # make dictionary of attributes
                # keys are:
                # QTL_ID,Name,Abbrev,PUBMED_ID,trait_ID,trait,FlankMarkers,
                # VTO_name,Map_Type,Significance,P-value,Model,
                # Test_Base,Variance, Bayes-value,PTO_name,gene_IDsrc,peak_cM,
                # CMO_name,gene_ID,F-Stat,LOD-score,Additive_Effect,
                # Dominance_Effect,Likelihood_Ratio,LS-means,Breed,
                # trait (duplicate with Name),Variance,Bayes-value,
                # F-Stat,LOD-score,Additive_Effect,Dominance_Effect,
                # Likelihood_Ratio,LS-means

                # deal with poorly formed attributes
                if re.search(r'"FlankMarkers";', attr):
                    attr = re.sub(r'FlankMarkers;', '', attr)
                attr_items = re.sub(r'"', '', attr).split(";")
                bad_attrs = set()
                for a in attr_items:
                    if not re.search(r'=', a):
                        # bad_attr_flag = True  # TODO unused
                        # remove this attribute from the list
                        bad_attrs.add(a)

                attr_set = set(attr_items) - bad_attrs
                attribute_dict = dict(item.split("=") for item in attr_set)

                qtl_num = attribute_dict.get('QTL_ID')
                if self.testMode and int(qtl_num) not in self.test_ids:
                    continue

                # make association between QTL and trait
                qtl_id = 'AQTL:' + str(qtl_num)
                model.addIndividualToGraph(qtl_id, None, geno.genoparts['QTL'])
                geno.addTaxon(taxon_id, qtl_id)

                trait_id = 'AQTLTrait:'+attribute_dict.get('trait_ID')

                # if pub is in attributes, add it to the association
                pub_id = None
                if 'PUBMED_ID' in attribute_dict.keys():
                    pub_id = attribute_dict.get('PUBMED_ID')
                    if re.match(r'ISU.*', pub_id):
                        pub_id = 'AQTLPub:' + pub_id.strip()
                        reference = Reference(g, pub_id)
                    else:
                        pub_id = 'PMID:' + pub_id.strip()
                        reference = Reference(
                            g, pub_id, Reference.ref_types['journal_article'])
                    reference.addRefToGraph()

                # Add QTL to graph
                assoc = G2PAssoc(
                    g, self.name, qtl_id, trait_id,
                    model.object_properties['is_marker_for'])
                assoc.add_evidence(eco_id)
                assoc.add_source(pub_id)
                if 'P-value' in attribute_dict.keys():
                    s = re.sub(r'<', '', attribute_dict.get('P-value'))
                    if ',' in s:
                        s = re.sub(r',', '.', s)
                    if s.isnumeric():
                        score = float(s)
                        assoc.set_score(score)

                assoc.add_association_to_graph()
                # TODO make association to breed
                # (which means making QTL feature in Breed background)

                # get location of QTL
                chromosome = re.sub(r'Chr\.', '', chromosome)
                chrom_id = makeChromID(chromosome, taxon_id, 'CHR')

                chrom_in_build_id = \
                    makeChromID(chromosome, build_id, 'MONARCH')
                geno.addChromosomeInstance(
                    chromosome, build_id, build_label, chrom_id)
                qtl_feature = Feature(g, qtl_id, None, geno.genoparts['QTL'])
                if start_bp == '':
                    start_bp = None
                qtl_feature.addFeatureStartLocation(
                    start_bp, chrom_in_build_id, strand,
                    [Feature.types['FuzzyPosition']])
                if stop_bp == '':
                    stop_bp = None
                qtl_feature.addFeatureEndLocation(
                    stop_bp, chrom_in_build_id, strand,
                    [Feature.types['FuzzyPosition']])
                qtl_feature.addTaxonToFeature(taxon_id)
                qtl_feature.addFeatureToGraph()

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        logger.warning("Bad attribute flags in this file")
        logger.info("Done with QTL genomic mappings for %s", taxon_id)
        return

    def _process_trait_mappings(self, raw, limit=None):
        """
        This method mapps traits from/to ...

        Triples created:

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        model = Model(g)

        # with open(raw, 'r') as csvfile:
        #     filereader = csv.reader(csvfile, delimiter=',')
        #     row_count = sum(1 for row in filereader)
        #     row_count = row_count - 1

        with open(raw, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip header line
            for row in filereader:
                line_counter += 1
                # need to skip the last line
                if len(row) < 8:
                    logger.info(
                        "skipping line %d: %s", line_counter, '\t'.join(row))
                    continue
                (vto_id, pto_id, cmo_id, ato_column, species, trait_class,
                 trait_type, qtl_count) = row

                ato_id = re.sub(r'ATO #', 'AQTLTrait:',
                                re.sub(r'\].*', '',
                                       re.sub(r'\[', '', ato_column)))
                ato_id = ato_id.strip()

                ato_label = re.sub(r'.*\]\s*', '', ato_column)

                # if species == 'Cattle':
                #     ato_id = re.sub(r'ATO:', 'AQTLTraitCattle:', ato_id)
                # elif species == 'Chicken':
                #     ato_id = re.sub(r'ATO:', 'AQTLTraitChicken:', ato_id)
                # elif species == 'Sheep':
                #     ato_id = re.sub(r'ATO:', 'AQTLTraitSheep:', ato_id)
                # elif species == 'Horse':
                #     ato_id = re.sub(r'ATO:', 'AQTLTraitHorse:', ato_id)
                # elif species == 'Pig':
                #     ato_id = re.sub(r'ATO:', 'AQTLTraitPig:', ato_id)
                # elif species == 'Rainbow trout':
                #     ato_id = re.sub(
                #       r'ATO:', 'AQTLTraitRainbowTrout:', ato_id)
                # else:
                #     logger.warning(
                #       'Unknown species %s foufnd in trait mapping file.',
                #       species)
                #     continue
                # print(ato_label)

                model.addClassToGraph(ato_id, ato_label.strip())

                if re.match(r'VT:.*', vto_id):
                    model.addClassToGraph(vto_id, None)
                    model.addEquivalentClass(ato_id, vto_id)
                if re.match(r'LPT:.*', pto_id):
                    model.addClassToGraph(pto_id, None)
                    model.addXref(ato_id, pto_id)
                if re.match(r'CMO:.*', cmo_id):
                    model.addClassToGraph(cmo_id, None)
                    model.addXref(ato_id, cmo_id)

        logger.info("Done with trait mappings")
        return

    # TODO consider static functions
    def _get_tax_by_common_name(self, common_name):

        tax_map = {
            'chicken': 9031,
            'cattle': 9913,
            'pig': 9823,
            'sheep': 9940,
            'horse': 9796,
            'rainbow_trout': 8022}

        return 'NCBITaxon:'+str(tax_map[common_name])

    def _map_build_by_abbrev(self, build):

        build_id = None
        build_map = {
            'Btau_4.6': {
                'UCSC': 'bosTau7',
                'NCBIGenome': 82,
                'NCBIAssembly': 'GCF_000003205.5'},
            'UMD_3.1': {
                'NCBIAssembly': 'GCF_000003055.6',
                'NCBIGenome': 82},
            'GG_4.0': {
                'NCBIAssembly': '317958',
                'UCSC': 'galGal4'},
            'SS_10.2': {
                'NCBIAssembly': '304498',
                'UCSC': 'susScr3'},
            'OAR_3.1': {
                'UCSC': 'oviAri3',
                'NCBIAssembly': 'GCF_000298735.1',
                'NCBIGenome':  83},
            'EquCab2.0': {
                'UCSC': 'equCab2',
                'NCBIAssembly': 'GCF_000002305.2'}
        }

        b = build_map.get(build)
        if b is not None:
            if 'UCSC' in b.keys():
                build_id = 'UCSC:'+b.get('UCSC')
            else:
                build_id = 'NCBIAssembly:'+b.get('NCBIAssembly')

        return build_id

    def _map_linkage_by_organism(self, organism):
        """
        Need to add appropriate linkage maps...but need more information before
        they become identified.
        This is not yet confirmed, thus is not utilized.
        :param organism:
        :return:

        """

        # TODO this hash is unused and would be an external translation table
        # tax_map = {
            # 'chicken': 'Wageningen-chicken',    # PMID:10645958
            # 'cattle': 'USDA-MARC',              # PMID:9074927
            # 'pig': 'USDA-MARC',                 # PMID:8743988
            # # -- this is an assumption
            # 'sheep': 'Wageningen-sheep',        # PMID:11435411
            # 'horse': 'Swinburne-Penedo',        # PMID:16314071 PMID:16093715
            # 'rainbow_trout': 'USDA-NCCCWA'}     # PMID:19019240 PMID:22101344
        # return

    def getTestSuite(self):
        import unittest
        from tests.test_animalqtl import AnimalQTLdbTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(AnimalQTLdbTestCase)

        return test_suite
