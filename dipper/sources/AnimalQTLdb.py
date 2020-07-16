import csv
import logging
import re
import gzip

from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

#       https://www.animalgenome.org/tmp/QTL_EquCab2.0.gff.txt.gz'
# mapDwnLd36738TDWS.txt.gz
AQDL = 'https://www.animalgenome.org/QTLdb'
LOG = logging.getLogger(__name__)


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

    GENEINFO = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO'
    GITDIP = 'https://raw.githubusercontent.com/monarch-initiative/dipper/master'

    gff_columns = [  # GFF files (QTL_*.gff.txt.gz) should contain these columns
        'SEQNAME',
        'SOURCE',
        'FEATURE',
        'START',
        'END',
        'SCORE',
        'STRAND',
        'FRAME',
        'ATTRIBUTE'
    ]

    qtl_columns = [  # columns in *_QTLdata.txt files
        'qtl_id',
        'qtl_symbol',
        'trait_name',
        'assotype',
        'empty',
        'chromosome',
        'position_cm',
        'range_cm',
        'flankmark_a2',
        'flankmark_a1',
        'peak_mark',
        'flankmark_b1',
        'flankmark_b2',
        'exp_id',
        'model_id',
        'test_base',
        'psig_level',
        'lod_score',
        'ls_mean',
        'p_values',
        'f_statistics',
        'variance',
        'bayes_value',
        'likelihood_ratio',
        'trait_id',
        'dom_effect',
        'add_effect',
        'pubmed_id',
        'gene_id',
        'gene_id_src',
        'gene_id_type',
        'empty2'
    ]

    gene_info_columns = [  # columns from *.gene_info.gz files
        'tax_id',
        'GeneID',
        'Symbol',
        'LocusTag',
        'Synonyms',
        'dbXrefs',
        'chromosome',
        'map_location',
        'description',
        'type_of_gene',
        'Symbol_from_nomenclature_authority',
        'Full_name_from_nomenclature_authority',
        'Nomenclature_status',
        'Other_designations',
        'Modification_date',
        'Feature_type'
    ]

    trait_mapping_columns = [  # columns in file 'trait_mappings'
                               # (column info only used for this file at the moment)
        'VT',
        'LPT',
        'CMO',
        'ATO',
        'Species',
        'Class',
        'Type',
        'QTL_Count'
    ]

    files = {
        # defaulting to this
        'cattle_bp': {
            'file': 'QTL_Btau_4.6.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_Btau_4.6.gff.txt.gz',
            'curie': 'cattleQTL',
            'columns': gff_columns,
        },
        'cattle_cm': {
            'file': 'cattle_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/cattle_QTLdata.txt',
            'curie': 'cattleQTL',
            'columns': qtl_columns,
        },
        'chicken_bp': {
            'file': 'QTL_GG_4.0.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_GG_4.0.gff.txt.gz',
            'curie': 'chickenQTL',
            'columns': gff_columns,
        },
        'chicken_cm': {
            'file': 'chicken_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/chicken_QTLdata.txt',
            'curie': 'chickenQTL',
            'columns': qtl_columns,
        },
        'pig_bp': {
            'file': 'QTL_SS_10.2.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_SS_10.2.gff.txt.gz',
            'curie': 'pigQTL',
            'columns': gff_columns,
        },
        'pig_cm': {
            'file': 'pig_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/pig_QTLdata.txt',
            'curie': 'pigQTL',
            'columns': qtl_columns,
        },
        'sheep_bp': {
            'file': 'QTL_OAR_3.1.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_OAR_3.1.gff.txt.gz',
            'curie': 'sheepQTL',
            'columns': gff_columns,
        },
        'sheep_cm': {
            'file': 'sheep_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/sheep_QTLdata.txt',
            'curie': 'sheepQTL',
            'columns': qtl_columns,
        },
        'horse_bp': {
            'file': 'QTL_EquCab2.0.gff.txt.gz',
            'url': AQDL + '/tmp/QTL_EquCab2.0.gff.txt.gz',
            'curie': 'horseQTL',
            'columns': gff_columns,
        },
        'horse_cm': {
            'file': 'horse_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/horse_QTLdata.txt',
            'curie': 'horseQTL',
            'columns': qtl_columns,
        },
        'rainbow_trout_cm': {
            'file': 'rainbow_trout_QTLdata.txt',
            'url': AQDL + '/export/KSUI8GFHOT6/rainbow_trout_QTLdata.txt',
            'curie': 'rainbow_troutQTL',
            'columns': qtl_columns,
        },

        #                  Gene_info from NCBI
        # to reasure TEC that when we see an integer
        # it is a gene identifier from NCBI for the species
        # misses will not block, but they will squawk,
        # The last three are homemade, (in DipperCache)

        # pig  # "Sus scrofa"  # NCBITaxon:9823
        'Sus_scrofa_info': {
            'file': 'Sus_scrofa.gene_info.gz',
            'url': GENEINFO + '/Mammalia/Sus_scrofa.gene_info.gz',
            'columns': gene_info_columns,
        },
        # cattle  # "Bos taurus"      # NCBITaxon:9913
        'Bos_taurus_info': {
            'file': 'Bos_taurus.gene_info.gz',
            'url': GENEINFO + '/Mammalia/Bos_taurus.gene_info.gz',
            'columns': gene_info_columns,
        },
        # chicken  # "Gallus gallus"  # NCBITaxon:9031
        'Gallus_gallus_info': {
            'file': 'Gallus_gallus.gene_info.gz',
            'url': GENEINFO + '/Non-mammalian_vertebrates/Gallus_gallus.gene_info.gz',
            'columns': gene_info_columns,
        },
        # horse  # "Equus caballus"  # NCBITaxon:9796
        'Equus_caballus_info': {
            # This file isn't on NCBI's ftp site, so need to use the cached URL instead.
            'file': 'Equus_caballus.gene_info.gz',
            'url': Source.DIPPERCACHE + '/Equus_caballus.gene_info.gz',
            'columns': gene_info_columns,
        },
        # sheep  # "Ovis aries"  # NCBITaxon:9940
        'Ovis_aries_info': {
            # This file isn't on NCBI's ftp site, so need to use the cached URL instead.
            'file': 'Ovis_aries.gene_info.gz',
            'url': Source.DIPPERCACHE + '/Ovis_aries.gene_info.gz',
            'columns': gene_info_columns,
        },
        # rainbow trout  # "Oncorhynchus mykiss"  # NCBITaxon:8022
        'Oncorhynchus_mykiss_info': {
            # This file isn't on NCBI's ftp site, so need to use the cached URL instead.
            'file': 'Oncorhynchus_mykiss.gene_info.gz',
            'url': Source.DIPPERCACHE + '/Oncorhynchus_mykiss.gene_info.gz',
            'columns': gene_info_columns,
        },
        ########################################
        # TODO add rainbow_trout_bp when available

        'trait_mappings': {
            'file': 'trait_mappings.csv',
            'url': AQDL + '/export/trait_mappings.csv',
            'columns': trait_mapping_columns,
        },
    }

    # AQTL ids
    test_ids = {
        28483, 29016, 29018, 8945, 29385, 12532, 31023, 14234, 17138, 1795, 1798, 32133
    }

    def __init__(self, graph_type, are_bnodes_skolemized, data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='animalqtldb',
            ingest_title='Animal QTL db',
            ingest_url='http://www.animalgenome.org/cgi-bin/QTLdb/index',
            ingest_logo='source-animalqtldb.png',
            license_url=None,
            data_rights="'" + AQDL + '/faq#32' + "'"
            # file_handle=None
        )

        self.gene_info = list()
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
            LOG.info("Only parsing first %s rows fo each file", str(limit))

        if self.test_only:
            self.test_mode = True
            graph = self.testgraph
        else:
            graph = self.graph

        trait_src_key = 'trait_mappings'
        traitmap = '/'.join((self.rawdir, self.files[trait_src_key]['file']))

        LOG.info("Parsing trait mapping  file %s", traitmap)
        self._process_trait_mappings(traitmap, trait_src_key, limit)

        geno = Genotype(graph)
        animals = ['chicken', 'pig', 'horse', 'rainbow_trout', 'sheep', 'cattle']

        for common_name in animals:
            txid_num = self.resolve(common_name).split(':')[1]
            taxon_label = self.localtt[common_name]
            taxon_curie = self.globaltt[taxon_label]
            taxon_num = taxon_curie.split(':')[1]
            txid_num = taxon_num  # for now
            taxon_word = taxon_label.replace(' ', '_')
            src_key = taxon_word + '_info'
            gene_info_file = '/'.join((
                self.rawdir, self.files[src_key]['file']))
            self.gene_info = list()
            LOG.info('Ingesting %s', gene_info_file)
            with gzip.open(gene_info_file, 'rt') as gi_gz:
                filereader = csv.reader(gi_gz, delimiter='\t')
                # skipping header checking, b/c not all of these gene_info files have
                # headers
                col = self.files[src_key]['columns']
                for row in filereader:
                    if re.match('^#', row[0][0]):
                        continue
                    if len(row) != len(col):
                        LOG.warning(
                            "Problem parsing in %s row %s\n"
                            "got %s cols but expected %s",
                            gene_info_file, row, len(row), len(col))
                    self.gene_info.append(row[col.index('GeneID')])
            LOG.info(
                'Gene Info for %s has %i entries', common_name, len(self.gene_info))
            build = None

            fname_bp = common_name + '_bp'
            if fname_bp in self.files:
                bpfile = self.files[fname_bp]['file']
                mch = re.search(r'QTL_([\w\.]+)\.gff.txt.gz', bpfile)
                if mch is None:
                    LOG.error("Can't match a gff build to " + fname_bp)
                else:
                    build = mch.group(1)
                    build_id = self.localtt[build]
                    LOG.info("Build UCSC label is: %s", build_id)

                    geno.addReferenceGenome(build_id, build, txid_num)

                if build_id is not None:
                    self._process_qtls_genomic_location(
                        '/'.join((self.rawdir, bpfile)),
                        fname_bp,
                        txid_num,
                        build_id,
                        build,
                        common_name,
                        limit
                    )

            fname_cm = common_name + '_cm'
            if fname_cm in self.files:
                cmfile = self.files[fname_cm]['file']
                self._process_qtls_genetic_location(
                    '/'.join((self.rawdir, cmfile)),
                    fname_cm,
                    txid_num,
                    common_name,
                    limit)

        LOG.info("Finished parsing")
        return

    def _process_qtls_genetic_location(
            self, raw, src_key, txid, common_name, limit=None):
        """
        This function processes

        Triples created:

        :param limit:
        :return:

        """
        aql_curie = self.files[src_key]['curie']

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        model = Model(graph)
        eco_id = self.globaltt['quantitative trait analysis evidence']
        taxon_curie = 'NCBITaxon:' + txid

        LOG.info("Processing genetic location for %s from %s", taxon_curie, raw)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            # no header in these files, so no header checking
            col = self.files[src_key]['columns']
            for row in reader:

                if len(row) != len(self.qtl_columns):
                    LOG.warning(
                        "Problem parsing %s line %i containing: \n%s\n"
                        "got %i cols but expected %i",
                        raw, reader.line_num, row, len(row), len(col))
                    continue
                else:
                    qtl_id = row[col.index('qtl_id')].strip()
                    qtl_symbol = row[col.index('qtl_symbol')].strip()
                    trait_name = row[col.index('trait_name')].strip()
                    # assotype = row[col.index('assotype')].strip()
                    # empty = row[col.index('empty')].strip()
                    chromosome = row[col.index('chromosome')].strip()
                    position_cm = row[col.index('position_cm')].strip()
                    range_cm = row[col.index('range_cm')].strip()
                    # flankmark_a2 = row[col.index('flankmark_a2')].strip()
                    # flankmark_a1 = row[col.index('flankmark_a1')].strip()
                    peak_mark = row[col.index('peak_mark')].strip()
                    # flankmark_b1 = row[col.index('flankmark_b1')].strip()
                    # flankmark_b2 = row[col.index('flankmark_b2')].strip()
                    # exp_id = row[col.index('exp_id')].strip()
                    # model_id = row[col.index('model_id')].strip()
                    # test_base = row[col.index('test_base')].strip()
                    # sig_level = row[col.index('sig_level')].strip()
                    # lod_score = row[col.index('lod_score')].strip()
                    # ls_mean = row[col.index('ls_mean')].strip()
                    p_values = row[col.index('p_values')].strip()
                    # f_statistics = row[col.index('f_statistics')].strip()
                    # variance = row[col.index('variance')].strip()
                    # bayes_value = row[col.index('bayes_value')].strip()
                    # likelihood_ratio = row[col.index('likelihood_ratio')].strip()
                    trait_id = row[col.index('trait_id')].strip()
                    # dom_effect = row[col.index('dom_effect')].strip()
                    # add_effect = row[col.index('add_effect')].strip()
                    pubmed_id = row[col.index('pubmed_id')].strip()
                    gene_id = row[col.index('gene_id')].strip()
                    gene_id_src = row[col.index('gene_id_src')].strip()
                    # gene_id_type = row[col.index('gene_id_type')].strip()
                    # empty2 = row[col.index('empty2')].strip()

                if self.test_mode and int(qtl_id) not in self.test_ids:
                    continue

                qtl_id = common_name + 'QTL:' + qtl_id.strip()
                trait_id = ':'.join((aql_curie, trait_id.strip()))

                # Add QTL to graph
                feature = Feature(graph, qtl_id, qtl_symbol, self.globaltt['QTL'])
                feature.addTaxonToFeature(taxon_curie)

                # deal with the chromosome
                chrom_id = makeChromID(chromosome, taxon_curie, 'CHR')

                # add a version of the chromosome which is defined as
                # the genetic map
                build_id = 'MONARCH:' + common_name.strip() + '-linkage'
                build_label = common_name + ' genetic map'
                geno.addReferenceGenome(build_id, build_label, taxon_curie)
                chrom_in_build_id = makeChromID(chromosome, build_id, 'MONARCH')
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
                        (start, stop) = [
                            int(float(x.strip())) for x in re.split(r'-', range_cm)]
                    else:
                        LOG.info(
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
                    [self.globaltt['FuzzyPosition']])
                feature.addFeatureEndLocation(
                    stop, chrom_in_build_id, None,
                    [self.globaltt['FuzzyPosition']])
                feature.addFeatureToGraph()

                # sometimes there's a peak marker, like a rsid.
                # we want to add that as a variant of the gene,
                # and xref it to the qtl.
                dbsnp_id = None
                if peak_mark != '' and peak_mark != '.' and \
                        re.match(r'rs', peak_mark.strip()):
                    dbsnp_id = 'dbSNP:' + peak_mark.strip()

                    model.addIndividualToGraph(
                        dbsnp_id, None, self.globaltt['sequence_alteration']
                    )

                    model.addXref(
                        qtl_id, dbsnp_id, xref_category=blv.terms.SequenceVariant.value
                    )

                gene_id = gene_id.replace('uncharacterized ', '').strip()
                if gene_id is not None and gene_id != '' and gene_id != '.'\
                        and re.fullmatch(r'[^ ]*', gene_id) is not None:

                    # we assume if no src is provided and gene_id is an integer,
                    # then it is an NCBI gene ... (okay, lets crank that back a notch)
                    if gene_id_src == '' and gene_id.isdigit() and \
                            gene_id in self.gene_info:
                        # LOG.info(
                        #    'Warm & Fuzzy saying %s is a NCBI gene for %s',
                        #    gene_id, common_name)
                        gene_id_src = 'NCBIgene'
                    elif gene_id_src == '' and gene_id.isdigit():
                        LOG.warning(
                            'Cold & Prickely saying %s is a NCBI gene for %s',
                            gene_id, common_name)
                        gene_id_src = 'NCBIgene'
                    elif gene_id_src == '':
                        LOG.error(
                            ' "%s" is a NOT NCBI gene for %s', gene_id, common_name)
                        gene_id_src = None

                    if gene_id_src == 'NCBIgene':
                        gene_id = 'NCBIGene:' + gene_id
                        # we will expect that these will get labels elsewhere
                        geno.addGene(gene_id, None)
                        # FIXME what is the right relationship here?
                        geno.addAffectedLocus(qtl_id, gene_id)

                        if dbsnp_id is not None:
                            # add the rsid as a seq alt of the gene_id
                            vl_id = '_:' + re.sub(
                                r':', '', gene_id) + '-' + peak_mark.strip()
                            geno.addSequenceAlterationToVariantLocus(dbsnp_id, vl_id)
                            geno.addAffectedLocus(vl_id, gene_id)

                # add the trait
                model.addClassToGraph(
                    trait_id,
                    trait_name,
                    class_category=blv.terms.PhenotypicFeature.value
                )

                # Add publication
                reference = None
                if re.match(r'ISU.*', pubmed_id):
                    pub_id = 'AQTLPub:' + pubmed_id.strip()
                    reference = Reference(graph, pub_id)
                elif pubmed_id != '':
                    pub_id = 'PMID:' + pubmed_id.strip()
                    reference = Reference(
                        graph, pub_id, self.globaltt['journal article']
                    )

                if reference is not None:
                    reference.addRefToGraph()

                # make the association to the QTL
                assoc = G2PAssoc(
                    graph, self.name, qtl_id, trait_id, self.globaltt['is marker for']
                )
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
                    scr = re.sub(r'<', '', p_values)
                    scr = re.sub(r',', '.', scr)  # international notation
                    if scr.isnumeric():
                        score = float(scr)
                        assoc.set_score(score)  # todo add score type
                # TODO add LOD score?
                assoc.add_association_to_graph()

                # make the association to the dbsnp_id, if found
                if dbsnp_id is not None:
                    # make the association to the dbsnp_id
                    assoc = G2PAssoc(
                        graph, self.name, dbsnp_id, trait_id,
                        self.globaltt['is marker for'])
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
                        scr = re.sub(r'<', '', p_values)
                        scr = re.sub(r',', '.', scr)
                        if scr.isnumeric():
                            score = float(scr)
                            assoc.set_score(score)  # todo add score type
                    # TODO add LOD score?

                    assoc.add_association_to_graph()

                # off by one - the following actually gives us (limit + 1) records
                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

        LOG.info("Done with QTL genetic info")
        return

    def _process_qtls_genomic_location(
            self, raw, src_key, txid, build_id, build_label, common_name, limit=None):
        """
        This method

        Triples created:

        :param limit:
        :return:
        """
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        line_counter = 0
        geno = Genotype(graph)
        # assume that chrs get added to the genome elsewhere

        taxon_curie = 'NCBITaxon:' + txid
        eco_id = self.globaltt['quantitative trait analysis evidence']
        LOG.info("Processing QTL locations for %s from %s", taxon_curie, raw)
        with gzip.open(raw, 'rt', encoding='ISO-8859-1') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            # no header in GFF, so no header checking
            col = self.files[src_key]['columns']
            for row in reader:
                line_counter += 1
                if re.match(r'^#', ' '.join(row)):
                    continue

                if len(row) != len(col):
                    LOG.warning("Problem parsing in %s row %s\n"
                                "got %s cols but expected %s",
                                raw, row, len(row), len(col))
                    continue
                else:
                    # Doing this non-positional mapping for consistency, but I'm not
                    # sure we need to do this for GFF, since columns in GFF are probably
                    # not going to change anytime soon.
                    chromosome = row[col.index('SEQNAME')].strip()
                    # qtl_source = row[col.index('SOURCE')].strip()
                    # qtl_type = row[col.index('FEATURE')].strip()
                    start_bp = row[col.index('START')].strip()
                    stop_bp = row[col.index('END')].strip()
                    # score = row[col.index('SCORE')].strip()
                    strand = row[col.index('STRAND')].strip()
                    # frame = row[col.index('FRAME')].strip()
                    attr = row[col.index('ATTRIBUTE')].strip()

                example = '''
Chr.Z   Animal QTLdb    Production_QTL  33954873      34023581...
QTL_ID=2242;Name="Spleen percentage";Abbrev="SPLP";PUBMED_ID=17012160;trait_ID=2234;
trait="Spleen percentage";breed="leghorn";"FlankMarkers=ADL0022";VTO_name="spleen mass";
MO_name="spleen weight to body weight ratio";Map_Type="Linkage";Model="Mendelian";
Test_Base="Chromosome-wise";Significance="Significant";P-value="<0.05";F-Stat="5.52";
Variance="2.94";Dominance_Effect="-0.002";Additive_Effect="0.01
                '''
                str(example)
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
                for attributes in attr_items:
                    if not re.search(r'=', attributes):
                        # remove this attribute from the list
                        bad_attrs.add(attributes)

                attr_set = set(attr_items) - bad_attrs
                attribute_dict = dict(item.split("=") for item in attr_set)

                qtl_num = attribute_dict.get('QTL_ID')
                if self.test_mode and int(qtl_num) not in self.test_ids:
                    continue
                # make association between QTL and trait based on taxon

                qtl_id = common_name + 'QTL:' + str(qtl_num)
                model.addIndividualToGraph(qtl_id, None, self.globaltt['QTL'])
                geno.addTaxon(taxon_curie, qtl_id)

                #
                trait_id = 'AQTLTrait:' + attribute_dict.get('trait_ID')

                # if pub is in attributes, add it to the association
                pub_id = None
                if 'PUBMED_ID' in attribute_dict.keys():
                    pub_id = attribute_dict.get('PUBMED_ID')
                    if re.match(r'ISU.*', pub_id):
                        pub_id = 'AQTLPub:' + pub_id.strip()
                        reference = Reference(graph, pub_id)
                    else:
                        pub_id = 'PMID:' + pub_id.strip()
                        reference = Reference(
                            graph, pub_id, self.globaltt['journal article'])
                    reference.addRefToGraph()

                # Add QTL to graph
                assoc = G2PAssoc(
                    graph, self.name, qtl_id, trait_id,
                    self.globaltt['is marker for'])
                assoc.add_evidence(eco_id)
                assoc.add_source(pub_id)
                if 'P-value' in attribute_dict.keys():
                    scr = re.sub(r'<', '', attribute_dict.get('P-value'))
                    if ',' in scr:
                        scr = re.sub(r',', '.', scr)
                    if scr.isnumeric():
                        score = float(scr)
                        assoc.set_score(score)

                assoc.add_association_to_graph()
                # TODO make association to breed
                # (which means making QTL feature in Breed background)

                # get location of QTL
                chromosome = re.sub(r'Chr\.', '', chromosome)
                chrom_id = makeChromID(chromosome, taxon_curie, 'CHR')

                chrom_in_build_id = makeChromID(chromosome, build_id, 'MONARCH')
                geno.addChromosomeInstance(
                    chromosome, build_id, build_label, chrom_id)
                qtl_feature = Feature(graph, qtl_id, None, self.globaltt['QTL'])
                if start_bp == '':
                    start_bp = None
                qtl_feature.addFeatureStartLocation(
                    start_bp, chrom_in_build_id, strand,
                    [self.globaltt['FuzzyPosition']])
                if stop_bp == '':
                    stop_bp = None
                qtl_feature.addFeatureEndLocation(
                    stop_bp, chrom_in_build_id, strand,
                    [self.globaltt['FuzzyPosition']])
                qtl_feature.addTaxonToFeature(taxon_curie)
                qtl_feature.addFeatureToGraph()

                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        # LOG.warning("Bad attribute flags in this file")  # what does this even mean??
        LOG.info("Done with QTL genomic mappings for %s", taxon_curie)
        return

    def _process_trait_mappings(self, raw, src_key, limit=None):
        """
        This method mapps traits from/to ...

        Triples created:

        :param limit:
        :return:
        """
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        col = self.files[src_key]['columns']

        with open(raw, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            header = next(reader, None)
            self.check_fileheader(col, header, src_key)
            for row in reader:
                # need to skip the last line
                if len(row) != len(col):
                    LOG.info("skipping line %d: %s", reader.line_num, '\t'.join(row))
                    continue
                if limit is not None and reader.line_num > limit:
                    break
                vto_id = row[col.index('VT')].strip()
                pto_id = row[col.index('LPT')].strip()
                cmo_id = row[col.index('CMO')].strip()
                ato_column = row[col.index('ATO')].strip()
                # species = row[col.index('Species')].strip()
                # trait_class = row[col.index('Class')].strip()
                # trait_type = row[col.index('Type')].strip()
                # qtl_count = row[col.index('QTL_Count')].strip()

                ato_id = re.sub(
                    r'ATO #', 'AQTLTrait:', re.sub(
                        r'\].*', '', re.sub(r'\[', '', ato_column)))
                ato_id = ato_id.strip()

                ato_label = re.sub(r'.*\]\s*', '', ato_column)

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

        LOG.info("Done with trait mappings")

    def getTestSuite(self):
        import unittest
        from tests.test_animalqtl import AnimalQTLdbTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(AnimalQTLdbTestCase)

        return test_suite
