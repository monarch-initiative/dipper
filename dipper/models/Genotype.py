import logging
import re
from dipper.models.Model import Model
from dipper.models.Family import Family
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.graph.Graph import Graph

__author__ = 'nlw'
logger = logging.getLogger(__name__)


class Genotype():
    """
    These methods provide convenient methods to
    add items related to a genotype and it's parts to a supplied graph.
    They follow the patterns set out in
    GENO https://github.com/monarch-initiative/GENO-ontology.
    For specific sequence features,
    we use the GenomicFeature class to create them.

    """

    # special genotype parts mapped to their
    # GENO and SO classes that we explicitly reference here
    genoparts = {
        'intrinsic_genotype': 'GENO:0000000',
        'extrinsic_genotype': 'GENO:0000524',
        'effective_genotype': 'GENO:0000525',
        'sex_qualified_genotype': 'GENO:0000645',
        'male_genotype': 'GENO:0000646',
        'female_genotype': 'GENO:0000647',
        'genomic_background': 'GENO:0000611',
        'unspecified_genomic_background': 'GENO:0000649',
        'genomic_variation_complement': 'GENO:0000009',
        'karyotype_variation_complement': 'GENO:0000644',
        'variant_single_locus_complement': 'GENO:0000030',
        'variant_locus': 'GENO:0000002',
        'reference_locus': 'GENO:0000036',
        'allele': 'GENO:0000512',
        'gene': 'SO:0000704',
        'QTL': 'SO:0000771',
        'transgene': 'SO:0000902',  # not really used any more
        'transgenic_insertion': 'SO:0001218',
        'pseudogene': 'SO:0000336',
        'cytogenetic marker': 'SO:0000341',
        'sequence_feature': 'SO:0000110',
        'sequence_alteration': 'SO:0001059',
        'insertion': 'SO:0000667',
        'deletion': 'SO:0000159',
        'substitution': 'SO:1000002',
        'duplication': 'SO:1000035',
        'translocation': 'SO:0000199',
        'inversion': 'SO:1000036',
        'tandem_duplication': 'SO:1000173',
        'point_mutation': 'SO:1000008',
        'population': 'PCO:0000001',  # population
        'family': 'PCO:0000020',  # family
        'wildtype': 'GENO:0000511',
        'reagent_targeted_gene': 'GENO:0000504',
        'targeted_gene_subregion': 'GENO:0000534',
        'targeted_gene_complement': 'GENO:0000527',
        'biological_region': 'SO:0001411',
        'missense_variant': 'SO:0001583',
        'transcript': 'SO:0000233',
        'polypeptide': 'SO:0000104',
        'cDNA': 'SO:0000756',
        'sequence_variant_causing_loss_of_function_of_polypeptide':
            'SO:1000118',
        'sequence_variant_causing_gain_of_function_of_polypeptide':
            'SO:1000125',
        'sequence_variant_causing_inactive_catalytic_site': 'SO:1000120',
        'sequence_variant_affecting_polypeptide_function': 'SO:1000117',
        'regulatory_transgene_feature': 'GENO:0000637',
        'coding_transgene_feature': 'GENO:0000638',
        'protein_coding_gene': 'SO:0001217',
        'ncRNA_gene': 'SO:0001263',
        'RNAi_reagent': 'SO:0000337',
        'heritable_phenotypic_marker': 'SO:0001500'
    }

    object_properties = {
        'is_mutant_of': 'GENO:0000440',
        'derives_from': 'RO:0001000',
        'has_alternate_part': 'GENO:0000382',
        'has_reference_part': 'GENO:0000385',
        'has_sex_agnostic_genotype_part': 'GENO:0000650',
        'in_taxon': 'RO:0002162',
        'has_zygosity': 'GENO:0000608',
        # is_seq_var_inst_of links a alternate locus (instance)
        # to a gene (class)
        'is_sequence_variant_instance_of': 'GENO:0000408',
        'targets_instance_of': 'GENO:0000414',
        'is_reference_instance_of': 'GENO:0000610',
        'has_part': 'BFO:0000051',
        # use has_member_with_allelotype when relating populations
        'has_member_with_allelotype': 'GENO:0000225',
        'is_allelotype_of': 'GENO:0000206',
        'has_genotype': 'GENO:0000222',
        'has_phenotype': 'RO:0002200',
        'has_gene_product': 'RO:0002205',
        'translates_to': 'RO:0002513',
        'is_targeted_expression_variant_of': 'GENO:0000443',
        'is_transgene_variant_of': 'GENO:0000444',
        'has_variant_part': 'GENO:0000382',
        # targeted_by isa between a (reagent-targeted gene) and a morpholino
        'targeted_by': 'GENO:0000634',
        # FIXME should derives_sequence_from_gene just be subsequence of?
        'derives_sequence_from_gene': 'GENO:0000639',
        'has_affected_locus': 'GENO:0000418'
    }

    annotation_properties = {
        # TODO change properties with
        # https://github.com/monarch-initiative/GENO-ontology/issues/21
        # FIXME
        # reference_nucleotide, reference_amino_acid, altered_nucleotide
        # results_in_amino_acid_change are FIXME Made up terms
        'reference_nucleotide': 'GENO:reference_nucleotide',
        'reference_amino_acid': 'GENO:reference_amino_acid',
        'altered_nucleotide': 'GENO:altered_nucleotide',
        'results_in_amino_acid_change': 'GENO:results_in_amino_acid_change'
    }

    zygosity = {
        'homoplasmic': 'GENO:0000602',
        'heterozygous': 'GENO:0000135',
        'indeterminate': 'GENO:0000137',
        'heteroplasmic': 'GENO:0000603',
        'hemizygous-y': 'GENO:0000604',
        'hemizygous-x': 'GENO:0000605',
        'homozygous': 'GENO:0000136',
        'hemizygous': 'GENO:0000606',
        'complex_heterozygous': 'GENO:0000402',
        'simple_heterozygous': 'GENO:0000458'
    }

    properties = object_properties.copy()
    properties.update(annotation_properties)

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)

        return

    def addGenotype(
            self, genotype_id, genotype_label, genotype_type=None,
            genotype_description=None):
        """
        If a genotype_type is not supplied,
        we will default to 'intrinsic_genotype'
        :param genotype_id:
        :param genotype_label:
        :param genotype_type:
        :param genotype_description:
        :return:

        """

        if genotype_type is None:
            genotype_type = self.genoparts['intrinsic_genotype']

        self.model.addIndividualToGraph(
            genotype_id, genotype_label, genotype_type, genotype_description)
        return

    def addAllele(
            self, allele_id, allele_label, allele_type=None,
            allele_description=None):
        """
        Make an allele object.
        If no allele_type is added, it will default to a geno:allele
        :param allele_id: curie for allele (required)
        :param allele_label: label for allele (required)
        :param allele_type: id for an allele type (optional,
        recommended SO or GENO class)
        :param allele_description: a free-text description of the allele
        :return:

        """

        # TODO should we accept a list of allele types?
        if allele_type is None:
            allele_type = self.genoparts['allele']  # TODO is this a good idea?
        self.model.addIndividualToGraph(
            allele_id, allele_label, allele_type, allele_description)

        return

    def addGene(
            self, gene_id, gene_label, gene_type=None,
            gene_description=None):
        if gene_type is None:
            gene_type = self.genoparts['gene']
        # genes are classes
        self.model.addClassToGraph(
            gene_id, gene_label, gene_type, gene_description)

        return

    def addConstruct(self, construct_id, construct_label, construct_type=None,
                     construct_description=None):
        # TODO add base type for construct
        # if (constrcut_type is None):
        #    constrcut_type=self.construct_base_type
        self.model.addIndividualToGraph(construct_id, construct_label,
                                        construct_type, construct_description)

        return

    def addDerivesFrom(self, child_id, parent_id):
        """
        We add a derives_from relationship between the child and parent id.
        Examples of uses include between:
        an allele and a construct or strain here,
        a cell line and it's parent genotype.  Adding the parent and child to
        the graph should happen outside of this function call to ensure graph
        integrity.
        :param child_id:
        :param parent_id:
        :return:

        """

        self.graph.addTriple(
            child_id, self.properties['derives_from'], parent_id)

        return

    def addSequenceDerivesFrom(self, child_id, parent_id):
        self.graph.addTriple(
            child_id, self.properties['derives_sequence_from_gene'], parent_id)
        return

    def addAlleleOfGene(self, allele_id, gene_id, rel_id=None):
        """
        We make the assumption here that if the relationship is not provided,
        it is a
        GENO:is_sequence_variant_instance_of.

        Here, the allele should be a variant_locus, not a sequence alteration.
        :param allele_id:
        :param gene_id:
        :param rel_id:
        :return:

        """
        if rel_id is None:
            rel_id = self.properties['is_sequence_variant_instance_of']
        self.graph.addTriple(allele_id, rel_id, gene_id)
        return

    def addAffectedLocus(self, allele_id, gene_id, rel_id=None):
        """
        We make the assumption here that if the relationship is not provided,
        it is a
        GENO:is_sequence_variant_instance_of.

        Here, the allele should be a variant_locus, not a sequence alteration.
        :param allele_id:
        :param gene_id:
        :param rel_id:
        :return:

        """
        if rel_id is None:
            rel_id = self.properties['has_affected_locus']
        self.graph.addTriple(allele_id, rel_id, gene_id)
        return

    def addGeneProduct(
            self, sequence_id, product_id,
            product_label=None, product_type=None):
        """
        Add gene/variant/allele has_gene_product relationship
        Can be used to either describe a gene to transcript relationship
        or gene to protein
        :param sequence_id:
        :param product_id:
        :param product_label:
        :param product_type:
        :return:

        """
        if product_label is not None and product_type is not None:
            self.model.addIndividualToGraph(
                product_id, product_label, product_type)
        self.graph.addTriple(
            sequence_id, self.properties['has_gene_product'], product_id)

        return

    def addPolypeptide(
            self, polypeptide_id, polypeptide_label=None,
            transcript_id=None, polypeptide_type=None, ):
        """
        :param polypeptide_id:
        :param polypeptide_label:
        :param polypeptide_type:
        :param transcript_id:
        :return:

        """

        if polypeptide_type is None:
            polypeptide_type = self.genoparts['polypeptide']
        self.model.addIndividualToGraph(
            polypeptide_id, polypeptide_label, polypeptide_type)
        if transcript_id is not None:
            self.graph.addTriple(
                transcript_id, self.properties['translates_to'],
                polypeptide_id)

        return

    def addPartsToVSLC(
            self, vslc_id, allele1_id, allele2_id, zygosity_id=None,
            allele1_rel=None, allele2_rel=None):
        """
        Here we add the parts to the VSLC.  While traditionally alleles
        (reference or variant loci) are traditionally added, you can add any
        node (such as sequence_alterations for unlocated variations) to a vslc
        if they are known to be paired.  However, if a sequence_alteration's
        loci is unknown, it probably should be added directly to the GVC.
        :param vslc_id:
        :param allele1_id:
        :param allele2_id:
        :param zygosity_id:
        :param allele1_rel:
        :param allele2_rel:
        :return:

        """

        # vslc has parts allele1/allele2

        # vslc = gu.getNode(vslc_id)  # TODO unused
        if allele1_id is not None:
            self.addParts(allele1_id, vslc_id, allele1_rel)
        if allele2_id is not None and allele2_id.strip() != '':
            self.addParts(allele2_id, vslc_id, allele2_rel)

        # figure out zygosity if it's not supplied
        if zygosity_id is None:
            if allele1_id == allele2_id:
                zygosity_id = self.zygosity['homozygous']
            else:
                zygosity_id = self.zygosity['heterozygous']

        if zygosity_id is not None:
            self.graph.addTriple(
                vslc_id, self.properties['has_zygosity'], zygosity_id)

        return

    def addVSLCtoParent(self, vslc_id, parent_id):
        """
        The VSLC can either be added to a genotype or to a GVC.
        The vslc is added as a part of the parent.
        :param vslc_id:
        :param parent_id:
        :return:
        """

        self.addParts(
            vslc_id, parent_id, self.properties['has_alternate_part'])

        return

    def addParts(self, part_id, parent_id, part_relationship=None):
        """
        This will add a has_part (or subproperty) relationship between
        a parent_id and the supplied part.
        By default the relationship will be BFO:has_part,
        but any relationship could be given here.
        :param part_id:
        :param parent_id:
        :param part_relationship:
        :return:

        """

        if part_relationship is None:
            part_relationship = self.properties['has_part']

        self.graph.addTriple(parent_id, part_relationship, part_id)

        return

    def addSequenceAlteration(
            self, sa_id, sa_label, sa_type=None, sa_description=None):
        if sa_type is None:
            sa_type = self.genoparts['sequence_alteration']
        self.model.addIndividualToGraph(
            sa_id, sa_label, sa_type, sa_description)

        return

    def addSequenceAlterationToVariantLocus(self, sa_id, vl_id):
        self.addParts(sa_id, vl_id, self.properties['has_alternate_part'])
        return

    def addGenomicBackground(
            self, background_id, background_label,
            background_type=None, background_description=None):
        if background_type is None:
            background_type = self.genoparts['genomic_background']
        self.model.addIndividualToGraph(
            background_id, background_label, background_type,
            background_description)

        return

    def addGenomicBackgroundToGenotype(
            self, background_id, genotype_id, background_type=None):
        if background_type is None:
            background_type = self.genoparts['genomic_background']
        self.model.addType(background_id, background_type)
        self.addParts(background_id, genotype_id,
                      self.object_properties['has_reference_part'])

        return

    def addTaxon(self, taxon_id, genopart_id):
        """
        The supplied geno part will have the specified taxon added with
        RO:in_taxon relation.
        Generally the taxon is associated with a genomic_background,
        but could be added to any genotype part (including a gene,
        regulatory element, or sequence alteration).
        :param taxon_id:
        :param genopart_id:

        :return:

        """
        self.graph.addTriple(
            genopart_id, self.properties['in_taxon'], taxon_id)

        return

    def addGeneTargetingReagentToGenotype(self, reagent_id, genotype_id):
        # for example, add a morphant reagent thingy to the genotype,
        # assuming it's a extrinsic_genotype
        self.graph.addTriple(
            genotype_id, self.properties['has_variant_part'], reagent_id)

        return

    def addGeneTargetingReagent(
            self, reagent_id, reagent_label, reagent_type, gene_id,
            description=None):
        """
        Here, a gene-targeting reagent is added.
        The actual targets of this reagent should be added separately.
        :param reagent_id:
        :param reagent_label:
        :param reagent_type:

        :return:

        """

        # TODO add default type to reagent_type
        self.model.addIndividualToGraph(
            reagent_id, reagent_label, reagent_type, description)

        self.graph.addTriple(
            reagent_id, self.object_properties['targets_instance_of'], gene_id)

        return

    def addReagentTargetedGene(
            self, reagent_id, gene_id, targeted_gene_id=None,
            targeted_gene_label=None, description=None):
        """
        This will create the instance of a gene that is targeted by a molecular
        reagent (such as a morpholino or rnai).
        If an instance id is not supplied,
        we will create it as an anonymous individual which is of the type
        GENO:reagent_targeted_gene.
        We will also add the targets relationship between the reagent and
        gene class.

        <targeted_gene_id> a GENO:reagent_targeted_gene
        rdf:label targeted_gene_label
        dc:description description
        <reagent_id> GENO:targets_instance_of <gene_id>

        :param reagent_id:
        :param gene_id:
        :param targeted_gene_id:
        :return:

        """

        # akin to a variant locus
        if targeted_gene_id is None:
            targeted_gene_id = '_' + gene_id + '-' + reagent_id
            targeted_gene_id = targeted_gene_id.replace(":", "")
        self.model.addIndividualToGraph(
            targeted_gene_id, targeted_gene_label,

            self.genoparts['reagent_targeted_gene'], description)

        if gene_id is not None:
            self.graph.addTriple(
                targeted_gene_id,
                self.object_properties['is_targeted_expression_variant_of'],
                gene_id)

        self.graph.addTriple(
            targeted_gene_id, self.properties['targeted_by'], reagent_id)

        return

    def addTargetedGeneSubregion(
            self, tgs_id, tgs_label, tgs_type=None, tgs_description=None):
        if tgs_type is None:
            tgs_type = self.genoparts['targeted_gene_subregion']

        self.model.addIndividualToGraph(
            tgs_id, tgs_label, tgs_type, tgs_description)

    def addMemberOfPopulation(self, member_id, population_id):
        self.graph.addTriple(
            population_id,
            self.properties['has_member_with_allelotype'],
            member_id)
        return

    def addTargetedGeneComplement(
            self, tgc_id, tgc_label, tgc_type=None, tgc_description=None):
        if tgc_type is None:
            tgc_type = self.genoparts['targeted_gene_complement']
        self.model.addIndividualToGraph(
            tgc_id, tgc_label, tgc_type, tgc_description)

        return

    def addGenome(self, taxon_id, taxon_label=None):
        if taxon_label is None:
            taxon_label = taxon_id
        genome_label = taxon_label+' genome'
        genome_id = self.makeGenomeID(taxon_id)
        self.model.addClassToGraph(
            genome_id, genome_label, Feature.types['genome'])

        return

    def addReferenceGenome(self, build_id, build_label, taxon_id):
        genome_id = self.makeGenomeID(taxon_id)
        self.model.addIndividualToGraph(
            build_id, build_label, Feature.types['reference_genome'])
        self.model.addType(build_id, genome_id)
        self.addTaxon(taxon_id, build_id)

        return

    def makeGenomeID(self, taxon_id):
        # scrub off the taxon prefix.  put it in base space
        # TODO: revisit as BNODE?

        genome_id = re.sub(r'.*\:', ':', taxon_id) + 'genome'

        return genome_id

    def addChromosome(
            self, chr, tax_id, tax_label=None, build_id=None,
            build_label=None):
        """
        if it's just the chromosome, add it as an instance of a SO:chromosome,
        and add it to the genome. If a build is included,
        punn the chromosome as a subclass of SO:chromsome, and make the
        build-specific chromosome an instance of the supplied chr.
        The chr then becomes part of the build or genome.
        """
        family = Family()
        # first, make the chromosome class, at the taxon level
        chr_id = makeChromID(str(chr), tax_id)
        if tax_label is not None:
            chr_label = makeChromLabel(chr, tax_label)
        else:
            chr_label = makeChromLabel(chr)
        genome_id = self.makeGenomeID(tax_id)
        self.model.addClassToGraph(
            chr_id, chr_label, Feature.types['chromosome'])
        self.addTaxon(tax_id, genome_id)  # add the taxon to the genome

        if build_id is not None:
            # the build-specific chromosome
            chrinbuild_id = makeChromID(chr, build_id)
            if build_label is None:
                build_label = build_id
            chrinbuild_label = makeChromLabel(chr, build_label)
            # add the build-specific chromosome as an instance of the chr class

            self.model.addIndividualToGraph(
                chrinbuild_id, chrinbuild_label, chr_id)

            # add the build-specific chromosome
            # as a member of the build (both ways)
            family.addMember(build_id, chrinbuild_id)
            family.addMemberOf(chrinbuild_id, build_id)

        return

    def addChromosomeClass(self, chrom_num, taxon_id, taxon_label):
        taxon = re.sub('NCBITaxon:', '', taxon_id)
        # the chrom class (generic) id
        chrom_class_id = makeChromID(chrom_num, taxon, 'CHR')
        chrom_class_label = makeChromLabel(chrom_num, taxon_label)
        self.model.addClassToGraph(
            chrom_class_id, chrom_class_label, Feature.types['chromosome'])

        return

    def addChromosomeInstance(
            self, chr_num, reference_id, reference_label, chr_type=None):
        """
        Add the supplied chromosome as an instance within the given reference
        :param chr_num:
        :param reference_id: for example, a build id like UCSC:hg19
        :param reference_label:
        :param chr_type: this is the class that this is an instance of.
        typically a genome-specific chr

        :return:

        """
        family = Family(self.graph)
        chr_id = makeChromID(str(chr_num), reference_id, 'MONARCH')
        chr_label = makeChromLabel(str(chr_num), reference_label)

        self.model.addIndividualToGraph(
            chr_id, chr_label, Feature.types['chromosome'])
        if chr_type is not None:
            self.model.addType(chr_id, chr_type)

        # add the build-specific chromosome
        # as a member of the build  (both ways)
        family.addMember(reference_id, chr_id)
        family.addMemberOf(chr_id, reference_id)

        return

    def make_variant_locus_label(self, gene_label, allele_label):
        if gene_label is None:
            gene_label = ''
        label = gene_label.strip()+'<' + allele_label.strip() + '>'

        return label

    def make_vslc_label(self, gene_label, allele1_label, allele2_label):
        """
        Make a Variant Single Locus Complement (VSLC) in monarch-style.
        :param gene_label:
        :param allele1_label:
        :param allele2_label:
        :return:
        """

        vslc_label = ''

        if gene_label is None and \
                allele1_label is None and allele2_label is None:
            logger.error("Not enough info to make vslc label")
            return None

        top = self.make_variant_locus_label(gene_label, allele1_label)
        bottom = ''
        if allele2_label is not None:
            bottom = self.make_variant_locus_label(gene_label, allele2_label)

        vslc_label = '/'.join((top, bottom))

        return vslc_label

    def make_experimental_model_with_genotype(
             self, genotype_id, genotype_label, taxon_id, taxon_label):

        animal_id = '-'.join((taxon_id, 'with', genotype_id))
        animal_id = re.sub(r':', '', animal_id)
        animal_id = '_:'+animal_id

        animal_label = ' '.join((genotype_label, taxon_label))
        self.model.addIndividualToGraph(animal_id, animal_label, taxon_id)
        self.graph.addTriple(
            animal_id, Genotype.object_properties['has_genotype'], genotype_id)
        return animal_id
