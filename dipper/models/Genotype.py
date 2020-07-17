import logging
import re
from dipper.models.Model import Model
from dipper.models.Family import Family
from dipper.graph.Graph import Graph
from dipper.models.GenomicFeature import makeChromID, makeChromLabel
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

__author__ = 'nlw'
LOG = logging.getLogger(__name__)


class Genotype():
    """
    These methods provide convenient methods to
    add items related to a genotype and it's parts to a supplied graph.
    They follow the patterns set out in
    GENO https://github.com/monarch-initiative/GENO-ontology.
    For specific sequence features,
    we use the GenomicFeature class to create them.

    """

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".format(graph))
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        return

    def addGenotype(
            self, genotype_id, genotype_label,
            genotype_type=None,
            genotype_description=None
    ):
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
            genotype_type = self.globaltt['intrinsic_genotype']

        self.model.addIndividualToGraph(
            genotype_id, genotype_label, genotype_type, genotype_description
        )

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
            allele_type = self.globaltt['allele']  # TODO is this a good idea?
        self.model.addIndividualToGraph(
            allele_id, allele_label, allele_type, allele_description
        )

    def addGene(
            self, gene_id, gene_label=None, gene_type=None, gene_description=None
    ):
        ''' genes are classes '''
        if gene_type is None:
            gene_type = self.globaltt['gene']
        self.model.addClassToGraph(
            gene_id, gene_label, gene_type, gene_description
        )

    def addConstruct(
            self, construct_id, construct_label, construct_type=None,
            construct_description=None,
            construct_category=None, construct_type_category=None):
        """
        :param construct_id:
        :param construct_label:
        :param construct_type:
        :param construct_description:
        :param construct_category: a biolink category CURIE for construct_id
        :param construct_type_category: a biolink category CURIE for construct_type
        :return:

        """
        # TODO add base type for construct
        # if (constrcut_type is None):
        #    construct_type=self.construct_base_type
        self.model.addIndividualToGraph(construct_id, construct_label,
                                        construct_type, construct_description,
                                        ind_category=construct_category,
                                        ind_type_category=construct_type_category)

    def addDerivesFrom(
            self, child_id, parent_id, child_category=None, parent_category=None
    ):
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
            child_id,
            self.globaltt['derives_from'],
            parent_id,
            subject_category=child_category,
            object_category=parent_category
        )

    def addSequenceDerivesFrom(self, child_id, parent_id,
                               child_category=None,
                               parent_category=None):
        self.graph.addTriple(
            child_id,
            self.globaltt['sequence_derives_from'],
            parent_id,
            subject_category=child_category,
            object_category=parent_category
        )

        return

    def addAlleleOfGene(self, allele_id, gene_id, rel_id=None):
        """
        We make the assumption here that if the relationship is not provided,
        it is a
        GENO:is_allele_of.

        Here, the allele should be a variant_locus, not a sequence alteration.
        :param allele_id:
        :param gene_id:
        :param rel_id:
        :return:

        """
        if rel_id is None:
            rel_id = self.globaltt["is_allele_of"]
        self.graph.addTriple(allele_id, rel_id, gene_id)

    def addAffectedLocus(
            self, allele_id, gene_id, rel_id=None):
        """
        We make the assumption here that if the relationship is not provided,
        it is a
        GENO:has_affected_feature.

        Here, the allele should be a variant_locus, not a sequence alteration.
        :param allele_id:
        :param gene_id:
        :param rel_id:
        :return:

        """
        if rel_id is None:
            rel_id = self.globaltt['has_affected_feature']
        self.graph.addTriple(allele_id, rel_id, gene_id)

    def addGeneProduct(
            self,
            sequence_id,
            product_id,
            product_label=None,
            product_type=None,
            sequence_category=None,
            product_category=None
    ):
        """
        Add gene/variant/allele has_gene_product relationship
        Can be used to either describe a gene to transcript relationship
        or gene to protein
        :param sequence_id:
        :param product_id:
        :param product_label:
        :param product_type:
        :param sequence_category: a biolink category CURIE for sequence_id [blv.terms.Gene].value
        :param product_category: a biolink category CURIE for product_id
        :return:

        """
        if product_label is not None and product_type is not None:
            self.model.addIndividualToGraph(
                product_id,
                product_label,
                product_type,
                ind_category=product_category
            )
        self.graph.addTriple(
            sequence_id,
            self.globaltt['has gene product'],
            product_id,
            subject_category=sequence_category,
            object_category=product_category
        )

    def addPolypeptide(
            self,
            polypeptide_id,
            polypeptide_label=None,
            transcript_id=None,
            polypeptide_type=None
    ):
        """
        :param polypeptide_id:
        :param polypeptide_label:
        :param polypeptide_type:
        :param transcript_id:
        :return:

        """
        if polypeptide_type is None:
            polypeptide_type = self.globaltt['polypeptide']
        self.model.addIndividualToGraph(
            polypeptide_id, polypeptide_label, polypeptide_type
        )
        if transcript_id is not None:
            self.graph.addTriple(
                transcript_id,
                self.globaltt['translates_to'],
                polypeptide_id
            )

    def addPartsToVSLC(
            self,
            vslc_id,
            allele1_id,
            allele2_id,
            zygosity_id=None,
            allele1_rel=None,
            allele2_rel=None
    ):
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

        if allele1_id is not None:
            self.addParts(allele1_id, vslc_id, allele1_rel)
        if allele2_id is not None and allele2_id.strip() != '':
            self.addParts(allele2_id, vslc_id, allele2_rel)

        # figure out zygosity if it's not supplied
        if zygosity_id is None:
            if allele1_id == allele2_id:
                zygosity_id = self.globaltt['homozygous']
            else:
                zygosity_id = self.globaltt['heterozygous']

        if zygosity_id is not None:
            self.graph.addTriple(vslc_id, self.globaltt['has_zygosity'], zygosity_id)


    def addVSLCtoParent(
            self,
            vslc_id,
            parent_id,
            part_category=None,
            parent_category=None
    ):
        """
        The VSLC can either be added to a genotype or to a GVC.
        The vslc is added as a part of the parent.
        :param vslc_id:
        :param parent_id:
        :param part_category: a biolink category CURIE for part
        :param parent_category: a biolink category CURIE for parent
        :return:
        """

        self.addParts(
            vslc_id,
            parent_id,
            self.globaltt['has_variant_part'],
            part_category=part_category,
            parent_category=parent_category
        )

    def addParts(
            self,
            part_id,
            parent_id,
            part_relationship=None,
            part_category=None,
            parent_category=None
    ):
        """
        This will add a has_part (or subproperty) relationship between
        a parent_id and the supplied part.
        By default the relationship will be BFO:has_part,
        but any relationship could be given here.
        :param part_id:
        :param parent_id:
        :param part_relationship:
        :param part_category: a biolink vocab curie for part_id
        :param parent_category: a biolink vocab curie for parent_id
        :return:

        """
        if part_relationship is None:
            part_relationship = self.globaltt['has_part']
        # Fail loudly if parent or child identifiers are None
        if parent_id is None:
            raise TypeError('Attempt to pass None as parent')
        elif part_id is None:
            raise TypeError('Attempt to pass None as child')
        elif part_relationship is None:
            part_relationship = self.globaltt['has_part']

        self.graph.addTriple(
            parent_id,
            part_relationship,
            part_id,
            subject_category=parent_category,
            object_category=part_category
        )

    def addSequenceAlteration(
            self,
            sa_id,
            sa_label,
            sa_type=None,
            sa_description=None
    ):

        if sa_type is None:
            sa_type = self.globaltt['sequence_alteration']

        self.model.addIndividualToGraph(
            sa_id,
            sa_label,
            sa_type,
            sa_description
        )

    def addSequenceAlterationToVariantLocus(self, sa_id, vl_id):
        self.addParts(sa_id, vl_id, self.globaltt['has_variant_part'])

    def addGenomicBackground(
            self,
            background_id,
            background_label,
            background_type=None,
            background_description=None
    ):
        if background_type is None:
            background_type = self.globaltt['genomic_background']
        self.model.addIndividualToGraph(
            background_id, background_label, background_type, background_description
        )

    def addGenomicBackgroundToGenotype(
            self, background_id, genotype_id, background_type=None
    ):
        if background_type is None:
            background_type = self.globaltt['genomic_background']
        self.model.addType(background_id, background_type)
        self.addParts(
            background_id, genotype_id, self.globaltt['has_reference_part']
        )

    def addTaxon(self, taxon_id, genopart_id, genopart_category=None):
        """
        The supplied geno part will have the specified taxon added with
        RO:in_taxon relation.
        Generally the taxon is associated with a genomic_background,
        but could be added to any genotype part (including a gene,
        regulatory element, or sequence alteration).
        :param taxon_id:
        :param genopart_id:
        :param genopart_category: a biolink term for genopart_id
        :return:

        """
        self.graph.addTriple(genopart_id, self.globaltt['in taxon'], taxon_id)

    def addGeneTargetingReagentToGenotype(self, reagent_id, genotype_id):

        """
        Add genotype has_variant_part reagent_id. For example, add a morphant
        reagent thingy to the genotype, assuming it's a extrinsic_genotype
        Also a triple to assign biolink categories to genotype and reagent.
        :param reagent_id
        :param genotype_id
        :return:

        """
        self.graph.addTriple(genotype_id, self.globaltt['has_variant_part'], reagent_id)

    def addGeneTargetingReagent(
            self,
            reagent_id,
            reagent_label,
            reagent_type,
            gene_id,
            description=None,
            reagent_category=None
    ):
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
            reagent_id,
            reagent_label,
            reagent_type,
            description,
            ind_category=reagent_category
        )

        self.graph.addTriple(reagent_id, self.globaltt['targets_gene'], gene_id)

    def addReagentTargetedGene(
            self,
            reagent_id,
            gene_id,
            targeted_gene_id=None,
            targeted_gene_label=None,
            description=None,
            reagent_category=None
    ):
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
        <reagent_id> GENO:targets_gene <gene_id>

        :param reagent_id:
        :param gene_id:
        :param targeted_gene_id:
        :param reagent_category: a biolink category CURIE for reagent_id
        :return:

        """

        # akin to a variant locus
        if targeted_gene_id is None:
            targeted_gene_id = '_' + gene_id + '-' + reagent_id
            targeted_gene_id = targeted_gene_id.replace(":", "")
        self.model.addIndividualToGraph(
            targeted_gene_id,
            targeted_gene_label,
            self.globaltt['reagent_targeted_gene'],
            description,
            ind_category=reagent_category
        )

        if gene_id is not None:
            self.graph.addTriple(
                targeted_gene_id, self.globaltt['is_expression_variant_of'], gene_id
            )

        self.graph.addTriple(
            targeted_gene_id, self.globaltt['is_targeted_by'], reagent_id
        )

    def addTargetedGeneSubregion(
            self, tgs_id, tgs_label, tgs_type=None, tgs_description=None):
        if tgs_type is None:
            tgs_type = self.globaltt['targeted_gene_subregion']

        self.model.addIndividualToGraph(tgs_id, tgs_label, tgs_type, tgs_description)

    def addMemberOfPopulation(self, member_id, population_id):
        self.graph.addTriple(
            population_id, self.globaltt['has_member_with_allelotype'], member_id
        )

    def addTargetedGeneComplement(
            self, tgc_id, tgc_label, tgc_type=None, tgc_description=None
    ):
        if tgc_type is None:
            tgc_type = self.globaltt['targeted_gene_complement']
        self.model.addIndividualToGraph(tgc_id, tgc_label, tgc_type, tgc_description)

    def addGenome(self, taxon_num, taxon_label=None, genome_id=None):
        ncbitaxon = 'NCBITaxon:' + taxon_num
        if taxon_label is None:
            if ncbitaxon in self.globaltcid:
                taxon_label = self.globaltcid[ncbitaxon]
            else:
                logging.warning('Add ' + ncbitaxon + ' to global translation table')
                taxon_label = taxon_num
        elif ncbitaxon in self.globaltcid and taxon_label != self.globaltcid[ncbitaxon]:
            logging.warning(
                '"' + self.globaltcid[ncbitaxon] + '" may need updating from "' +
                taxon_label + '" in global translation table')
            logging.warning(
                '"' + taxon_label + '": " ' + self.globaltcid[ncbitaxon] + '"' +
                ' may need to be added to a local translation table')

        genome_label = taxon_label + ' genome'

        if genome_id is None:
            genome_id = self.makeGenomeID(taxon_num)

        self.model.addClassToGraph(genome_id, genome_label, self.globaltt['genome'])

    def addReferenceGenome(self, build_id, build_label, taxon_id):
        genome_id = self.makeGenomeID(taxon_id)
        self.model.addIndividualToGraph(
            build_id,
            build_label,
            self.globaltt['reference_genome'],
            blv.terms['GenomeBuild']
        )
        self.model.addType(
            build_id, genome_id, subject_category=blv.terms['GenomeBuild']
        )
        if re.match(r'[0-9]+', taxon_id):
            taxon_id = 'NCBITaxon:' + taxon_id

        self.addTaxon(taxon_id, build_id, genopart_category=blv.terms['GenomeBuild'])

    @staticmethod
    def makeGenomeID(taxon_id):
        # scrub off the taxon prefix.  put it in base space
        # TODO: revisit as yet another BNODE?
        # should never be called if a real genome iri exists
        # should create the opaque bode and label together
        # genome_id = re.sub(r'.*\:', '_:', taxon_id) + 'genome'
        genome_id = '_:' + taxon_id + 'genome'
        return genome_id

    def addChromosome(
            self, chrom, tax_id, tax_label=None, build_id=None, build_label=None):
        """
        if it's just the chromosome, add it as an instance of a SO:chromosome,
        and add it to the genome. If a build is included,
        punn the chromosome as a subclass of SO:chromsome, and make the
        build-specific chromosome an instance of the supplied chr.
        The chr then becomes part of the build or genome.
        """
        family = Family(self.graph)
        # first, make the chromosome class, at the taxon level
        chr_id = makeChromID(str(chrom), tax_id)
        if tax_label is not None:
            chr_label = makeChromLabel(chrom, tax_label)
        else:
            chr_label = makeChromLabel(chrom)
        genome_id = self.makeGenomeID(tax_id)
        self.model.addClassToGraph(chr_id, chr_label, self.globaltt['chromosome'])
        self.addTaxon(tax_id, genome_id)  # add the taxon to the genome

        if build_id is not None:
            # the build-specific chromosome
            chrinbuild_id = makeChromID(chrom, build_id)
            if build_label is None:
                build_label = build_id
            chrinbuild_label = makeChromLabel(chrom, build_label)
            # add the build-specific chromosome as an instance of the chr class

            self.model.addIndividualToGraph(chrinbuild_id, chrinbuild_label, chr_id)

            # add the build-specific chromosome
            # as a member of the build (both ways)
            family.addMember(
                build_id, chrinbuild_id, group_category=blv.terms['GenomeBuild']
            )
            family.addMemberOf(
                chrinbuild_id, build_id, group_category=blv.terms['GenomeBuild']
            )

    def addChromosomeClass(self, chrom_num, taxon_id, taxon_label):
        taxon = re.sub('NCBITaxon:', '', taxon_id)
        # the chrom class (generic) id
        chrom_class_id = makeChromID(chrom_num, taxon, 'CHR')
        chrom_class_label = makeChromLabel(chrom_num, taxon_label)
        self.model.addClassToGraph(
            chrom_class_id, chrom_class_label, self.globaltt['chromosome']
        )

    def addChromosomeInstance(
            self, chr_num, reference_id, reference_label, chr_type=None
    ):
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
            chr_id, chr_label, self.globaltt['chromosome']
        )
        if chr_type is not None:
            self.model.addType(chr_id, chr_type)

        # add the build-specific chromosome
        # as a member of the build  (both ways)
        family.addMember(
            reference_id, chr_id, group_category=blv.terms['GenomeBuild']
        )
        family.addMemberOf(chr_id, reference_id)

    @staticmethod
    def make_variant_locus_label(gene_label, allele_label):
        if gene_label is None:
            gene_label = ''
        label = gene_label.strip() + '<' + allele_label.strip() + '>'

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

        if gene_label is None and allele1_label is None and allele2_label is None:
            LOG.error("Not enough info to make vslc label")
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
        animal_id = '_:' + animal_id

        animal_label = ' '.join((genotype_label, taxon_label))
        self.model.addIndividualToGraph(animal_id, animal_label, taxon_id)

        self.graph.addTriple(
            animal_id, self.globaltt['has_genotype'], genotype_id
        )
        return animal_id
