__author__ = 'nlw'

from rdflib import RDF
import logging

from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper.models.GenomicFeature import Feature,makeChromID,makeChromLabel
import re

logger = logging.getLogger(__name__)

class Genotype():
    """
    These methods provide convenient methods to add items related to a genotype and it's parts to a supplied graph.
    They follow the patterns set out in GENO https://github.com/monarch-initiative/GENO-ontology.
    For specific sequence features, we use the GenomicFeature class to create them.
    """

    # special genotype parts mapped to their GENO and SO classes that we explicitly reference here
    genoparts = {
        'intrinsic_genotype': 'GENO:0000000',
        'extrinsic_genotype': 'GENO:0000524',
        'effective_genotype': 'GENO:0000525',
        'genomic_background': 'GENO:0000611',
        'genomic_variation_complement': 'GENO:0000009',
        'variant_single_locus_complement': 'GENO:0000030',
        'variant_locus': 'GENO:0000002',
        'reference_locus' : 'GENO:0000036',
        'allele': 'GENO:0000008',
        'gene': 'SO:0000704',
        'QTL': 'SO:0000771',
        'transgene': 'SO:0000902',
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
        'population' : 'GENO:0000110',  #collection of organisms; consider OBI:population or EFO:population
        'wildtype' : 'GENO:0000511',
        'reagent_targeted_gene': 'GENO:0000504',
        'targeted_gene_subregion' : 'GENO:0000534',
        'targeted_gene_complement' : 'GENO:0000527',
        'biological_region' : 'SO:0001411',
        'missense_variant': 'SO:0001583',
        'transcript': 'SO:0000233',
        'polypeptide': 'SO:0000104',
        'cDNA': 'SO:0000756',
        'sequence_variant_causing_loss_of_function_of_polypeptide': 'SO:1000118',
        'sequence_variant_causing_gain_of_function_of_polypeptide': 'SO:1000125',
        'sequence_variant_causing_inactive_catalytic_site': 'SO:1000120',
        'sequence_variant_affecting_polypeptide_function': 'SO:1000117'
    }

    object_properties = {
        'is_mutant_of': 'GENO:0000440',
        'derives_from': 'RO:0001000',
        'has_alternate_part': 'GENO:0000382',
        'has_reference_part': 'GENO:0000385',
        'in_taxon': 'RO:0002162',
        'has_zygosity': 'GENO:0000608',
        'is_sequence_variant_instance_of': 'GENO:0000408',  #links a alternate locus (instance) to a gene (class)
        'targets_instance_of': 'GENO:0000414',
        'is_reference_instance_of': 'GENO:0000610',
        'has_part': 'BFO:0000051',
        'has_member_with_allelotype': 'GENO:0000225',  #use this when relating populations
        'is_allelotype_of': 'GENO:0000206',
        'has_genotype': 'GENO:0000222',
        'has_phenotype': 'RO:0002200',
        'transcribed_to': 'RO:0002205',
        'translates_to': 'RO:0002513'
    }

    annotation_properties = {
        'reference_nucleotide': 'GENO:reference_nucleotide', #Made up term
        'reference_amino_acid': 'GENO:reference_amino_acid', #Made up term
        'altered_nucleotide': 'GENO:altered_nucleotide', #Made up term
        'results_in_amino_acid_change': 'GENO:results_in_amino_acid_change' #Made up term
    }

    zygosity = {
        'homoplasmic': 'GENO:0000602',
        'heterozygous': 'GENO:0000135',
        'indeterminate': 'GENO:0000137',
        'heteroplasmic': 'GENO:0000603',
        'hemizygous-y': 'GENO:0000604',
        'hemizygous-x': 'GENO:0000605',
        'homozygous': 'GENO:0000136',
        'hemizygous': 'GENO:0000606'
    }

    properties = object_properties.copy()
    properties.update(annotation_properties)

    def __init__(self, graph):

        self.gu = GraphUtils(curie_map.get())


        self.graph = graph

        self.gu.loadProperties(self.graph,self.object_properties,self.gu.OBJPROP)

        return

    def addGenotype(self, genotype_id, genotype_label, genotype_type=None, genotype_description=None):
        """
        If a genotype_type is not supplied, we will default to 'intrinsic_genotype'
        :param graph:
        :param genotype_id:
        :param genotype_label:
        :param genotype_type:
        :return:
        """
        if genotype_type is None:
            genotype_type = self.genoparts['intrinsic_genotype']

        self.gu.addIndividualToGraph(self.graph, genotype_id, genotype_label, genotype_type, genotype_description)

        return

    def addAllele(self, allele_id, allele_label, allele_type=None, allele_description=None):
        '''
        Make an allele object. If no allele_type is added, it will default to a geno:allele
        :param allele_id: curie for allele (required)
        :param allele_label: label for allele (required)
        :param allele_type: id for an allele type (optional, recommended SO or GENO class)
        :param allele_description: a free-text description of the allele
        :return:
        '''
        #TODO should we accept a list of allele types?
        if (allele_type is None):
            allele_type = self.genoparts['allele']  #TODO is this a good idea?
        self.gu.addIndividualToGraph(self.graph, allele_id, allele_label, allele_type, allele_description)

        return

    def addGene(self, gene_id, gene_label, gene_type=None, gene_description=None):
        if (gene_type is None):
            gene_type = self.genoparts['gene']
        #genes are classes
        self.gu.addClassToGraph(self.graph, gene_id, gene_label, gene_type, gene_description)

        return

    def addConstruct(self, construct_id, construct_label, construct_type=None, construct_description=None):
        #todo add base type for construct
        #if (constrcut_type is None):
        #    constrcut_type=self.construct_base_type
        self.gu.addIndividualToGraph(self.graph, construct_id, construct_label, construct_type, construct_description)

        return

    def addDerivesFrom(self, child_id, parent_id):
        """
        We add a derives_from relationship between the child and parent id.  Examples of uses include between:
        an allele and a construct or strain here, a cell line and it's parent genotype.  Adding the
        parent and child to the graph should happen outside of this function call to
        ensure graph integrity.
        :param allele_id:
        :param construct_id:
        :return:
        """
        rel = self.gu.getNode(self.properties['derives_from'])
        self.graph.add((self.gu.getNode(child_id), rel, self.gu.getNode(parent_id)))

        return


    def addAlleleOfGene(self, allele_id, gene_id, rel_id=None):
        """
        We make the assumption here that if the relationship is not provided, it is a
        GENO:is_sequence_variant_instance_of
        :param allele_id:
        :param gene_id:
        :param rel_id:
        :return:
        """
        if (rel_id is None):
            rel_id = self.properties['is_sequence_variant_instance_of']
        self.graph.add((self.gu.getNode(allele_id),self.gu.getNode(rel_id),self.gu.getNode(gene_id)))

        return

    def addTranscript(self, variant_id, transcript_id, transcript_label=None, transcript_type=None):
        """
        Add gene/variant/allele transcribes_to relationship
        :param graph:
        :param variant_id:
        :param transcript_id:
        :param transcript_label:
        :return:
        """
        self.gu.addIndividualToGraph(self.graph, transcript_id, transcript_label, transcript_type)
        self.gu.addTriple(self.graph, variant_id, self.properties['transcribed_to'], transcript_id)

        return

    def addPolypeptide(self, polypeptide_id, polypeptide_label=None, transcript_id=None, polypeptide_type=None, ):
        """
        :param polypeptide_id:
        :param polypeptide_label:
        :param polypeptide_type:
        :param transcript_id:
        :return:
        """
        if polypeptide_type is None:
            polypeptide_type = self.genoparts['polypeptide']
        self.gu.addIndividualToGraph(self.graph, polypeptide_id, polypeptide_label, polypeptide_type)
        if transcript_id is not None:
            self.gu.addTriple(self.graph, transcript_id, self.properties['translates_to'], polypeptide_id)

        return


    def addPartsToVSLC(self, vslc_id, allele1_id, allele2_id, zygosity_id=None):
        """
        Here we add the parts to the VSLC.  While traditionally alleles (reference or variant loci) are
        traditionally added, you can add any node (such as sequence_alterations for unlocated variations)
        to a vslc if they are known to be paired.  However, if a sequence_alteration's loci is unknown,
        it probably should be added directly to the GVC.
        :param vslc_id:
        :param allele1_id:
        :param allele2_id:
        :param zygosity_id:
        :return:
        """

        #vslc has parts allele1/allele2
        gu = self.gu
        has_zygosity = gu.getNode(self.properties['has_zygosity'])

        vslc = gu.getNode(vslc_id)
        if allele1_id is not None:
            self.addParts(allele1_id, vslc_id)
        if allele2_id is not None:
            self.addParts(allele2_id, vslc_id)

        #figure out zygosity if it's not supplied
        if (zygosity_id is None):
            if (allele1_id == allele2_id):
                zygosity_id = self.zygosity['homozygous']
            else:
                zygosity_id = self.zygosity['heterozygous']

        if (zygosity_id is not None):
            self.graph.add((vslc, has_zygosity, gu.getNode(zygosity_id)))

        return

    def addVSLCtoParent(self, vslc_id, parent_id):
        """
        The VSLC can either be added to a genotype or to a GVC.  The vslc is added as a part of the parent.
        :param vslc_id:
        :param parent_id:
        :return:
        """
        self.addParts(vslc_id, parent_id, self.properties['has_alternate_part'])

        return

    def addParts(self, part_id, parent_id, part_relationship=None):
        """
        This will add a has_part (or subproperty) relationship between a parent_id and the supplied part.
        By default the relationship will be BFO:has_part, but any relationship could be given here.
        :param variant_part:
        :param parent_id:
        :param part_relationship:
        :return:
        """
        gu = self.gu
        if part_relationship is None:
            has_part = gu.getNode(self.properties['has_part'])
        else:
            has_part = gu.getNode(part_relationship)

        self.graph.add((gu.getNode(parent_id), has_part, gu.getNode(part_id)))

        return

    def addSequenceAlteration(self, sa_id, sa_label, sa_type=None, sa_description=None):
        if sa_type is None:
            sa_type = self.genoparts['sequence_alteration']
        self.gu.addIndividualToGraph(self.graph, sa_id, sa_label, sa_type, sa_description)

        return

    def addSequenceAlterationToVariantLocus(self, sa_id, vl_id):
        self.addParts(sa_id, vl_id, self.properties['has_alternate_part'])
        return

    def addGenomicBackground(self, background_id, background_label, background_type=None, background_description=None):
        if background_type is None:
            background_type = self.genoparts['genomic_background']
        self.gu.addIndividualToGraph(self.graph, background_id, background_label, background_type, background_description)

        return


    def addGenomicBackgroundToGenotype(self, background_id, genotype_id):
        gu = self.gu

        self.graph.add((gu.getNode(background_id), RDF['type'], gu.getNode(self.genoparts['genomic_background'])))
        self.graph.add((gu.getNode(genotype_id), gu.getNode(self.properties['has_reference_part']), gu.getNode(background_id)))

        return

    def addTaxon(self, taxon_id, genopart_id):
        """
        The supplied geno part will have the specified taxon added with RO:in_taxon relation.
        Generally the taxon is associated with a genomic_background, but could be added to any
        genotype part (including a gene, regulatory element, or sequence alteration).
        :param taxon_id:
        :param genopart_id:
        :return:
        """
        in_taxon = self.gu.getNode(self.properties['in_taxon'])
        s = self.gu.getNode(genopart_id)
        self.graph.add((s, in_taxon, self.gu.getNode(taxon_id)))

        return

    def addGeneTargetingReagentToGenotype(self, reagent_id, genotype_id):
        #for example, add a morphant reagent thingy to the genotype, assuming it's a extrinsic_genotype


        return

    def addGeneTargetingReagent(self, reagent_id, reagent_label, reagent_type, description=None):
        """
        Here, a gene-targeting reagent is added.  The actual targets of this reagent should be added separately.
        :param reagent_id:
        :param reagent_label:
        :param reagent_type:
        :return:
        """
        #TODO add default type to reagent_type
        self.gu.addIndividualToGraph(self.graph, reagent_id, reagent_label, reagent_type, description)

        return

    def addReagentTargetedGene(self, reagent_id, gene_id, targeted_gene_id=None, targeted_gene_label=None,
                               description=None):
        """
        This will create the instance of a gene that is targeted by a molecular reagent (such as a morpholino or rnai).
        If an instance id is not supplied, we will create it as an anonymous individual which is of the
        type GENO:reagent_targeted_gene.  We will also add the targets relationship between the reagent and gene class.

        <targeted_gene_id> a GENO:reagent_targeted_gene
            rdf:label targeted_gene_label
            dc:description description
        <reagent_id> GENO:targets_instance_of <gene_id>

        :param reagent_id:
        :param gene_id:
        :param targeted_gene_id:
        :return:
        """

        #TODO is this a bad to assume the reagent is targeting at GENE specifically?
        # we are assuming that the reagent targets a GENE as opposed to any genomic feature.
        # but maybe that's not right, because i bet some reagents might target sequence_alterations.
        gu = self.gu
        targets = gu.getNode(self.properties['targets_instance_of'])
        self.graph.add((gu.getNode(reagent_id), targets, gu.getNode(gene_id)))

        if (targeted_gene_id is None):
            targeted_gene_id = '_' + gene_id + '-' + reagent_id
        self.gu.addIndividualToGraph(self.graph, targeted_gene_id, targeted_gene_label,
                                     self.genoparts['reagent_targeted_gene'], description)

        return

    def addTargetedGeneSubregion(self, tgs_id, tgs_label, tgs_type=None, tgs_description=None):
        if tgs_type is None:
            tgs_type = self.genoparts['targeted_gene_subregion']
        self.gu.addIndividualToGraph(self.graph, tgs_id, tgs_label, tgs_type, tgs_description)

        return


    def addTargetedGeneComplement(self, tgc_id, tgc_label, tgc_type=None, tgc_description=None):
        if tgc_type is None:
            tgc_type = self.genoparts['targeted_gene_complement']
        self.gu.addIndividualToGraph(self.graph, tgc_id, tgc_label, tgc_type, tgc_description)

        return

    def addMemberOfPopulation(self,member_id,population_id):
        self.graph.add((self.gu.getNode(population_id),
                        self.gu.getNode(self.properties['has_member_with_allelotype']),
                        self.gu.getNode(member_id)))

        return

    def addGenome(self,taxon_id,taxon_label=None):
        if taxon_label is None:
            taxon_label = taxon_id
        genome_label = taxon_label+' genome'
        genome_id = self.makeGenomeID(taxon_id)
        self.gu.addClassToGraph(self.graph,genome_id,genome_label,Feature.types['genome'])
        self.addTaxon(taxon_id,genome_id)

        return

    def addReferenceGenome(self,build_id,build_label,taxon_id):
        genome_id = self.makeGenomeID(taxon_id)
        self.gu.addIndividualToGraph(self.graph,build_id,build_label,Feature.types['reference_genome'])
        self.gu.addType(self.graph,build_id,genome_id)
        self.gu.addSubclass(self.graph,genome_id,build_id)

        return

    def makeGenomeID(self,taxon_id):
        #scrub off the taxon prefix.  put it in base space

        genome_id = re.sub('.*\:',':',taxon_id)+'genome'

        return genome_id

    def addChromosome(self,chr,tax_id,tax_label=None,build_id=None,build_label=None):
        #if it's just the chromosome, add it as an instance of a SO:chromosome, and add it to the genome.
        # if a build is included, punn the chromosome as a subclass of SO:chromosome, and
        # make the build-specific chromosome an instance of the supplied chr.  The chr then becomes part of the
        # build or genome.

        #first, make the chromosome class, at the taxon level
        chr_id = makeChromID(str(chr),tax_id)
        if tax_label is not None:
            chr_label = makeChromLabel(chr,tax_label)
        else:
            chr_label = makeChromLabel(chr)
        genome_id = self.makeGenomeID(tax_id)
        self.gu.addClassToGraph(self.graph,chr_id,chr_label,Feature.types['chromosome'])
        #add it as a member of the genome (both ways)
        self.gu.addMember(self.graph,genome_id,chr_id)
        self.gu.addMemberOf(self.graph,chr_id,genome_id)

        self.addTaxon(tax_id,genome_id)  #add the taxon to the genome


        if build_id is not None:
            chrinbuild_id = makeChromID(chr,build_id)  #the build-specific chromosome
            if build_label is None:
                build_label = build_id
            chrinbuild_label = makeChromLabel(chr,build_label)
            #add the build-specific chromosome as an instance of the chr class
            self.gu.addIndividualToGraph(self.graph,chrinbuild_id,chrinbuild_label,chr_id)

            #add the build-specific chromosome as a member of the build  (both ways)
            self.gu.addMember(self.graph,build_id,chrinbuild_id)
            self.gu.addMemberOf(self.graph,chrinbuild_id,build_id)

        return
