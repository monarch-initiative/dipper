__author__ = 'nlw'

from rdflib import Graph, Literal, RDF, OWL, URIRef
from rdflib.namespace import RDFS, DC

from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
from models.Assoc import Assoc
import curie_map
from models.GenomicFeature import Feature


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
        'point_mutation': 'SO:1000008'
    }

    relationship = {
        'is_mutant_of': 'GENO:0000440',
        'derives_from': 'RO:0001000',
        'has_alternate_part': 'GENO:0000382',
        'has_reference_part': 'GENO:0000385',
        'in_taxon' : 'RO:0002162',
        'has_zygosity': 'GENO:0000608',
        'is_sequence_variant_instance_of': 'GENO:0000408',
        'targets_instance_of': 'GENO:0000414',
        'is_reference_instance_of' : 'GENO:0000610',
        'has_part' : 'BFO:0000051'
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


    def __init__(self, graph):

        self.gu = GraphUtils(curie_map.get())

        self.graph = graph

        self.gu.loadObjectProperties(self.graph,self.relationship)

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

    def addAlleleDerivesFrom(self, allele_id, construct_or_strain_id):
        """
        We add a derives_from relationship between an allele and a construct or strain here.  Adding the
        allele and constructs to the graph should happen before (or after) this function call to
        ensure graph integrity.
        :param allele_id:
        :param construct_id:
        :return:
        """
        rel = self.gu.getNode(self.relationship['derives_from'])
        self.graph.add((self.gu.getNode(allele_id), rel, self.gu.getNode(construct_or_strain_id)))

        return


    def addAlleleOfGene(self, allele_id, gene_id, rel_id=None):
        """
        We make the assumption here that if the relationship is not provided, it is an alteration (not a reference) part
        :param allele_id:
        :param gene_id:
        :param rel_id:
        :return:
        """
        if (rel_id is None):
            rel_id = self.relationship['has_alternate_part']
        self.addParts(allele_id, gene_id, rel_id)

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
        has_zygosity = gu.getNode(self.relationship['has_zygosity'])

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
        self.addParts(vslc_id, parent_id, self.relationship['has_alternate_part'])

        return

    def addParts(self, part_id, parent_id, part_relationship=None):
        """
        This will add a has_part (or subproperty) relationship between a parent_id and the supplied part.
        By default the relationship will be GENO:has_part, but any relationship could be given here.
        :param variant_part:
        :param parent_id:
        :param part_relationship:
        :return:
        """
        gu = self.gu
        if part_relationship is None:
            has_part = gu.getNode(self.relationship['has_part'])
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
        self.addParts(sa_id, vl_id, self.relationship['has_alternate_part'])
        return

    def addGenomicBackgroundToGenotype(self, background_id, genotype_id):
        gu = self.gu

        self.graph.add((gu.getNode(background_id), RDF['type'], gu.getNode(self.genoparts['genomic_background'])))
        self.graph.add((gu.getNode(genotype_id), gu.getNode(self.relationship['has_reference_part']), gu.getNode(background_id)))

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
        in_taxon = self.gu.getNode(self.relationship['in_taxon'])
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
        targets = gu.getNode(self.relationship['targets_instance_of'])
        self.graph.add((gu.getNode(reagent_id), targets, gu.getNode(gene_id)))

        if (targeted_gene_id is None):
            targeted_gene_id = '_' + gene_id + '-' + reagent_id
        self.gu.addIndividualToGraph(self.graph, targeted_gene_id, targeted_gene_label,
                                     self.genoparts['reagent_targeted_gene'], description)

        return


