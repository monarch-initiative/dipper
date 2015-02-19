__author__ = 'nicole'

from rdflib import Graph, RDF, OWL, URIRef

from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
from models.Assoc import Assoc
from conf import curie_map


class Genotype():
    '''
    This defines Genotype objects and their parts
    You should create one "Genotype" model/object per genotype that you are looping through,
    so that each genotype itself is a graph, then merge all the genotypes together
    to create a big merged graph for a single resource.
    '''

    #special genotype parts mapped to their GENO and SO classes that we explicitly reference here
    genoparts = {
        'intrinsic_genotype' : 'GENO:0000000',
        'effective_genotype' : 'GENO:0000525',
        'genomic_background' : 'GENO:0000010',
        'genomic_variation_complement' : 'GENO:0000009',
        'variant_single_locus_complement' : 'GENO:0000030',
        'alternate_locus' : 'GENO:0000512',
        'allele' : 'GENO:0000008',
        'gene' : 'SO:0000704',  #GENO:0000014 ?
        'QTL' : 'SO:0000771',
        'transgene' : 'SO:0000902',
        'pseudogene' : 'SO:0000336',
        'cytogenetic marker' : 'SO:0000341',  #chr band
        'sequence_feature' : 'SO:0000110',
        'sequence_alteration' : 'SO:0001059',
        'insertion' : 'SO:0000667',
        'deletion' : 'SO:0000159',
        'substitution' : 'SO:1000002',
        'duplication' : 'SO:1000035',
        'translocation' : 'SO:0000199',
        'inversion' : 'SO:1000036',
        'tandem_duplication' : 'SO:1000173',
        'point_mutation' : 'SO:1000008'
    }

    relationship = {
        'is_mutant_of' : 'GENO:0000440',
        'derives_from' : 'RO:0001000',
        'has_alternate_part' : 'GENO:0000382',
        'has_reference_part' : 'GENO:0000385',
        'in_taxon' : 'RO:0000216',
        'has_zygosity' : 'GENO:0000608',   #what exactly "has zygosity"?  is it the allele?  genotype?
        'is_sequence_variant_instance_of' : 'GENO:0000408',
    }

    zygosity = {
        'homoplasmic' : 'GENO:0000602',
        'heterozygous' : 'GENO:0000135',
        'indeterminate' : 'GENO:0000137',
        'heteroplasmic' : 'GENO:0000603',
        'hemizygous-y' : 'GENO:0000604',
        'hemizygous-x' : 'GENO:0000605',
        'homozygous' : 'GENO:0000136',
        'hemizygous' : 'GENO:0000606'
    }


    def __init__(self, genotype_id, genotype_label):

        self.cm = curie_map.get()

        self.gu = GraphUtils(self.cm)
        self.cu = CurieUtil(self.cm)

        self.g = Graph()

        #create the genotype and add it to the graph
        #TODO this makes the assumption that it is an intrinsic_genotype; need to generalize
        self.geno = URIRef(self.cu.get_uri(genotype_id))
        self.addNode(genotype_id, genotype_label, self.genoparts['intrinsic_genotype'])

        #add all the relationships as ObjectProperties into the graph
        for rel in self.relationship.keys():
            self.g.add((URIRef(self.cu.get_uri(self.relationship[rel])), RDF['type'], Assoc.OWLPROP))

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
        if (allele_type is None):
            allele_type = self.genoparts['allele']
        self.addNode(allele_id, allele_label, allele_type, allele_description)

        return

    def addGene(self, gene_id, gene_label, gene_type=None, gene_description=None):
        if (gene_type is None):
            gene_type = self.genoparts['gene']
        self.addNode(gene_id, gene_label, gene_type, gene_description)

        return

    def addConstruct(self, construct_id, construct_label, construct_type=None, construct_description=None):
        #todo add base type for construct
        #if (constrcut_type is None):
        #    constrcut_type=self.construct_base_type
        self.addNode(construct_id, construct_label, construct_type, construct_description)

        return

    def addAlleleDerivesFromConstruct(self, allele_id, construct_id):
        rel = self.cu.get_uri(self.relationship['derives_from'])
        self.g.add((URIRef(self.cu.get_uri(allele_id)), URIRef(rel), URIRef(self.cu.get_uri(construct_id))))
        return


    def addNode(self, id, label, type=None, description=None):
        self.gu.addClassToGraph(self.g,id,label,type,description)
        return

    def addAlleleOfGene(self, allele_id, gene_id, rel_id=None):
        rel=rel_id
        if (rel_id is None):
            rel=self.relationship['is_mutant_of']
        self.g.add((URIRef(self.cu.get_uri(allele_id)), URIRef(self.cu.get_uri(rel)), URIRef(self.cu.get_uri(gene_id))))

        return

    def addAlleleToGenotype(self, allele_id, genotype_id):
        #TODO perhaps here is where we'll build the other genotype parts?
        #for now, just keep it simple and add the allele to the genotype atomically
        self.g.add((URIRef(self.cu.get_uri(allele_id)), URIRef(OWL['hasPart']), URIRef(self.cu.get_uri(genotype_id))))
        return

    def getGraph(self):
        return self.g