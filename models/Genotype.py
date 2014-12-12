__author__ = 'nicole'

from rdflib import BNode, Graph, Literal, RDF, OWL, extras, Namespace, URIRef
from rdflib.namespace import FOAF, RDFS, DC


from models.Assoc import Assoc
import re
import urllib

# The first of many associations
#This one is specific for making a disease-to-phenotype
class Genotype():
    #FIXME init might not require all these elements

    #TODO turn this into a proper dictionary
    #classes and relationships
    has_disposition = "GENO:0000208"
    intrinsic_genotype = "http://purl.obolibrary.org/obo/GENO_0000000"
    allele_base_type = 'http://purl.obolibrary.org/obo/GENO_0000008'
    gene_base_type = 'http://purl.obolibrary.org/obo/SO_0000704'
    derives_from = 'http://purl.obolibrary.org/obo/RO_0001000'
    is_mutant_of = URIRef('http://purl.obolibrary.org/obo/GENO_0000440')
    has_zygosity=URIRef('http://purl.obolibrary.org/obo/GENO_0000400')


    CLASS = "owl:Class"
    IND = "owl:NamedIndividual"
    monarch_namespace = Namespace('http://www.monarchinitiative.org/')

    prefixes = {'OMIM', 'http://www.omim.org/',
                'DECIPHER', 'http://www.decipher.org/',
                'ORPHANET', 'http://www.orpha.net/',
                'HP', 'http://purl.obolibrary.org/obo/'}

    def __init__(self, genotype_id, genotype_label, description=None):

        self.g = Graph()

        #create the genotype and add it to the graph
        #TODO this makes the assumption that it is an intrinsic_genotype; need to generalize
        self.geno = URIRef(genotype_id)
        self.addNode(genotype_id, genotype_label, self.intrinsic_genotype)

        return

    def addAllele(self, allele_id, allele_label, allele_type=None, allele_description=None):
        if (allele_type is None):
            allele_type = self.allele_base_type
        self.addNode(allele_id, allele_label, allele_type, allele_description)

        return

    def addGene(self, gene_id, gene_label, gene_type=None, gene_description=None):
        if (gene_type is None):
            gene_type = self.gene_base_type
        self.addNode(gene_id, gene_label, gene_type, gene_description)

        return

    def addConstruct(self, construct_id, construct_label, construct_type=None, construct_description=None):
        #todo add base type for construct
        #if (constrcut_type is None):
        #    constrcut_type=self.construct_base_type
        self.addNode(construct_id, construct_label, construct_type, construct_description)

        return

    def addAlleleDerivesFromConstruct(self, allele_id, construct_id):
        self.g.add((URIRef(allele_id), URIRef(self.derives_from), URIRef(construct_id)))
        return


    def addNode(self, id, label, type=None, description=None):
        n = URIRef(id)

        self.g.add((n, RDF['type'], OWL['class']))
        self.g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(type)
            self.g.add((n, OWL['subclassOf'], t))
        if (description is not None):
            self.g.add((n, DC['description'], Literal(description)))
        return

    def addAlleleOfGene(self, allele_id, gene_id, rel_id=None):
        rel=rel_id
        if (rel_id is None):
            rel=self.is_mutant_of
        self.g.add((URIRef(allele_id), URIRef(rel), URIRef(gene_id)))

        return

    def addAlleleToGenotype(self, allele_id, genotype_id):
        #TODO perhaps here is where we'll build the other genotype parts?
        #for now, just keep it simple and add the allele to the genotype atomically
        self.g.add((URIRef(genotype_id), OWL['hasPart'], URIRef(allele_id)))
        return

    def getGraph(self):
        return self.g