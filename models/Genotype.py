__author__ = 'nicole'

from rdflib import BNode, Graph, Literal, RDF, OWL, extras, Namespace, URIRef
from rdflib.namespace import FOAF, RDFS, DC


from models.Assoc import Assoc
import re
import urllib
from utils.CurieUtil import CurieUtil
from models.Assoc import Assoc

# The first of many associations
#This one is specific for making a disease-to-phenotype
class Genotype():
    #FIXME init might not require all these elements

    #TODO turn this into a proper dictionary
    #classes and relationships
    curie_map = {
        'GENO': 'http://purl.obolibrary.org/obo/GENO_',
        'SO' : 'http://purl.obolibrary.org/obo/SO_'
    }

    genoparts = {
        'intrinsic_genotype' : 'GENO:0000000',
        'allele_base_type' : 'GENO:0000008',
        'gene_base_type' : 'SO:0000704'
    }

    #relationships
    relationship = {
        'is_mutant_of' : 'GENO:0000440',
        'derives_from' : 'RO:0001000'
    }


    def __init__(self, genotype_id, genotype_label, curie_map):

        self.curie_map.update(curie_map)
        self.cu = CurieUtil(self.curie_map)

        self.g = Graph()

        #create the genotype and add it to the graph
        #TODO this makes the assumption that it is an intrinsic_genotype; need to generalize
        self.geno = URIRef(self.cu.get_uri(genotype_id))
        self.addNode(genotype_id, genotype_label, self.genoparts['intrinsic_genotype'])

        for rel in self.relationship.keys():
            self.g.add((URIRef(self.cu.get_uri(self.relationship[rel])), RDF['type'], Assoc.OWLPROP))

        return

    def addAllele(self, allele_id, allele_label, allele_type=None, allele_description=None):
        if (allele_type is None):
            allele_type = self.genoparts['allele_base_type']
        self.addNode(allele_id, allele_label, allele_type, allele_description)

        return

    def addGene(self, gene_id, gene_label, gene_type=None, gene_description=None):
        if (gene_type is None):
            gene_type = self.genoparts['gene_base_type']
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
        n = URIRef(self.cu.get_uri(id))

        self.g.add((n, RDF['type'], Assoc.OWLCLASS))
        self.g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(self.cu.get_uri(type))
            self.g.add((n, Assoc.OWLSUBCLASS, t))
        if (description is not None):
            self.g.add((n, DC['description'], Literal(description)))
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
        self.g.add((URIRef(self.cu.get_uri(genotype_id)), OWL['hasPart'], URIRef(self.cu.get_uri(allele_id))))
        return

    def getGraph(self):
        return self.g