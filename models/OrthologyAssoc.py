__author__ = 'nicole'

from rdflib import URIRef
from rdflib.namespace import RDF

from utils.CurieUtil import CurieUtil
from models.Assoc import Assoc
from conf import curie_map


class OrthologyAssoc(Assoc):

    relationships = {
        'orthologous' : 'RO:HOM0000017', #in orthology relationship with
        'least_diverged_orthologous' : 'RO:HOM0000020', #in 1 to 1 orthology relationship with
        'homologous': 'RO:HOM0000019', ## in 1 to 1 homology relationship with
        'paralogous' : 'RO:HOM0000011', ## in paralogy relationship with (generic)
        'in_paralogous' : 'RO:HOM0000023', #in in-paralogy relationship with
        'ohnologous' : 'RO:HOM0000022', #in ohnology relationship with
        'xenologous': 'RO:HOM0000018', ## in xenology relationship with
        'has_member' : 'RO:0002351'
        }


    def __init__(self,assoc_id, gene1, gene2, pub, evidence_code):
        self.cu = CurieUtil(curie_map.get())
        self.annot_id = assoc_id
        self.gene1 = gene1
        self.gene2 = gene2
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.relationships['orthologous']  # default
        self.pub_list = None

        self.setSubject(gene1)
        self.setObject(gene2)

        return

    def addOrthologyAssociationToGraph(self,g):

        self.addAssociationToGraph(g)

        return

    def addGeneFamilyToGraph(self,g,family_id):
        '''
        Make an association between a group of genes and a grouping class.
        We make the assumption that the genes in the association are part of the supplied
        family_id
        :param family_id:
        :param g: the graph to modify
        :return:
        '''

        f = URIRef(self.cu.get_uri(family_id))
        a = URIRef(self.cu.get_uri(self.gene1))
        b = URIRef(self.cu.get_uri(self.gene2))
        p = URIRef(self.cu.get_uri(self.relationships['has_member']))

        #make the assumption that the genes have already been added as classes previously
        #add the family grouping as an instance of a gene family?

        gene_family = 'DATA:3148'  #http://edamontology.org/data_3148
        g.add((f, RDF['type'], URIRef(self.cu.get_uri(gene_family))))  #Instance?

        #is a gene family annotation

        #add each gene to the family
        g.add((f, p, a))
        g.add((f, p, b))


        return
