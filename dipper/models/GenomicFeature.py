__author__ = 'nicole'

import re

from rdflib import URIRef,Literal
from rdflib.namespace import RDF,XSD

from dipper.curie import curie_map
from dipper.utils.GraphUtils import GraphUtils
from dipper.utils.CurieUtil import CurieUtil
from dipper.models.Assoc import Assoc


class Feature() :
    '''
    Dealing with genomic features here.  By default they are all faldo:Regions.
    We use SO for typing genomic features.
    At the moment, RO:has_subsequence is the default relationship between the regions, but this should be tested/verified.

    TODO: the graph additions are in the addXToFeature functions, but should be separated.
    TODO: this will need to be extended to properly deal with fuzzy positions in faldo.
    '''

    relationships = {
        'location' : 'faldo:location',
        'gene_product_of' : 'RO:0002204',
        'has_gene_product' : 'RO:0002205',
        'is_about' : 'IAO:00000136',
        'has_subsequence' : 'RO:0002524',
        'is_subsequence_of' : 'RO:0002525',
    }

    types = {
        'region' : 'faldo:Region',
        'begin' : 'faldo:begin',
        'end' : 'faldo:end',
        'chromosome' : 'SO:0000340',
        'chromosome_arm' : 'SO:0000105',
        'chromosome_band' : 'SO:0000341',
        'chromosome_part' : 'SO:0000830'

    }

    def __init__(self,id,label,type,description=None):
        self.id = id
        self.label = label
        self.type = type
        self.description = description
        self.gu = GraphUtils(curie_map.get())
        self.cu = CurieUtil(curie_map.get())
        return

    def addFeatureToGraph(self,graph):
        '''
        We make the assumption here that all features are instances, and also that they
        are regions
        :param graph:
        :return:
        '''
        self.gu.addIndividualToGraph(graph,self.id,self.label,self.type,self.description)
        if re.match('SO',self.type):
            n = URIRef(self.cu.get_uri(self.id))
            region = URIRef(self.cu.get_uri(self.types['region']))
            graph.add((n,RDF['type'],region))
        return

    def addCoordinatesOfFeature(self,graph,start,stop):
        '''
        Given the start and stop coordinates, it will add two triples for the feature, using
        faldo:begin and faldo:end as int literals
        :param graph:
        :param start:
        :param stop:
        :return:
        '''
        self.start = start
        self.stop = stop
        n = URIRef(self.cu.get_uri(self.id))
        begin = URIRef(self.cu.get_uri(self.types['begin']))
        end = URIRef(self.cu.get_uri(self.types['end']))
        graph.add((n,begin,Literal(start,datatype=XSD['integer'])))
        graph.add((n,end,Literal(stop,datatype=XSD['integer'])))
        return

    def addSubsequenceOfFeature(self,graph,parentid):
        '''
        This will add a triple like:
        feature subsequence_of parent
        :param graph:
        :param parentid:
        :return:
        '''
        n = URIRef(self.cu.get_uri(self.id))
        p = URIRef(self.cu.get_uri(parentid))
        subsequence=URIRef(self.cu.get_uri(self.relationships['is_subsequence_of']))
        graph.add((n,subsequence,p))
        return


    def addTaxonToFeature(self,graph,taxonid):
        '''
        Given the taxon id, this will add the following triple:
        feature in_taxon taxonid
        :param graph:
        :param id:
        :param taxonid:
        :return:
        '''
        self.taxon = taxonid
        n = URIRef(self.cu.get_uri(self.id))
        t = URIRef(self.cu.get_uri(self.taxon))
        intaxon=URIRef(self.cu.get_uri(Assoc.relationships['in_taxon']))
        graph.add((n,intaxon,t))
        return

def makeChromID(chrom, taxon=None):
    '''
    This will take a chromosome number and a NCBI taxon number,
    and create a unique identifier for the chromosome.  These identifiers
    are made in the @base space like:
    Homo sapiens (9606) chr1 ==> :9606chr1
    Mus musculus (10090) chrX ==> :10090chrX

    :param chrom: the chromosome (preferably without any chr prefix)
    :param taxon: the numeric portion of the taxon id
    :return:
    '''
    if (taxon is None):
        print('WARN: no taxon for this chrom.  you may have conflicting ids')
        taxon = ''
    # replace any chr-like prefixes with blank to standardize
    c = re.sub('ch(r?)[omse]*', '', chrom)
    id = ('').join((':', taxon, 'chr', c))
    return id

