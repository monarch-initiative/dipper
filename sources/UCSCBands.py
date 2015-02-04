import os
from stat import *
import urllib
from urllib import request
import re
import time
from datetime import datetime
import gzip,os.path
import json
from rdflib import Graph, Literal, URIRef, Namespace, BNode
from rdflib.namespace import RDF, RDFS, OWL, DC,XSD
import unicodedata

from sources.Source import Source
from models.D2PAssoc import D2PAssoc
from models.DispositionAssoc import DispositionAssoc
#from models.Features import Feature,makeChromID
from models.Dataset import Dataset
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
import config
import curie_map

class UCSCBands(Source):
    '''
    This will take the UCSC defintions of cytogenic bands and create the nested structures to enable
    overlap and containment queries.
    For example,
    13q21.31 ==>  13q21.31,  13q21.3,  13q21,  13q2,  13q, 13
    We leverage the Faldo model here for region definitions, and map each of the chromosomal parts to SO.
    At the moment, this only computes the bands for Human.
    We differentiate the species by adding the species id to the identifier prior to the chromosome number.
    The identifers therefore are generically created like:
    <species number>chr<num><band>
    We will create triples for a given band like:
    :9606chr1p36.33 rdf[type] SO:chromosome_band, faldo:Region
    :9606chr1p36 subsequence_of :9606chr1p36.3
    :9606chr1p36 faldo:location [ faldo:begin 0 faldo:end 2300000]

    where any band in the file is an instance of a chr_band (or a more specific type), is a subsequence
    of it's containing region, and is located in the specified coordinates.

    we determine the containing regions of the band by parsing the band-string; since each alphanumeric
    is a significant "place", we can split it with the shorter strings being parents of the longer string

    TODO: will can then sort the locations of the annotated bands, and propagate them to the
    intermediate/parent regions

    TODO: any species by commandline argument

    TODO: abstract this out into a model
    '''

    files = {
        'human_bands' : {
            'file' : '9606cytoBand.txt.gz',
            'url' : 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz',
            'taxon' : '9606'
        },
    }


    relationships = {
        'gene_product_of' : 'RO:0002204',
        'has_gene_product' : 'RO:0002205',
        'is_about' : 'IAO:00000136',
        'has_subsequence' : 'RO:0002524',
        'is_subsequence_of' : 'RO:0002525',
    }


    def __init__(self):
        Source.__init__(self, 'ucscbands')

        self.load_bindings()

        self.dataset = Dataset('ucscbands', 'UCSC Cytogenic Bands', 'http://hgdownload.cse.ucsc.edu')

        #data-source specific warnings (will be removed when issues are cleared)
        #print()

        return

    def fetch(self, is_dl_forced):

        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                ('/').join((self.rawdir,file['file'])),
                                is_dl_forced)
            self.dataset.setFileAccessUrl(file['url'])
            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        self.dataset.setVersion(filedate)

        return


    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows")

        print("Parsing files...")

        self._get_chrbands(limit,self.files['human_bands'])
        #TODO as a post-processing step, we must propagate the coordinates to the upper-level features

        self.load_core_bindings()
        self.load_bindings()

        print("Done parsing files.")

        return

    def _get_chrbands(self,limit,data):
        '''
        :param limit:
        :return:
        '''
        cu = CurieUtil(curie_map.get())

        print("INFO: Processing Chr bands from",)
        line_counter=0
        myfile=('/').join((self.rawdir,data['file']))
        print("FILE:",myfile)
        taxon = data['taxon']
        mybands = {}

        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                if (re.match('^#',line)):
                    continue

                #chr13	4500000	10000000	p12	stalk
                (chrom,start,stop,band,rtype) = line.split('\t')
                line_counter += 1

                cid = makeChromID(chrom,taxon)
                tax_id = (':').join(('NCBITaxon',taxon))
                region_type_id = self._map_type_of_region(rtype)

                #add the chr
                cfeature = Feature(cid,chrom,self._map_type_of_region('chromosome'))
                cfeature.addFeatureToGraph(self.graph) #fixme - should probably only add this at the end?
                cfeature.addTaxonToFeature(self.graph,tax_id)
                #add the chr to the hashmap of coordinates
                if chrom not in mybands.keys():
                    mybands[chrom] = {'min' : 0, 'max' : 0 }

                #TODO add genome build as reference to the chromosomes (or at minimum to the dataset)
                #add the region and it's location
                maploc_id = cid+band
                bfeature = Feature(maploc_id,chrom+band,region_type_id,rtype)
                bfeature.addFeatureToGraph(self.graph)
                bfeature.addCoordinatesOfFeature(self.graph,start,stop)

                #get the parent bands, and make them unique
                parents = list(self._make_parent_bands(band,set()))
                #alphabetical sort will put them in smallest to biggest
                parents.sort(reverse=True)
                #print('parents of',chrom,band,':',parents)

                #add the parents to the graph, in hierarchical order
                #TODO this is somewhat inefficient due to re-adding upper-level nodes when iterating over the file
                for i in range(len(parents)):
                    pid = cid+parents[i]
                    if (re.match('[pq]$',parents[i])):
                        rti = self._map_type_of_region('chromosome_arm')
                    else:
                        #FIXME there really ought to be a feature that is broader than a band, but part of an arm
                        rti = self._map_type_of_region('chromosome_part')

                    pfeature = Feature(pid,chrom+parents[i],rti)
                    pfeature.addFeatureToGraph(self.graph)
                    #add the relationships to the parent
                    if (i < len(parents)-1):
                        ppid = cid+parents[i+1]
                        pfeature.addSubsequenceOfFeature(self.graph,ppid)
                    else:
                        #add the last one (p or q usually) as attached to the chromosome
                        pfeature.addSubsequenceOfFeature(self.graph,cid)

                #connect the band here to the first one in the parent list
                bfeature.addSubsequenceOfFeature(self.graph,cid+parents[0])

                #Here, we add the parents to a hashmap of chr bands to propagate the chromosomal coords
                for p in parents:
                    k = chrom+p
                    sta=int(start)
                    sto=int(stop)
                    if k not in mybands.keys():
                        b = {'min' : min(sta,sto), 'max' : max(sta,sto)}
                        mybands[k] = b
                    else:
                        b = mybands.get(k)
                        b['min'] = min(sta,sto,b['min'])
                        b['max'] = max(sta,sto,b['max'])
                        mybands[k] = b
                        #also, set the max for the chrom
                        c = mybands.get(chrom)
                        c['max'] = max(sta,sto,c['max'])
                        mybands[chrom] = c

                if (limit is not None and line_counter > limit):
                    break

        #add the band locations to the graph.
        for b in mybands.keys():
            bid = makeChromID(b,taxon)
            bfeature = Feature(bid,None,None)  #for now, hopefully it won't overwrite
            bfeature.addCoordinatesOfFeature(self.graph,mybands.get(b)['min'],mybands.get(b)['max'])
        return

    def _make_parent_bands(self,band,child_bands):
        '''
        #this will determine the grouping bands that it belongs to, recursively
        #13q21.31 ==>  13, 13q, 13q2, 13q21, 13q21.3, 13q21.31

        :param chrom:
        :param band:
        :param child_bands:
        :return:
        '''
        m=re.match('([pq]\d+(?:\.\d+)?)',band)
        if (len(band) > 0):
            if (m):
                p=str(band[0:len(band)-1])
                p = re.sub('\.$','',p)
                if p is not None:
                    child_bands.add(p)
                    self._make_parent_bands(p,child_bands)
        else:
            child_bands = set()
        return child_bands


    def _map_type_of_region(self,type):
        '''
        Note that interband seems specific for polytene chromosomes, but not for
        negatively stained regions of regular chromosomes.  Assigning to chr_band for now.
        :param type:
        :return:
        '''
        so_id = 'SO:0000830'
        type_to_so_map = {
            'acen' : 'SO:0001795', #TODO check using regional centromere
            'gvar' : 'SO:0000628', #chromosomal structural element
            'stalk' : 'SO:0000830', #FIXME using chromosome part for now
            'gneg' : 'SO:0000341', #FIXME.  using chr_band
            'gpos100' : 'SO:0000341',  #FIXME.  using chr_band
            'gpos25' : 'SO:0000341',   #FIXME.  using chr_band
            'gpos50' : 'SO:0000341',  #FIXME.  using chr_band
            'gpos75' : 'SO:0000341',  #FIXME.  using chr_band
            'chromosome' : 'SO:0000340',
            'chromosome_arm' : 'SO:0000105',
            'chromosome_band' : 'SO:0000341',
            'chromosome_part' : 'SO:0000830'
        }

        if (type in type_to_so_map):
            so_id = type_to_so_map.get(type)
        else:
            print("WARN: unmapped code",type,". Defaulting to chr_part 'SO:0000830'.")

        return so_id

