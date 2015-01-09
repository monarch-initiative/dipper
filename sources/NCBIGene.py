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
from sources.Source import Source

from models.D2PAssoc import D2PAssoc
from models.DispositionAssoc import DispositionAssoc
from models.Dataset import Dataset
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
import config

class NCBIGene(Source):
    '''
     OMIM is an unusual source.  We can get lots of the disease-gene associations, including allelic variants
     from their ftp site, which is obtainable anonymously.  However, more detailed information is available
     via their API.  So, we pull the basic files from their ftp site, extract the omim identifiers,
     then query their API in batch.  (Note this requires an apiKey, which is not stored in the repo,
     but in a separate conf.json file.)
    '''

    files = {
        'gene_info' : {
            'file' : 'gene_info.gz',
            'url' : 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'
        },
        #'gene_history' : {
        #    'file': 'gene_history.gz',
        #    'url' : 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz'
        #},
        #'gene2pubmed' : {
        #    'file': 'gene2pubmed.gz',
        #    'url' : 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz'
        #},
    }

    prefixes = {
        'NCBIGene' : 'http://ncbi.nlm.nih.gov/gene/',
        'faldo' : 'http://biohackathon.org/resource/faldo#',
        'NCBITaxon' : 'http://ncbi.nlm.nih.gov/taxonomy/',
        'SO' : 'http://purl.obolibrary.org/obo/SO_'
    }

    curie_map = {}


    def __init__(self):
        Source.__init__(self, 'ncbigene')

        self.curie_map.update(D2PAssoc.curie_map)
        self.curie_map.update(DispositionAssoc.curie_map)
        self.curie_map.update(self.prefixes)

        self.load_bindings()

        self.dataset = Dataset('ncbigene', 'National Center for Biotechnology Information', 'http://ncbi.nih.nlm.gov/gene')
        #data-source specific warnings (will be removed when issues are cleared)
        #print()

        return

    def fetch(self):
        #this is fetching the standard files, not from the API/REST service
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],('/').join((self.rawdir,file['file'])))
            self.dataset.setFileAccessUrl(file['url'])
            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        self.dataset.setVersion(filedate)

        return


    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes: (Nothing)
        :return: None
        '''
        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.curie_map.keys():
            v=self.curie_map[k]
            self.graph.bind(k, Namespace(v))
        return

    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows")

        print("Parsing files...")


        g = self.graph

        gu = GraphUtils(self.curie_map)
        cu = CurieUtil(self.curie_map)

        self._get_gene_info(limit)

        ##### Write it out #####
        filewriter = open(self.outfile, 'w')
        self.load_core_bindings()
        self.load_bindings()

        print("Finished parsing files. Writing turtle to",self.outfile)
        print(g.serialize(format="turtle").decode(),file=filewriter)
        filewriter.close()


        return

    def _get_gene_info(self,limit):
        exact=URIRef(Assoc.relationships['hasExactSynonym'])
        related=URIRef(Assoc.relationships['hasRelatedSynonym'])
        intaxon=URIRef(Assoc.relationships['in_taxon'])

        #an omim-specific thing here; from the omim.txt.gz file, get the omim numbers
        #not unzipping the file
        print("INFO: Processing Gene records")
        line_counter=0
        myfile=('/').join((self.rawdir,self.files['gene_info']['file']))
        print("FILE:",myfile)
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                line_counter += 1
                if (re.match('^#',line)):
                    continue
                (tax_num,gene_num,symbol,locustag,
                 synonyms,xrefs,chr,map_loc,desc,
                 gtype,authority_symbol,name,
                 nomenclature_status,other_designations,modification_date) = line.split('\t')


                gene_id = (':').join(('NCBIGene',gene_num))
                tax_id = (':').join(('NCBITaxon',tax_num))
                gene_type_id = self._map_type_of_gene(gtype)

                gu = GraphUtils(self.curie_map)
                cu = CurieUtil(self.curie_map)
                n = URIRef(cu.get_uri(gene_id))
                t = URIRef(cu.get_uri(gene_type_id))
                gu.addClassToGraph(self.graph,gene_id,symbol,None,desc)

                #we are making genes classes, not instances.  so add it as a subclass here
                self.graph.add((n,Assoc.OWLSUBCLASS,t))
                self.graph.add((n,exact,Literal(name)))
                for s in synonyms:
                    self.graph.add((n,related,Literal(s)))
                self.graph.add((n,intaxon,URIRef(cu.get_uri(tax_id))))


                #deal with coordinate information:
                begin = URIRef(cu.get_uri('faldo:begin'))
                end = URIRef(cu.get_uri('faldo:end'))
                region = URIRef(cu.get_uri('faldo:Region'))
                loc = URIRef(cu.get_uri('faldo:location'))
                #make blank nodes for the regions, and positions

                generegion = BNode((':').join((chr,map_loc)))
                self.graph.add((n,loc,generegion))
                self.graph.add((generegion,RDF['type'],region))
                self.graph.add((generegion,begin,Literal(map_loc)))
                self.graph.add((generegion,end,Literal(map_loc)))

                #todo add reference and strand info

                if (limit is not None and line_counter > limit):
                    break

        return

    def _map_type_of_gene(self,type):
        so_id = 'SO:0000704'
        type_to_so_map = {
            'ncRNA': 'SO:0001263',
            'other': 'SO:0000704',
            'protein-coding': 'SO:0001217',
            'pseudo': 'SO:0000336',
            'rRNA': 'SO:0001637',
            'snRNA': 'SO:0001268',
            'snoRNA': 'SO:0001267',
            'tRNA': 'SO:0001272',
            'unknown': 'SO:0000704'
        }

        if (type in type_to_so_map):
            so_id = type_to_so_map.get(type)
        else:
            print("WARN: unmapped code",type,". Defaulting to 'SO:0000704'.")

        return so_id


        return so_id