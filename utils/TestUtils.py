from rdflib import Graph
from sources import Source
from sources.HPOAnnotations import HPOAnnotations
from sources.ZFIN import ZFIN
from sources.OMIM import OMIM
from sources.BioGrid import BioGrid
from sources.MGI import MGI
from sources.IMPC import IMPC
from sources.Panther import Panther
from sources.NCBIGene import NCBIGene
from sources.UCSCBands import UCSCBands
from sources.CTD import CTD
from sources.GeneReviews import GeneReviews
import os, hashlib


class TestUtils:

    def __init__(self, source=None):
        # instantiate source object
        self.source = source
        self.graph = Graph()
        if (source is not None):
            self.source.load_bindings()

        return

    def query_graph(self, query):
        query_result = self.graph.query(query)

        for row in query_result:
            print(row)

        return

    def check_query_syntax(self, query):
        self.source.graph.query(query)
        return

    def load_graph_from_turtle(self):
        file = self.source.outdir+'/'+self.source.name+'.ttl'
        if not os.path.exists(file):
            raise Exception("file:"+file+" does not exist")
        # load turtle file into graph
        self.graph.parse(file, format="turtle")

        return

    def get_file_md5(self, directory, file, blocksize=2**20):
        # reference: http://stackoverflow.com/questions/
        #            1131220/get-md5-hash-of-big-files-in-python

        md5 = hashlib.md5()
        with open(os.path.join(directory, file), "rb") as f:
            while True:
                buffer = f.read(blocksize)
                if not buffer:
                    break
                md5.update(buffer)

        return md5.hexdigest()
