import rdflib
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
import os


class TestUtils:

    def __init__(self, source):
        # instantiate source object
        self.source = source
        self.graph = Graph()
        self._load_graph_from_turtle()
        return

    def query_graph(self, query):
        query_result = self.graph.query(query)

        for row in query_result:
            print("%s %s" % row)

        return

    def _load_graph_from_turtle(self):
        file = self.source.outdir+'/'+self.source.name+'.ttl'
        if not os.path.exists(file):
            raise Exception("file:"+file+" does not exist")
        # load turtle file into graph
        self.graph.parse(file, format="turtle")

        return