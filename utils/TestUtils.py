from rdflib import Graph
import os


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