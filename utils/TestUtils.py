from rdflib import Graph
import os


class TestUtils:

    def __init__(self):
        # instantiate source object
        self.graph = Graph()

        return

    def query_graph(self, query):
        query_result = self.graph.query(query)

        for row in query_result:
            print(", ".join(row))

        return

    def check_query_syntax(self, query, source):
        source.load_bindings()
        source.graph.query(query)
        return

    def load_graph_from_turtle(self, source):
        file = source.outdir+'/'+source.name+'.ttl'
        if not os.path.exists(file):
            raise Exception("file:"+file+" does not exist")
        # load turtle file into graph
        self.graph.parse(file, format="turtle")

        return