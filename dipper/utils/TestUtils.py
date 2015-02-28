from rdflib import Graph
import os
import logging
import sys

logger = logging.getLogger(__name__)


class TestUtils:

    def __init__(self, graph=None):
        # instantiate source object
        self.graph = graph
        if self.graph is None:
            self.graph = Graph()

        return

    def query_graph(self, query, is_formatted=True):
        query_result = self.graph.query(query)
        output = []
        for row in query_result:
            if is_formatted:
                result_set = []
                for val in row:
                    if val is None:
                        val = 'null'
                    result_set.append(val)
                    output.append(", ".join(result_set))
            else:
                output.append(row)

        return output

    def check_query_syntax(self, query, source):
        source.load_bindings()
        source.graph.query(query)
        return

    def load_graph_from_turtle(self, source):
        file = source.outdir+'/'+source.name+'.ttl'
        if not os.path.exists(file):
            logger.error("file: %s does not exist", file)
            sys.exit(1)
        # load turtle file into graph
        self.graph.parse(file, format="turtle")

        return