import logging
import io
from dipper.graph.RDFGraph import RDFGraph
from unittest.mock import patch, mock_open

logger = logging.getLogger(__name__)


class TestUtils:

    def test_graph_equality(self, turtlish, graph):
        """

        :param turtlish: String of triples in turtle
                         format without prefix header
        :param graph: Graph object to test against
        :return: Boolean, True if graphs contain same
                          set of triples
        """
        turtle_graph = RDFGraph()
        turtle_graph.bind_all_namespaces()
        prefixes = "\n".join(
            ["@prefix {}: <{}> .".format(n[0], n[1])
            for n in turtle_graph.namespace_manager.namespaces()]
        )

        turtle_string = prefixes + turtlish
        mock_file = io.StringIO(turtle_string)
        turtle_graph.parse(mock_file, format="turtle")
        turtle_triples = set(list(turtle_graph))
        ref_triples = set(list(graph))
        equality = turtle_triples == ref_triples
        if not equality:
            logger.warning("Triples do not match\n"
                           "Left hand difference: {}\n"
                           "Right hand difference:{}".format(
                turtle_triples - ref_triples,
                ref_triples - turtle_triples
            ))
        return equality
