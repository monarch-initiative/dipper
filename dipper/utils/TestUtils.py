import logging
import io
from dipper.graph.RDFGraph import RDFGraph
from rdflib import URIRef, RDF

LOG = logging.getLogger(__name__)


class TestUtils:

    @staticmethod
    def test_graph_equality(turtlish, graph):
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
            ["@prefix {}: <{}> .".format(
                n[0], n[1]) for n in turtle_graph.namespace_manager.namespaces()]
        )

        turtle_string = prefixes + turtlish
        mock_file = io.StringIO(turtle_string)
        turtle_graph.parse(mock_file, format="turtle")

        TestUtils.remove_ontology_axioms(graph)

        turtle_triples = set(list(turtle_graph))
        ref_triples = set(list(graph))
        equality = turtle_triples == ref_triples
        if not equality:
            LOG.warning(
                "Triples do not match\n"
                "\tLeft hand difference: %s\n"
                "\tRight hand difference: %s",
                sorted(turtle_triples - ref_triples),
                sorted(ref_triples - turtle_triples)
            )
        return equality

    @staticmethod
    def remove_ontology_axioms(graph):
        """
        Given an rdflib graph, remove any triples
        connected to an ontology node:
        {} a owl:Ontology
        :param graph: RDFGraph
        :return: None
        """
        ontology_iri = URIRef("http://www.w3.org/2002/07/owl#Ontology")

        for subject in graph.subjects(RDF.type, ontology_iri):
            for predicate, obj in graph.predicate_objects(subject):
                graph.remove((subject, predicate, obj))
            graph.remove((subject, RDF.type, ontology_iri))

