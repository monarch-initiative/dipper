from rdflib.graph import ConjunctiveGraph
from rdflib.namespace import RDF, OWL


def main():
    make_property_graph()


def make_property_graph():
    graph = ConjunctiveGraph()
    output_graph = ConjunctiveGraph()

    ontologies = [
        'https://raw.githubusercontent.com/monarch-initiative/SEPIO-ontology/master/src/ontology/sepio.owl',
        'https://raw.githubusercontent.com/monarch-initiative/GENO-ontology/develop/src/ontology/geno.owl',
        'https://raw.githubusercontent.com/oborel/obo-relations/master/ro.owl'
    ]

    for ontology in ontologies:
        print("parsing: " + ontology)
        graph.parse(ontology, format='xml')

    # Get object properties
    query = """
                SELECT ?property
                WHERE {
                    ?property a owl:ObjectProperty .
                }
            """
    query_result = graph.query(query)
    output_graph = add_property_to_graph(
        query_result, output_graph, OWL['ObjectProperty'])

    # Get annotation properties
    query = """
               SELECT ?property
               WHERE {
                    ?property a owl:AnnotationProperty .
                }
            """
    query_result = graph.query(query)
    output_graph = add_property_to_graph(
        query_result, output_graph, OWL['AnnotationProperty'])

    # Get data properties
    query = """
                SELECT ?property
                WHERE {
                    ?property a owl:DatatypeProperty .
                }
            """
    query_result = graph.query(query)
    output_graph = add_property_to_graph(
        query_result, output_graph, OWL['DatatypeProperty'])
    output_graph.serialize('property-graph.ttl', format="turtle")

    return


def add_property_to_graph(results, graph, property_type):
    for row in results:
        graph.add((row[0], RDF['type'], property_type))
    return graph


if __name__ == "__main__":
    main()

