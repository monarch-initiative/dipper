from rdflib.graph import ConjunctiveGraph
from rdflib.namespace import RDF, OWL
import argparse
import re


def main():
    parser = argparse.ArgumentParser(description='description',
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')
    args = parser.parse_args()
    property_list = get_properties_from_input(args.input)
    merged_graph = make_property_graph(property_list)

    # merge graphs
    merged_graph.parse(args.input, format="turtle")

    merged_graph.serialize(args.output, format="turtle")


def get_properties_from_input(file):
    input_graph = ConjunctiveGraph()
    input_graph.parse(file, format="turtle")

    query = """
                SELECT DISTINCT ?property
                WHERE {
                    ?subject ?property ?object .
                }
            """
    query_result = input_graph.query(query)
    # collapse to single list
    property_set = set()
    for row in query_result:
        property_set.add(row[0])

    return property_set


def make_property_graph(properties):
    graph = ConjunctiveGraph()
    output_graph = ConjunctiveGraph()

    ontologies = [
        'https://raw.githubusercontent.com/monarch-initiative/SEPIO-ontology/master/src/ontology/sepio.owl',
        'https://raw.githubusercontent.com/monarch-initiative/GENO-ontology/develop/src/ontology/geno.owl',
        'https://raw.githubusercontent.com/oborel/obo-relations/master/ro.owl',
        'http://purl.obolibrary.org/obo/iao.owl',
        'http://data.monarchinitiative.org/owl/ero.owl',
        'https://raw.githubusercontent.com/jamesmalone/OBAN/master/ontology/oban_core.ttl',
        'http://purl.obolibrary.org/obo/pco.owl',
        'http://purl.obolibrary.org/obo/xco.owl'
    ]

    for ontology in ontologies:
        print("parsing: " + ontology)
        if re.search(r'\.owl', ontology):
            graph.parse(ontology, format='xml')
        elif re.search(r'\.ttl', ontology):
            graph.parse(ontology, format='turtle')
        else:
            graph.parse(ontology)

    # Get object properties
    query = """
                SELECT ?property
                WHERE {
                    ?property a owl:ObjectProperty .
                }
            """
    query_result = graph.query(query)
    output_graph = add_property_to_graph(
        query_result, output_graph, OWL['ObjectProperty'], properties)

    # Get annotation properties
    query = """
               SELECT ?property
               WHERE {
                    ?property a owl:AnnotationProperty .
                }
            """
    query_result = graph.query(query)
    output_graph = add_property_to_graph(
        query_result, output_graph, OWL['AnnotationProperty'], properties)

    # Get data properties
    query = """
                SELECT ?property
                WHERE {
                    ?property a owl:DatatypeProperty .
                }
            """
    query_result = graph.query(query)
    output_graph = add_property_to_graph(
        query_result, output_graph, OWL['DatatypeProperty'], properties)

    return output_graph


def add_property_to_graph(results, graph, property_type, property_list):
    for row in results:
        if row[0] in property_list:
            graph.add((row[0], RDF['type'], property_type))
    return graph


if __name__ == "__main__":
    main()

