from rdflib.graph import ConjunctiveGraph
from rdflib.namespace import RDF, OWL, DC
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
    parser.add_argument('--format', '-f', type=str, default="turtle",
                        help='format of rdf file (turtle, n3, rdf/xml)')
    args = parser.parse_args()
    property_list = get_properties_from_input(args.input, args.format)
    merged_graph = make_property_graph(property_list)

    # merge graphs
    merged_graph.parse(args.input, format="turtle")

    merged_graph.serialize(args.output, format="turtle")


def get_properties_from_input(file, format):
    input_graph = ConjunctiveGraph()
    input_graph.parse(file, format=format)

    # collapse to single list
    property_set = set()
    for row in input_graph.predicates():
        property_set.add(row)

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
    output_graph = add_property_to_graph(
        graph.subjects(RDF['type'], OWL['ObjectProperty']),
        output_graph, OWL['ObjectProperty'], properties)

    # Get annotation properties
    output_graph = add_property_to_graph(
        graph.subjects(RDF['type'], OWL['AnnotationProperty']),
        output_graph, OWL['AnnotationProperty'], properties)

    # Get data properties
    output_graph = add_property_to_graph(
        graph.subjects(RDF['type'], OWL['DatatypeProperty']),
        output_graph, OWL['DatatypeProperty'], properties)

    output_graph.add((DC['source'], RDF['type'], OWL['ObjectProperty']))

    return output_graph


def add_property_to_graph(results, graph, property_type, property_list):
    for row in results:
        if row in property_list:
            graph.add((row, RDF['type'], property_type))
    return graph


if __name__ == "__main__":
    main()

