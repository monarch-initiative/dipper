from rdflib.graph import ConjunctiveGraph, URIRef
from rdflib.namespace import RDF, OWL, DCTERMS
from rdflib import util as rdflib_util
from xml.sax import SAXParseException
import argparse
import re
import logging

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='description',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        '--input', '-i', type=str, required=True,
        help='Location of input file')

    parser.add_argument(
        '--output', '-o', type=str, required=True,
        help='Location of output file')

    parser.add_argument(
        '--input_format', '-f', type=str, default="turtle",
        help='format of source rdf file (turtle, nt, rdf/xml)')

    parser.add_argument(
        '--output_format', '-g', type=str, default="turtle",
        help='format of target rdf file (turtle, nt, rdf/xml)')

    args = parser.parse_args()
    property_list = get_properties_from_input(args.input, args.input_format)
    merged_graph = make_property_graph(property_list, args)

    # merge graphs
    merged_graph.parse(args.input, format=args.input_format)

    merged_graph.serialize(args.output, format=args.output_format)


def get_properties_from_input(file, input_format):
    input_graph = ConjunctiveGraph()
    input_graph.parse(file, format=input_format)

    # collapse to single list
    property_set = list()
    for row in input_graph.predicates():
        property_set.append(row)

    return set(property_set)


def make_property_graph(properties, args):
    graph = ConjunctiveGraph()
    output_graph = ConjunctiveGraph()

    GH = 'https://raw.githubusercontent.com'
    OBO = 'https://purl.obolibrary.org/obo'
    ontologies = [
        OBO + '/sepio.owl',
        OBO + '/geno.owl',
        OBO + '/iao.owl',
        OBO + '/pco.owl',
        OBO + '/xco.owl',
        OBO + '/ro.owl',
        GH + '/jamesmalone/OBAN/master/ontology/oban_core.ttl',
    ]

    for ontology in ontologies:
        print("parsing: " + ontology)
        try:
            graph.parse(ontology, format=rdflib_util.guess_format(ontology))
        except SAXParseException as e:
            logger.error(e)
            logger.error('Retrying: ' + ontology)
            graph.parse(ontology, format="turtle")
        except OSError as e:  # URLError:
            # simple retry
            logger.error(e)
            logger.error('Retrying: ' + ontology)
            graph.parse(ontology, format=rdflib_util.guess_format(ontology))

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

    # Hardcoded properties
    output_graph.add(
        (URIRef('https://monarchinitiative.org/MONARCH_cliqueLeader'),
            RDF['type'], OWL['AnnotationProperty']))

    output_graph.add(
        (URIRef('https://monarchinitiative.org/MONARCH_anonymous'),
            RDF['type'], OWL['AnnotationProperty']))

    # Check monarch data triple
    data_url = "https://data.monarchinitiative.org/ttl/{0}".format(
        re.sub(r".*/", "", args.input))
    new_url = "https://data.monarchinitiative.org/ttl/{0}".format(
        re.sub(r".*/", "", args.output))
    if (URIRef(data_url), RDF.type, OWL['Ontology']) in output_graph:
        output_graph.remove(URIRef(data_url), RDF.type, OWL['Ontology'])

    output_graph.add((URIRef(new_url), RDF.type, OWL['Ontology']))

    for row in output_graph.predicates(
            DC['source'], OWL['AnnotationProperty']):
        if row == RDF['type']:
            output_graph.remove(
                (DC['source'], RDF['type'], OWL['AnnotationProperty']))

    output_graph.add((DC['source'], RDF['type'], OWL['ObjectProperty']))

    return output_graph


def add_property_to_graph(results, graph, property_type, property_list):
    for row in results:
        if row in property_list:
            graph.add((row, RDF['type'], property_type))
    return graph


if __name__ == "__main__":
    main()
