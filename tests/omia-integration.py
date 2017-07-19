from rdflib.graph import ConjunctiveGraph, URIRef
from rdflib import util as rdflib_util
import argparse
import logging

"""
OMIA requires integration across three sources:
OMIA, OMIM, and NCBI

There are some odd bugs, see https://github.com/monarch-initiative/dipper/issues/417

This script ensures we are getting the full set of model_of relationships
"""
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
EXPECTED_PAIRS = 175


def main():
    parser = argparse.ArgumentParser(
        description='OMIA integration test',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        '--input', '-i', type=str, required=True,
        help='Location of input ttl file')

    args = parser.parse_args()

    graph = ConjunctiveGraph()
    graph.parse(args.input, format=rdflib_util.guess_format(args.input))

    model_of = URIRef('http://purl.obolibrary.org/obo/RO_0003301')

    models = graph.subject_objects(model_of)
    model_len = len(list(models))

    if model_len < EXPECTED_PAIRS:
        logger.error("Not enough model_of predicates in graph:"
                     " {} expected {} check omia log for"
                     " warnings".format(model_len, EXPECTED_PAIRS))
        exit(1)
    else:
        logger.info("PASSED")

if __name__ == "__main__":
    main()

