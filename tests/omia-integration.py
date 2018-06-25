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
    model_len = len(set(list(models)))

    if model_len < EXPECTED_PAIRS:
        logger.error(
            "Not enough model_of predicates in graph: found {}, "
            "expected {} check omia log for warnings".format(
                model_len, EXPECTED_PAIRS))
        exit(1)
    # else:
    #    logger.info(
    #        "Found {} model_of predicates in graph, expected at least: {}".format(
    #            model_len, EXPECTED_PAIRS))



    breed = 'https://monarchinitiative.org/model/OMIA-breed:758' 
    disease = 'http://purl.obolibrary.org/obo/OMIM_305100'

    omim_diseases = graph.objects(
        subject = URIRef(breed),
        predicate = model_of
    )

    if list(omim_diseases) != [URIRef(disease)]:
        logger.error("Missing breed to omim triple for ".format(breed))
        logger.error(list(omim_diseases))
        exit(1)

    logger.info("PASSED")


if __name__ == "__main__":
    main()
