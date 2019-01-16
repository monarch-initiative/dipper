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
LOG = logging.getLogger(__name__)
EXPECTED_PAIRS = 175


def main():
    parser = argparse.ArgumentParser(
        description='OMIA integration test',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        '--input', '-i', type=str, required=True, help='Location of input ttl file')

    args = parser.parse_args()

    graph = ConjunctiveGraph()
    graph.parse(args.input, format=rdflib_util.guess_format(args.input))

    # "is model of": "RO:0003301"
    # is_model_of = URIRef('OBO:RO_0003301')
    is_model_of = URIRef('http://purl.obolibrary.org/obo/RO_0003301')

    # if we curie_map & globaltt here we could ...
    # (pfx lcl) = globaltt["is model of"].split(':')
    # iri = curie_map[pfx] + '_'.join((pfx, lcl))
    # is_model_of = URIRef(iri)

    models = graph.subject_objects(is_model_of)
    model_len = len(set(list(models)))

    if model_len < EXPECTED_PAIRS:
        LOG.error(
            "Not enough <RO:is model of> predicates in graph: found {}, "
            "expected {} check omia log for warnings".format(
                model_len, EXPECTED_PAIRS))
        exit(1)
    # else:
    #    LOG.info(
    #        "Found {} model_of predicates in graph, expected at least: {}".format(
    #            model_len, EXPECTED_PAIRS))

    breed = 'https://monarchinitiative.org/model/OMIA-breed:758'
    disease = 'http://omim.org/entry/305100'

    omim_diseases = graph.objects(
        subject=URIRef(breed),
        predicate=is_model_of
    )

    if list(omim_diseases) != [URIRef(disease)]:
        LOG.error("Missing breed to omim triple for %s", breed)
        LOG.error(list(omim_diseases))
        exit(1)

    LOG.info("PASSED")


if __name__ == "__main__":
    main()
