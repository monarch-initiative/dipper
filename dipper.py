#!/usr/bin/env python3

__author__ = 'nlw'

import argparse
import logging
import unittest
import importlib

from tests.test_general import GeneralGraphTestCase
from dipper.utils.TestUtils import TestUtils

test_suite = unittest.TestLoader().loadTestsFromTestCase(GeneralGraphTestCase)


def main():
    source_to_class_map = {
        'hpoa': 'HPOAnnotations',  # ~3 min
        'zfin': 'ZFIN',
        'omim': 'OMIM',  # full file takes ~15 min, due to required throttling
        'biogrid': 'BioGrid',  # interactions file takes <10 minutes
        'mgi': 'MGI',
        'impc': 'IMPC',
        'panther': 'Panther',  # this takes a very long time, ~1hr to map 7 species-worth of associations
        'ncbigene': 'NCBIGene',  # takes about 4 minutes to process 2 species
        'ucscbands': 'UCSCBands',
        'ctd': 'CTD',
        'genereviews': 'GeneReviews',
        'eom': 'EOM',  # Takes about 5 seconds.
        'coriell': 'Coriell',
        'clinvar': 'ClinVar',
        'monochrom': 'Monochrom',
        'kegg': 'KEGG',
        'animalqtldb': 'AnimalQTLdb',
        'ensembl': 'Ensembl',
        'hgnc': 'HGNC',
        'orphanet': 'Orphanet',
        'omia': 'OMIA',
        'flybase': 'FlyBase',
        'mmrrc' : 'MMRRC',
        'wormbase': 'WormBase',
        'mpd': 'MPD',
        'gwascatalog': 'GWASCatalog'
        # 'monarch': Monarch,
        # 'go': GeneOntology
    }

    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Dipper: Data Ingestion'
                                                 ' Pipeline for SciGraph',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--sources', type=str, required=True,
                        help='comma separated list of sources')
    parser.add_argument('-l', '--limit', type=int, help='limit number of rows')
    parser.add_argument('--parse_only', action='store_true',
                        help='parse files without writing'),
    parser.add_argument('--fetch_only', action='store_true',
                        help='fetch sources without parsing')
    parser.add_argument('-f', '--force', action='store_true',
                        help='force re-download of files')
    parser.add_argument('--no_verify', help='ignore the verification step',
                        action='store_true')
    parser.add_argument('--query', help='enter in a sparql query', type=str)
    parser.add_argument('-q', '--quiet', help='turn off info logging',
                        action="store_true")
    parser.add_argument('--debug', help='turn on debug logging',
                        action="store_true")
    parser.add_argument('--skip_tests', help='skip any testing', action="store_true")

    # BNodes can't be visualized in Protege, so you can materialize them for testing purposes with this flag
    parser.add_argument('-nb', '--no_bnodes', help="convert blank nodes into identified nodes", action="store_true")

    # TODO this preconfiguration should probably live in the conf.json, and the same filter be applied to all sources
    parser.add_argument('-t', '--taxon', type=str,
                        help='Add a taxon constraint on a source. Enter 1+ NCBITaxon numbers, comma delimited\n'
                             'Implemented taxa per source\n'
                             'NCBIGene: 9606,10090,7955\n'
                             'Panther: 9606,10090,10116,7227,7955,6239,8355\n'
                             'BioGrid: 9606,10090,10116,7227,7955,6239,8355\n'
                             'UCSCBands: 9606')
    parser.add_argument('-o', '--test_only', help='only process and output the pre-configured test subset',
                        action="store_true")

    parser.add_argument('--format', help='serialization format: turtle (default), xml, n3, nt, raw', type=str)


    args = parser.parse_args()
    tax_ids = None
    if args.taxon is not None:
        tax_ids = list(map(int, args.taxon.split(',')))

    taxa_supported = ['Panther', 'NCBIGene', 'BioGrid', 'UCSCBands']

    formats_supported = ['xml', 'n3', 'turtle', 'nt', 'ttl', 'raw']

    if args.quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        if args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

    if args.no_bnodes is True:
        logger.info("Will materialize all BNodes into BASE space")

    if args.query is not None:
        test_query = TestUtils()
        for source in args.sources.split(','):
            source = source.lower()
            mysource = source_to_class_map[source]()

            # import source lib
            module = "dipper.sources.{0}".format(mysource)
            imported_module = importlib.import_module(module)
            source_class = getattr(imported_module, mysource)

            test_query.check_query_syntax(args.query, source_class)
            test_query.load_graph_from_turtle(source_class)

        print(test_query.query_graph(args.query, True))
        exit(0)

    # run initial tests
    if (args.no_verify or args.skip_tests) is not True:
        unittest.TextTestRunner(verbosity=2).run(test_suite)

    # set serializer
    if args.format is not None:
        if args.format in formats_supported:
            if args.format == 'ttl':
                args.format = 'turtle'
        else:
            logger.error("You have specified an invalid serializer: %s", args.format)

            exit(0)
    else:
        args.format = 'turtle'

    # iterate through all the sources
    for source in args.sources.split(','):
        logger.info("\n******* %s *******", source)
        source = source.lower()
        src = source_to_class_map[source]

        # import source lib
        module = "dipper.sources.{0}".format(src)
        imported_module = importlib.import_module(module)
        source_class = getattr(imported_module, src)
        mysource = None
        if src in taxa_supported:
            mysource = source_class(tax_ids)
        else:
            mysource = source_class()
        if args.parse_only is False:
            mysource.fetch(args.force)

        mysource.settestonly(args.test_only)
        mysource.setnobnodes(args.no_bnodes)

        # run tests first
        if (args.no_verify or args.skip_tests) is not True:
            suite = mysource.getTestSuite()
            if suite is None:
                logger.warn("No tests configured for this source: %s", source)
            else:
                unittest.TextTestRunner(verbosity=2).run(suite)
        else:
            logger.info("Skipping Tests for source: %s", source)

        if args.test_only is False and args.fetch_only is False:
            mysource.parse(args.limit)
            mysource.write(format=args.format)

        # if args.no_verify is not True:

        #    status = mysource.verify()
        #    if status is not True:
        #        logger.error('Source %s did not pass verification tests.', source)
        #        exit(1)
        # else:
        #    logger.info('skipping verification step')
        logger.info('***** Finished with %s *****', source)
    # load configuration parameters
    # for example, keys

    logger.info("All done.")


if __name__ == "__main__":
    main()

###########################
