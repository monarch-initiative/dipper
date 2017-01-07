#!/usr/bin/env python3
import argparse
import json
import yaml
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Collapse results of scrape-impc.py and manual mappings from impc_procedures.yaml
    Note the manual map exists due to some procedures being served as pdf and not
    parsable by our web scraper.  There are also duplicate pages for certain iDs,
    for example:
    {"JAX_LDT_001": "https://www.mousephenotype.org/impress/parameters/159/12"},
    {"JAX_LDT_001": "https://www.mousephenotype.org/impress/protocol/159/12"},
    {"JAX_LDT_001": "https://www.mousephenotype.org/impress/parameters/159"}

    In these cases we take the protocol page
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--yaml', '-y', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')

    args = parser.parse_args()
    output_fh = open(args.output, 'w')
    procedure_list = json.load(open(args.input, 'r'))
    param_map = yaml.load(open(args.yaml, 'r'))

    for procedure_map in procedure_list:
        for code in procedure_map:
            if code not in param_map:
                param_map[code] = procedure_map[code]
            elif procedure_map[code] != param_map[code]:
                if re.search(r'protocol', procedure_map[code]):
                    param_map[code] = procedure_map[code]
                elif re.search(r'protocol', param_map[code]):
                    logger.info("Found dupe, keeping {0} over {1}".
                                format(param_map[code], procedure_map[code]))

    json.dump(param_map, output_fh)
    output_fh.close()

if __name__ == "__main__":
    main()
