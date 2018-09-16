#! /usr/bin/env python3

'''
    Isolate subset of ClinVar XML for a TestSet based on Various IDs

'''
import os
import re
# import yaml
import gzip
import logging
import argparse
import xml.etree.ElementTree as ET


LOG = logging.getLogger(__name__)

# The name of the ingest we are doing
IPATH = re.split(r'/', os.path.realpath(__file__))
(INAME, DOTPY) = re.split(r'\.', IPATH[-1].lower())
RPATH = '/' + '/'.join(IPATH[1:-2])
files = {
    'f1': {
        'file': 'ClinVarFullRelease_00-latest.xml.gz',
        'url': 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz'}
}

# handle arguments for IO
ARGPARSER = argparse.ArgumentParser()

# INPUT
ARGPARSER.add_argument(
    '-f', '--filename', default=files['f1']['file'],
    help="input filename. default: '" + files['f1']['file'] + "'")

ARGPARSER.add_argument(
    '-i', '--inputdir', default=RPATH + '/raw/clinvarxml_alpha',
    help="input path. default: '" + RPATH + '/raw/clinvarxml_alpha' "'")

# OUTPUT
ARGPARSER.add_argument(
    '-d', "--destination", default=RPATH + '/raw/clinvarxml_alpha',
    help='output path. default: "' + RPATH + '/raw/clinvarxml_alpha')

ARGPARSER.add_argument(
    '-o', "--output", default=INAME + '.xml',
    help='file name to write to')

ARGPARSER.add_argument(
    '-x', "--xpath", default=RPATH + '/',
    help='the path to emmit')

ARGS = ARGPARSER.parse_args()

FILENAME = ARGS.inputdir + '/' + ARGS.filename

OUTPUT = ARGS.destination + '/' + ARGS.output

element = []
#LOG.warning('RCV has ' + str(len(RCV)))

# ReleaseSet/ClinVarSet/
# ReferenceClinVarAssertion/MeasureSet

#######################################################
# main loop over xml
# taken in chunks composed of ClinVarSet stanzas

with gzip.open(FILENAME, 'rt') as fh:
    tree = ET.iterparse(fh, events=('start', 'end'))
    # event, root = iter(tree).next()
    for event, element in tree:
        if event == 'start' and element.tag == 'ReleaseSet':
            ReleaseSet = element # there is only one root
        elif event == 'end' and element.tag == 'ClinVarSet':
            ClinVarSet = element
            # Have a ClinVarSet to search based on the inclusion of some xpath
            items = ClinVarSet.findall(ARGS.xpath)
            #print(items)
            # for people to look at
            for item in items:
                print(str(ET.tostring(item)))
            ReleaseSet.remove(ClinVarSet)
        else:
            continue
