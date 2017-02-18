#! /usr/bin/env python3

'''
    Isolate subset of ClinVar XML for a TestSet based on Variant IDs

'''
import os
import re
# import sys
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
    help="path to input file. default: '" + RPATH + '/raw/clinvarxml_alpha' "'")

# OUTPUT
ARGPARSER.add_argument(
    '-d', "--destination", default=RPATH + '/raw/clinvarxml_alpha',
    help='directory to write into. default: "' + RPATH + '/raw/clinvarxml_alpha')

ARGPARSER.add_argument(
    '-o', "--output", default=INAME + '_subset.xml',
    help='file name to write to')

ARGS = ARGPARSER.parse_args()

FILENAME = ARGS.inputdir + '/' + ARGS.filename

OUTPUT = ARGS.destination + '/' + ARGS.output

# these variant IDs are taken from the original tab file parser
# I don't know how they were chosen
VARIANTS = [
    4288, 4289, 4290, 4291, 4297, 5240, 5241, 5242, 5243, 5244, 5245, 5246,
    7105, 8877, 9295, 9296, 9297, 9298, 9449, 10072, 10361, 10382, 12528,
    12529, 12530, 12531, 12532, 14353, 14823, 15872, 17232, 17233, 17234,
    17235, 17236, 17237, 17238, 17239, 17284, 17285, 17286, 17287, 18179,
    18180, 18181, 18343, 18363, 31951, 37123, 38562, 94060, 98004, 98005,
    98006, 98008, 98009, 98194, 98195, 98196, 98197, 98198, 100055, 112885,
    114372, 119244, 128714, 130558, 130559, 130560, 130561, 132146, 132147,
    132148, 144375, 146588, 147536, 147814, 147936, 152976, 156327, 161457,
    162000, 167132]

#######################################################
# main loop over xml
# taken in chunks composed of ClinVarSet stanzas

with gzip.open(FILENAME, 'rt') as fh:
    tree = ET.iterparse(fh, events=('start', 'end'))
    # event, root = iter(tree).next()
    for event, element in tree:
        if event == 'start' and element.tag == 'ReleaseSet':
            ReleaseSet = element
            continue
        elif  event == 'end' and element.tag == 'ClinVarSet':
            ClinVarSet = element
        else:
            continue

        # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/MeasureSet/@ID
        RCVAssertion = ClinVarSet.find('ReferenceClinVarAssertion')
        rcv_variant = int(RCVAssertion.find('MeasureSet').get('ID'))

        if rcv_variant not in VARIANTS:
            ClinVarSet.clear()  # can not clear itself
            ReleaseSet.remove(ClinVarSet)
            # remove('<ClinVarSet/>')
        else:
            print(rcv_variant)

print('writing to: ' + OUTPUT)
# tree.write(OUTPUT, encoding='utf8')

with open(OUTPUT, 'wb') as output:
    output.write(ET.tostring(ReleaseSet))
