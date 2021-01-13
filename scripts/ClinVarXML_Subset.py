#! /usr/bin/env python3

'''
    Isolate subset of ClinVar XML for a TestSet based on Various IDs

        - RCV  ClinVar:RCV000123456
        - VCV  ClinVarVariant:13659


'''
import os
import re
# import yaml
import gzip
import logging
import argparse
import xml.etree.ElementTree as ElementTree


LOG = logging.getLogger(__name__)

# The name of the ingest we are doing
IPATH = re.split(r'/', os.path.realpath(__file__))
(INAME, DOTPY) = re.split(r'\.', IPATH[-1].lower())
RPATH = '/' + '/'.join(IPATH[1:-2])
FTP = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml'
files = {
    'f1': {
        'file': 'ClinVarFullRelease_00-latest.xml.gz',
        'url': '/'.join((FTP, 'ClinVarFullRelease_00-latest.xml.gz'))}
}

# handle arguments for IO
ARGPARSER = argparse.ArgumentParser()

# INPUT
ARGPARSER.add_argument(
    '-f', '--filename', default=files['f1']['file'],
    help="input filename. default: '" + files['f1']['file'] + "'")

ARGPARSER.add_argument(
    '-i', '--inputdir', default=RPATH + '/raw/clinvar',
    help="input path. default: '" + RPATH + '/raw/clinvar' "'")

# OUTPUT
ARGPARSER.add_argument(
    '-d', "--destination", default=RPATH + '/raw/clinvar',
    help='output path. default: "' + RPATH + '/raw/clinvar')

ARGPARSER.add_argument(
    '-o', "--output", default=INAME + '.xml',
    help='file name to write to')

ARGPARSER.add_argument(
    '-t', "--testfile", default=RPATH + '/raw/test/clinvar_curie.txt',
    help='file w/list of CURIEs to capture')

ARGS = ARGPARSER.parse_args()

FILENAME = ARGS.inputdir + '/' + ARGS.filename

OUTPUT = ARGS.destination + '/' + ARGS.output

VARIANT = []    # ClinVarVariant:13659     (not VCV000013659)
DISEASE = []    # OMIM: ...
GENE = []       # NCBIGene:  ...
RCV = []        # ClinVar:RCV000783093

LOG.warning("Reading curies from %s", ARGS.testfile)

with open(ARGS.testfile) as f:
    for line in f:
        line = line.partition('#')[0].strip()  # no comment
        if line != "":
            (prfx, lcl_id) = re.split(r':', line, 2)
            if prfx == 'ClinVar':
                RCV.append(lcl_id.strip())

            elif prfx == 'ClinVarVariant':
                VARIANT.append(lcl_id.strip())

            elif prfx == 'NCBIGene':
                GENE.append(lcl_id.strip())

            elif prfx == 'OMIM':
                DISEASE.append(lcl_id.strip())

LOG.warning('RCV has ' + str(len(RCV)))
LOG.warning('GENE has ' + str(len(GENE)))
LOG.warning('DISEASE has ' + str(len(DISEASE)))
LOG.warning('VARIANT has ' + str(len(VARIANT)))

#######################################################
# main loop over xml
# taken in chunks composed of ClinVarSet stanzas

with gzip.open(FILENAME, 'rt') as fh:
    tree = ElementTree.iterparse(fh, events=('start', 'end'))
    # event, root = iter(tree).next()
    for event, element in tree:
        if event == 'start' and element.tag == 'ReleaseSet':
            ReleaseSet = element
            continue
        elif event == 'end' and element.tag == 'ClinVarSet':
            ClinVarSet = element
        else:
            continue

        # Clinvar sets to keep based on the inclusion of some identifier
        # it could be a ClinVarVariant ID or a RCV or an OMIM or a NCBIGene
        keep = False
        rcv_acc = rcv_variant = None

        # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/MeasureSet/@ID
        # these look like '13659'   (not VCV000013659)
        RCVAssertion = ClinVarSet.find('./ReferenceClinVarAssertion')

        # rcv_variants are not unique to a single ClinVarSet
        RCVAMeaSet = RCVAssertion.find('./MeasureSet')
        if RCVAMeaSet is not None:
            rcv_variant = RCVAMeaSet.get('ID')
        if rcv_variant is not None and rcv_variant in VARIANT:
            keep = True
            continue

        # ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/ClinVarAccession/@Acc
        # these look like 'RCV000783093'
        rcv_acc = RCVAssertion.find('./ClinVarAccession').get('Acc')

        if rcv_acc is not None and rcv_acc in RCV:
            keep = True
            continue

        # Disease
        for RCV_TraitSet in RCVAssertion.findall('TraitSet'):
            for RCV_Trait in RCV_TraitSet.findall('Trait[@Type="Disease"]'):
                for RCV_TraitXRef in RCV_Trait.findall('XRef[@DB="OMIM"]'):
                    rcv_disease_id = RCV_TraitXRef.get('ID')
                    if rcv_disease_id is not None and rcv_disease_id in DISEASE:
                        LOG.warning(
                            'Keeping %s for DISEASE %s', rcv_acc, rcv_disease_id)
                        keep = True
                        continue
        # Gene
        # /RCV/MeasureSet/Measure/MeasureRelationship[@Type]/XRef[@DB="Gene"]/@ID
        RCV_MeasureSet = RCVAssertion.find('./MeasureSet')
        if RCV_MeasureSet is None:
            # GenotypeSet
            pass
        elif not keep:  # no point in searchingif we are keeping already
            for rms in RCV_MeasureSet:
                RCV_Measure = rms.findall('./Measure')
                RCV_MRel = rms.find('./MeasureRelationship')
                if RCV_MRel is not None and not keep:
                    RCV_Gene = RCV_MRel.find('./XRef[@DB="Gene"]')
                    if RCV_Gene is not None:
                        rcv_ncbigene_id = RCV_Gene.get('ID')
                        if rcv_ncbigene_id is not None \
                                and rcv_ncbigene_id in GENE:
                            LOG.warning(rcv_ncbigene_id + ' in GENE')
                            keep = True
                            continue

        # check if we caught anything
        if not keep:
            ClinVarSet.clear()  # clinvar set can not clear itself
            ReleaseSet.remove(ClinVarSet)


print('writing to: ' + OUTPUT)

# for people to look at
with open(OUTPUT, 'wb') as output:
    output.write(ElementTree.tostring(ReleaseSet))

# to feed to dipper
with gzip.open(OUTPUT + '.gz', 'wb') as output:
    output.write(ElementTree.tostring(ReleaseSet))
