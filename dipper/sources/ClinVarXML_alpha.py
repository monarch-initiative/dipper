#! /usr/bin/env python3

'''
    First pass at converting ClinVar XML into
    RDF triples to be ingested by SciGraph.
    These triples comform to the core of the
    SEPIO Evidence & Provenance model 2016 Apr

    creating a test set.
        get a full dataset   default ClinVarFullRelease_00-latest.xml.gz
        get a list of RCV    default CV_test_RCV.txt
        put the input files the raw directory
        write the test set back to the raw directory
    ./scripts/ClinVarXML_Subset.sh | gzip > raw/clinvarxml_alpha/ClinVarTestSet.xml.gz

    parsing a test set  (producing plain blank nodes)
    ./dipper/sources/ClinVarXML_alpha.py -f ClinVarTestSet.xml.gz -o ClinVarTestSet_`datestamp`.nt -s False

    parsing a test set  (Skolemizing blank nodes  i.e. for Protege)
    ./dipper/sources/ClinVarXML_alpha.py -f ClinVarTestSet.xml.gz -o ClinVarTestSet_`datestamp`.nt

    For while we are still required to redundantly conflate the owl properties
    in with the data files.

    python3 ./scripts/add-properties2turtle.py --input ./out/ClinVarTestSet_`datestamp`.nt --output ./out/ClinVarTestSet_`datestamp`.nt --format nt

'''
import yaml
import os
import re
import gzip
import hashlib
import logging
import argparse
import xml.etree.ElementTree as ET

# import Requests
# from dipper.sources.Source import Source
# from dipper import curie_map  # hangs on to stale data?

LOG = logging.getLogger(__name__)

# The name of the ingest we are doing
IPATH = re.split(r'/', os.path.realpath(__file__))
(INAME, DOTPY) = re.split(r'\.', IPATH[-1].lower())
RPATH = '/' + '/'.join(IPATH[1:-3])
files = {
    'f1': {
        'file': 'ClinVarFullRelease_00-latest.xml.gz',
        'url': 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz'}
}

# regular expression to limit what is found in the CURIE identifier
# it is ascii centric and may(will) not pass some valid utf8 curies
CURIERE = re.compile(r'^.*:[A-Za-z0-9_][A-Za-z0-9_.]*[A-Za-z0-9_]*$')

# handle arguments for IO
ARGPARSER = argparse.ArgumentParser()

# INPUT
ARGPARSER.add_argument(
    '-f', '--filename', default=files['f1']['file'],
    help="input filename. default: '" + files['f1']['file'] + "'")

ARGPARSER.add_argument(
    '-i', '--inputdir', default=RPATH + '/raw/' + INAME,
    help="path to input file. default: '" + RPATH + '/raw/' + INAME + "'")


ARGPARSER.add_argument(
    '-l', "--localtt",
    default=RPATH + '/translationtable/' + INAME + '.yaml',
    help="'spud'\t'potato'   default: " +
    RPATH + '/translationtable/' + INAME + '.yaml')

ARGPARSER.add_argument(
    '-g', "--globaltt",
    default=RPATH + '/translationtable/global_terms.yaml',
    help="'potato'\t'PREFIX:p123'   default: " +
    RPATH + '/translationtable/global_term.yaml')

# OUTPUT '/dev/stdout' would be my first choice
ARGPARSER.add_argument(
    '-d', "--destination", default=RPATH + '/out',
    help='directory to write into. default: "' + RPATH + '/out"')

ARGPARSER.add_argument(
    '-o', "--output", default=INAME + '.nt',
    help='file name to write to. default: ' + INAME + '.nt')

ARGPARSER.add_argument(
    '-s', '--skolemize', default=True,
    help='default: True. False keeps plain blank nodes  "_:xxx"')

# TODO validate IO arguments
ARGS = ARGPARSER.parse_args()

BASENAME = re.sub(r'\.xml.gz$', '', ARGS.filename)
FILENAME = ARGS.inputdir + '/' + ARGS.filename

# avoid clobbering existing output until we are finished
OUTFILE = ARGS.destination + '/TMP_' + ARGS.output + '_PART'
try:
    os.remove(OUTFILE)
except FileNotFoundError:
    # no problem
    LOG.info("fresh start for %s", OUTFILE)

OUTTMP = open(OUTFILE, 'a')
OUTPUT = ARGS.destination + '/' + ARGS.output

# catch and release input for future study
REJECT = ARGS.inputdir + '/' + BASENAME + '_REJECT.xml'
IGNORE = ARGS.inputdir + '/' + BASENAME + '_IGNORE.txt'
try:
    os.remove(REJECT)
except FileNotFoundError:
    # no problem
    LOG.info("fresh start for %s", REJECT)
REJECT = open(REJECT, 'a')

# default to /dev/stdout if anything amiss

# CURIEMAP = curie_map.get()
# print("\n".join(CURIEMAP))
# this still fails miserably and returns a copy
# other than the one in this tree that I am updating. i.e.
# print(CURIEMAP['SEPIO']) -> key error
# after it is added to the .yaml

# hardcoding this while my loading from curie_map.yaml is wonky
CURIEMAP = {
    '':     'https://monarchinitiative.org',
    '_':    'https://monarchinitiative.org/.well-known/genid/',
    'MONARCH':  'https://monarchinitiative.org/MONARCH_',
    'MonarchData': 'https://data.monarchinitiative.org/ttl/',
    'dc':   'http://purl.org/dc/elements/1.1/',
    'rdf':  'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
    'foaf': 'http://xmlns.com/foaf/0.1/',
    'owl':  'http://www.w3.org/2002/07/owl#',
    'BFO':  'http://purl.obolibrary.org/obo/BFO_',
    'ECO':  'http://purl.obolibrary.org/obo/ECO_',
    'ERO':  'http://purl.obolibrary.org/obo/ERO_',
    'dbSNP': 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
    'GENO': 'http://purl.obolibrary.org/obo/GENO_',
    'GO':   'http://purl.obolibrary.org/obo/GO_',
    'RO':   'http://purl.obolibrary.org/obo/RO_',
    'MP':   'http://purl.obolibrary.org/obo/MP_',
    'OBAN': 'http://purl.org/oban/',
    'OMIM': 'http://purl.obolibrary.org/obo/OMIM_',
    'OBO':  'http://purl.obolibrary.org/obo/',
    'OIO':  'http://www.geneontology.org/formats/oboInOwl#',
    'IAO':  'http://purl.obolibrary.org/obo/IAO_',
    'Orphanet': 'http://www.orpha.net/ORDO/Orphanet_',
    'MedGen':   'http://www.ncbi.nlm.nih.gov/medgen/',
    'NCBITaxon': 'http://purl.obolibrary.org/obo/NCBITaxon_',
    'NCBIGene': 'http://www.ncbi.nlm.nih.gov/gene/',
    'MmusDv':   'http://purl.obolibrary.org/obo/MmusDv_',
    'SEPIO':    'http://purl.obolibrary.org/obo/SEPIO_',
    'SO':   'http://purl.obolibrary.org/obo/SO_',
    'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/',
    'ClinVarSubmitters': 'http://www.ncbi.nlm.nih.gov/clinvar/submitters/',
    'ClinVarVariant':    'http://www.ncbi.nlm.nih.gov/clinvar/variation/',
    'ClinVar':           'http://www.ncbi.nlm.nih.gov/clinvar/',
}

# def fetch():
# wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

# def parse()


# Buffer to store the triples below a MONARCH_association
# before we decide to whether to keep or not"
rcvtriples = []

# Buffer to store non redundant triples between RCV sets
releasetriple = set()


def make_spo(sub, prd, obj):
    '''
    Decorates the three given strings as a line of ntriples

    '''
    # To establish string as a curie and expand,
    # we use a global curie_map(.yaml)
    # sub are allways uri  (unless a bnode)
    # prd are allways uri (unless prd is 'a')
    # should fail loudly if curie does not exist
    if prd == 'a':
        prd = 'rdf:type'

    (subcuri, subid) = re.split(r':', sub)
    (prdcuri, prdid) = re.split(r':', prd)
    objt = ''

    # object is a curie or bnode or literal [string|number]

    match = re.match(CURIERE, obj)
    objcuri = None
    if match is not None:
        try:
            (objcuri, objid) = re.split(r':', obj)
        except ValueError:
            match = None
    if match is not None and objcuri in CURIEMAP:
        objt = CURIEMAP[objcuri] + objid
        # allow unexpanded bnodes in object
        if objcuri != '_' or CURIEMAP[objcuri] != '_:B':
            objt = '<' + objt + '>'
    elif obj.isnumeric():
        objt = '"' + obj + '"'
    else:
        # Literals may not contain the characters ", LF, CR '\'
        # except in their escaped forms. internal quotes as well.
        obj = obj.strip('"').replace('\\', '\\\\').replace('"', '\'')
        obj = obj.replace('\n', '\\n').replace('\r', '\\r')
        objt = '"' + obj + '"'

    # allow unexpanded bnodes in subject
    if subcuri is not None and subcuri in CURIEMAP and \
            prdcuri is not None and prdcuri in CURIEMAP:
        subjt = CURIEMAP[subcuri] + subid
        if subcuri != '_' or CURIEMAP[subcuri] != '_:B':
            subjt = '<' + subjt + '>'

        return subjt + ' <' + CURIEMAP[prdcuri] + prdid + '> ' + objt + ' .'
    else:
        LOG.error('Cant work with: ', subcuri, subid,  prdcuri, prdid, objt)
        return None


def write_spo(sub, prd, obj):
    '''
        write triples to a buffer incase we decide to drop them
    '''
    rcvtriples.append(make_spo(sub, prd, obj))


def scv_link(scv_sig, rcv_trip):
    '''
    Creates links between SCV based on their pathonicty/significance calls

    # GENO:0000840 - GENO:0000840 --> equivalent_to SEPIO:0000098
    # GENO:0000841 - GENO:0000841 --> equivalent_to SEPIO:0000098
    # GENO:0000843 - GENO:0000843 --> equivalent_to SEPIO:0000098
    # GENO:0000844 - GENO:0000844 --> equivalent_to SEPIO:0000098
    # GENO:0000840 - GENO:0000844 --> inconsistent_with SEPIO:0000101
    # GENO:0000841 - GENO:0000844 --> inconsistent_with SEPIO:0000101
    # GENO:0000841 - GENO:0000843 --> inconsistent_with SEPIO:0000101
    # GENO:0000840 - GENO:0000841 --> consistent_with SEPIO:0000099
    # GENO:0000843 - GENO:0000844 --> consistent_with SEPIO:0000099
    # GENO:0000840 - GENO:0000843 --> contradicts SEPIO:0000100
    '''

    sig = {  # 'arbitrary scoring scheme increments as powers of two'
        'GENO:0000840': 1,   # pathogenic
        'GENO:0000841': 2,   # likely pathogenic
        'GENO:0000844': 4,   # likely benign
        'GENO:0000843': 8,   # benign
        'GENO:0000845': 16,  # uncertain significance
    }

    lnk = {  # specific result from diff in 'arbitrary scoring scheme'
        0: 'SEPIO:0000098',
        1: 'SEPIO:0000099',
        2: 'SEPIO:0000101',
        3: 'SEPIO:0000101',
        4: 'SEPIO:0000099',
        6: 'SEPIO:0000101',
        7: 'SEPIO:0000100',
        8: 'SEPIO:0000126',
        12: 'SEPIO:0000126',
        14: 'SEPIO:0000126',
        15: 'SEPIO:0000126',
    }
    keys = sorted(scv_sig.keys())
    for scv_a in keys:
        scv_av = scv_sig.pop(scv_a)
        for scv_b in scv_sig.keys():
            link = lnk[abs(sig[scv_av] - sig[scv_sig[scv_b]])]
            rcv_trip.append(make_spo(scv_a, link, scv_b))
            rcv_trip.append(make_spo(scv_b, link, scv_a))
    return


# return a deterministic digest of input
# the 'b' is an experiment forcing the first char to be
# non numeric but valid hex
# which is in no way required for RDF
# but can help when using the identifier in other contexts
# which do not allow identifiers to begin with a digit
def digest_id(wordage):
    return 'b' + hashlib.sha1(wordage.encode('utf-8')).hexdigest()[0:15]


# Global translation table
# Translate labels found in ontologies
# to the terms they are for
GTT = {}
with open(ARGS.globaltt) as fh:
    GTT = yaml.safe_load(fh)

# Local translation table
# Translate external strings found in datasets
# to specific labels found in ontologies
LTT = {}
with open(ARGS.localtt) as fh:
    LTT = yaml.safe_load(fh)

# Overide the given Skolem IRI for our blank nodes
# with an unresovable alternative.
if ARGS.skolemize is False:
    CURIEMAP['_'] = '_:'


def resolve(label, local_tt):
    '''
    composite mapping
    given f(x) and g(x)    here:  GTT & LTT respectivly
    return g(x)|g(f(x))|f(x) in order of preference
    TODO consider returning x on fall through
    : return label's mapping

    '''

    if label is not None and label in GTT:
        term_id = GTT[label]
    elif label is not None and label in LTT:
        label = LTT[label]
        if label in GTT:
            term_id = GTT[label]
        else:
            LOG.warning(
                'Translated but do not have a term_id for label: ' + label)
            term_id = label
    else:
        LOG.error('Do not have any mapping for label: ' + label)

        term_id = None
    return term_id


# Seed releasetriple to avoid union with the empty set
# <MonarchData: + ARGS.output> <a> <owl:Ontology>
releasetriple.add(
    make_spo('MonarchData:' + ARGS.output, 'a', 'owl:Ontology'))

rjct_cnt = tot_cnt = 0
#######################################################
# main loop over xml
# taken in chunks composed of ClinVarSet stanzas
with gzip.open(FILENAME, 'rt') as fh:
    TREE = ET.iterparse(fh)  # w/o specifing events it defaults to 'end'
    for event, element in TREE:
        if element.tag != 'ClinVarSet':
            ReleaseSet = element
            continue
        else:
            ClinVarSet = element
            tot_cnt += 1

        if ClinVarSet.find('RecordStatus').text != 'current':
            LOG.warning(
                ClinVarSet.get('ID') + " is not current")

        # collect svc significance calls within a rcv
        pathocalls = {}

        # collect a list of othernames for this variant
        rcv_synonyms = []
        rcv_dbsnps = []

        # There is only one RCV per ClinVarSet
        rcv_variant_id = rcv_variant_type = rcv_variant_label = None
        rcv_disease_db = rcv_disease_id = rcv_disease_label = None
        rcv_disease_curi = rcv_ncbigene_id = rcv_gene_symbol = None

        RCVAssertion = ClinVarSet.find('./ReferenceClinVarAssertion')
        rcv_created = RCVAssertion.get('DateCreated')
        rcv_updated = RCVAssertion.get('DateLastUpdated')
        rcv_id = RCVAssertion.get('ID')
        # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/ClinVarAccession/@Acc
        # 162,466  2016-Mar
        rcv_acc = RCVAssertion.find('./ClinVarAccession').get('Acc')

        # I do not expect we care as we shouldn't keep the RCV.
        if RCVAssertion.find('./RecordStatus').text != 'current':
            LOG.warning(
                rcv_acc + " <is not current on>")  # + rs_dated)

        # # # Child elements
        #
        # /RCV/Assertion
        # /RCV/AttributeSet
        # /RCV/Citation
        # /RCV/ClinVarAccession
        # /RCV/ClinicalSignificance
        # /RCV/MeasureSet
        # /RCV/ObservedIn
        # /RCV/RecordStatus
        # /RCV/TraitSet

        #######################################################################
        # Our Genotype/Subject is a sequence alteration / Variant
        # which apparently was Measured

        # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/MeasureSet/@ID
        # 162,466  2016-Mar
        # 366,566  2017-Mar

        # are now up to three types
        # <GenotypeSet ID="424700" Type="CompoundHeterozygote">
        # <MeasureSet  ID="242681" Type="Variant">
        # <Measure     ID="46900"  Type="single nucleotide variant">

        RCV_MeasureSet = RCVAssertion.find('./MeasureSet')
        # Note: it is a "set" but have only seen a half dozen with two,
        # all of type:  copy number gain  SO:0001742

        if RCV_MeasureSet is None:
            #  201705 introduced GenotypeSet a CompoundHeterozygote
            #  with multiple variants
            RCV_GenotypeSet = RCVAssertion.find('./GenotypeSet')
            rcv_variant_supertype = RCV_GenotypeSet.get('Type')
            for RCV_MeasureSet in RCV_GenotypeSet.findall('./MeasureSet'):
                if rcv_variant_id is not None:
                    rcv_variant_id += ',' + RCV_MeasureSet.get('ID')
                else:
                    rcv_variant_id = RCV_MeasureSet.get('ID')
        else:
            rcv_variant_id = RCV_MeasureSet.get('ID')
            rcv_variant_supertype = RCV_MeasureSet.get('Type')

        for RCV_Measure in RCV_MeasureSet.findall('./Measure'):

            if rcv_variant_supertype == "Variant":
                rcv_variant_type = resolve(RCV_Measure.get('Type'), LTT)
            elif rcv_variant_supertype == "Haplotype":
                rcv_variant_type = resolve('haplotype', LTT)
            elif rcv_variant_supertype == "CompoundHeterozygote":
                rcv_variant_type = resolve(
                    'variant single locus complement', LTT)
                # this resolve('has_zygosity', LTT)
                # resolve('complex heterozygous', LTT)
            # get blessing for this
            elif rcv_variant_supertype == "Phase unknown":
                rcv_variant_type = resolve(RCV_Measure.get('Type'), LTT)
            else:
                rcv_variant_id = None
                LOG.warning(
                    rcv_acc + " UNKNOWN VARIANT SUPERTYPE / TYPE \n" +
                    rcv_variant_supertype + " / " + RCV_Measure.get('Type'))
                continue

            RCV_VariantName = RCV_Measure.find(
                './Name/ElementValue[@Type="Preferred"]')
            if RCV_VariantName is not None:
                rcv_variant_label = RCV_VariantName.text
            # else:
            #    LOG.warning(
            #       rcv_acc + " VARIANT MISSING LABEL")

            # XRef[@DB="dbSNP"]/@ID
            for RCV_dbSNP in \
                    RCV_Measure.findall('./XRef[@DB="dbSNP"]'):

                rcv_dbsnps.append(RCV_dbSNP.get('ID'))

            # this xpath works but is not supported by ElementTree.
            # ./AttributeSet/Attribute[starts-with(@Type, "HGVS")]
            for RCV_Synonym in \
                    RCV_Measure.findall('./AttributeSet/Attribute[@Type]'):
                if RCV_Synonym.get('Type') is not None and \
                        RCV_Synonym.text is not None and \
                        re.match(r'^HGVS', RCV_Synonym.get('Type')):
                    rcv_synonyms.append(RCV_Synonym.text)

            # /RCV/MeasureSet/Measure/Name/ElementValue/[@Type="Preferred"]
            # /RCV/MeasureSet/Measure/MeasureRelationship[@Type]/XRef[@DB="Gene"]/@ID

            # RCV_Variant = RCV_Measure.find(
            #    './MeasureRelationship[@Type="variant in gene"]')

            # 540074 genes overlapped by variant
            # 176970 within single gene
            # 24746 within multiple genes by overlap
            # 5698 asserted, but not computed
            # 439 near gene, upstream
            # 374 variant in gene
            # 54 near gene, downstream

            RCV_Variant = RCV_Measure.find('./MeasureRelationship')

            if RCV_Variant is None:  # try letting them all through
                LOG.info(ET.tostring(RCV_Measure).decode('utf-8'))
            else:
                rcv_variant_relationship_type = RCV_Variant.get('Type')
                # if rcv_variant_relationship_type is not None:
                #    LOG.warning(
                #        rcv_acc +
                #        ' rcv_variant_relationship_type ' +
                #        rcv_variant_relationship_type)

                # XRef[@DB="Gene"]/@ID
                RCV_Gene = RCV_Variant.find('./XRef[@DB="Gene"]')
                if rcv_ncbigene_id is None and RCV_Gene is not None:
                    rcv_ncbigene_id = RCV_Gene.get('ID')
                # elif rcv_ncbigene_id is None:
                #    LOG.warning(rcv_acc + " VARIANT MISSING NCBIGene ID")

                # Symbol/ElementValue[@Type="Preferred"]
                RCV_Symbol = RCV_Variant.find(
                    './Symbol/ElementValue[@Type="Preferred"]')
                if rcv_gene_symbol is None and RCV_Symbol is not None:
                    rcv_gene_symbol = RCV_Symbol.text

                if rcv_gene_symbol is None:
                    LOG.warning(rcv_acc + " VARIANT MISSING Gene Symbol")

        #######################################################################
        # the Object is the Disease, here is called a "trait"
        # reluctantly starting with the RCV disease
        # not the SCV traits as submitted due to time constraints

        for RCV_TraitSet in RCVAssertion.findall('./TraitSet'):
            # /RCV/TraitSet/Trait[@Type="Disease"]/@ID
            # 144,327   2016-Mar

            # /RCV/TraitSet/Trait[@Type="Disease"]/XRef/@DB
            #     29 Human Phenotype Ontology
            #     82 EFO
            #    659 Gene
            #  53218 Orphanet
            #  57356 OMIM
            # 142532 MedGen

            RCV_TraitName = RCV_TraitSet.find(
                './Trait[@Type="Disease"]/Name/ElementValue[@Type="Preferred"]')

            if RCV_TraitName is not None:
                rcv_disease_label = RCV_TraitName.text
            # else:
            #    LOG.warning(rcv_acc + " MISSING DISEASE NAME")

            # Prioritize OMIM
            for RCV_Trait in RCV_TraitSet.findall('./Trait[@Type="Disease"]'):
                if rcv_disease_db is not None:
                    break
                for RCV_TraitXRef in RCV_Trait.findall(
                        './XRef[@DB="OMIM"]'):
                    rcv_disease_db = RCV_TraitXRef.get('DB')
                    rcv_disease_id = RCV_TraitXRef.get('ID')
                    break

            # Accept Orphanet if no OMIM
            if rcv_disease_db is None or rcv_disease_id is None:
                for RCV_Trait in \
                        RCV_TraitSet.findall('./Trait[@Type="Disease"]'):
                    if rcv_disease_db is not None:
                        break
                    for RCV_TraitXRef in RCV_Trait.findall(
                            './XRef[@DB="Orphanet"]'):
                        rcv_disease_db = RCV_TraitXRef.get('DB')
                        rcv_disease_id = RCV_TraitXRef.get('ID')
                        break

            # Otherwise go with MedGen
            if rcv_disease_db is None or rcv_disease_id is None:
                for RCV_Trait in \
                        RCV_TraitSet.findall('./Trait[@Type="Disease"]'):
                    if rcv_disease_db is not None:
                        break
                    for RCV_TraitXRef in RCV_Trait.findall(
                            './XRef[@DB="MedGen"]'):
                        rcv_disease_db = RCV_TraitXRef.get('DB')
                        rcv_disease_id = RCV_TraitXRef.get('ID')
                        break

            # See if there are any leftovers. Possibilities include:
            # EFO, Gene, Human Phenotype Ontology
            if rcv_disease_db is None:
                for RCV_Trait in\
                        RCV_TraitSet.findall('./Trait[@Type="Disease"]'):
                    for RCV_TraitXRef in RCV_Trait.findall('./XRef'):
                        LOG.warning(
                            rcv_acc + " UNKNOWN DISEASE DB:\t" +
                            RCV_TraitXRef.get('DB') + ":" +
                            RCV_TraitXRef.get('ID'))
                        # 82372 MedGen
                        #    58 EFO
                        #     1 Human Phenotype Ontology
                        break

        # Check that we have enough info from the RCV
        # to justify parsing the related SCVs
        if rcv_disease_db is None or rcv_disease_id is None or \
                rcv_disease_label is None or rcv_variant_id is None or \
                rcv_variant_type is None or rcv_variant_label is None:
            LOG.info('%s is under specified. SKIPPING', rcv_acc)
            rjct_cnt += 1
            # Write this Clinvar set out so we can know what we are missing
            print(
                # minidom.parseString(
                #    ET.tostring(
                #        ClinVarSet)).toprettyxml(
                #           indent="   "), file=REJECT)
                #  too slow. doubles time
                ET.tostring(ClinVarSet).decode('utf-8'), file=REJECT)
            ClinVarSet.clear()
            continue

        # start anew
        del rcvtriples[:]

        rcv_disease_curi = rcv_disease_db + ':' + rcv_disease_id
        rcv_variant_id = 'ClinVarVariant:' + rcv_variant_id

        if rcv_ncbigene_id is not None and rcv_ncbigene_id.isnumeric():
            rcv_ncbigene_curi = 'NCBIGene:' + str(rcv_ncbigene_id)
            #           RCV only TRIPLES
            term_id = resolve(rcv_variant_relationship_type, LTT)
            if term_id is not None:
                # <rcv_variant_id> <GENO:0000418> <scv_ncbigene_id>
                write_spo(
                    rcv_variant_id,
                    term_id,
                    rcv_ncbigene_curi)
            else:
                # LOG.warning(
                # 'Check relationship type: ' + rcv_variant_relationship_type)
                continue

            # <scv_ncbigene_id><rdfs:label><scv_gene_symbol>
            if rcv_gene_symbol is not None:
                write_spo(rcv_ncbigene_curi, 'rdfs:label', rcv_gene_symbol)

        #######################################################################
        # Descend into each SCV grouped with the current RCV
        #######################################################################

        # keep a collection of a SCV's associations and patho significance call
        # when this RCV's set is complete, interlink based on patho call

        pathocalls = {}

        for SCV_Assertion in ClinVarSet.findall('./ClinVarAssertion'):

            # /SCV/AdditionalSubmitters
            # /SCV/Assertion
            # /SCV/AttributeSet
            # /SCV/Citation
            # /SCV/ClinVarAccession
            # /SCV/ClinVarSubmissionID
            # /SCV/ClinicalSignificance
            # /SCV/Comment
            # /SCV/CustomAssertionScore
            # /SCV/ExternalID
            # /SCV/MeasureSet
            # /SCV/ObservedIn
            # /SCV/RecordStatus
            # /SCV/StudyDescription
            # /SCV/StudyName
            # /SCV/TraitSet

            # init
            # scv_review = scv_significance = None
            # scv_assertcount += 1

            scv_id = SCV_Assertion.get('ID')
            monarch_id = digest_id(rcv_id + scv_id)
            monarch_assoc = 'MONARCH:' + monarch_id

            ClinVarAccession = SCV_Assertion.find('./ClinVarAccession')
            scv_acc = ClinVarAccession.get('Acc')
            scv_accver = ClinVarAccession.get('Version')
            scv_orgid = ClinVarAccession.get('OrgID')
            scv_updated = ClinVarAccession.get('DateUpdated')
            SCV_SubmissionID = SCV_Assertion.find('./ClinVarSubmissionID')
            if SCV_SubmissionID is not None:
                scv_submitter = SCV_SubmissionID.get('submitter')

            # blank node identifiers
            _evidence_id = '_:' + digest_id(monarch_id + '_evidence')
            write_spo(_evidence_id, 'rdfs:label', monarch_id + '_evidence')

            _assertion_id = '_:' + digest_id(monarch_id + '_assertion')
            write_spo(_assertion_id, 'rdfs:label', monarch_id + '_assertion')

            #                   TRIPLES
            # <monarch_assoc><rdf:type><OBAN:association>  .
            write_spo(monarch_assoc, 'rdf:type', 'OBAN:association')
            # <monarch_assoc>
            #   <OBAN:association_has_subject>
            #       <ClinVarVariant:rcv_variant_id>
            write_spo(
                monarch_assoc, 'OBAN:association_has_subject', rcv_variant_id)
            # <ClinVarVariant:rcv_variant_id><rdfs:label><rcv_variant_label>  .
            write_spo(rcv_variant_id, 'rdfs:label', rcv_variant_label)
            # <ClinVarVariant:rcv_variant_id><rdf:type><rcv_variant_type>  .
            write_spo(rcv_variant_id, 'rdf:type', rcv_variant_type)
            if rcv_variant_supertype == "CompoundHeterozygote":
                write_spo(
                   rcv_variant_id,
                   resolve('has_zygosity', LTT),
                   resolve('complex heterozygous', LTT))

            # <ClinVarVariant:rcv_variant_id><GENO:0000418>

            # RCV/MeasureSet/Measure/AttributeSet/XRef[@DB="dbSNP"]/@ID
            # <ClinVarVariant:rcv_variant_id><OWL:sameAs><dbSNP:rs>
            for rcv_variant_dbsnp_id in rcv_dbsnps:
                write_spo(
                    rcv_variant_id,
                    'OIO:hasdbxref',
                    'dbSNP:' + rcv_variant_dbsnp_id)
            rcv_dbsnps = []
            # <ClinVarVariant:rcv_variant_id><in_taxon><human>
            write_spo(
                rcv_variant_id, resolve('in taxon', LTT),
                resolve('Homo sapiens', LTT))

            # /RCV/MeasureSet/Measure/AttributeSet/Attribute[@Type="HGVS.*"]
            for syn in rcv_synonyms:
                write_spo(rcv_variant_id, 'OIO:hasExactSynonym', syn)
            rcv_synonyms = []
            # <monarch_assoc><OBAN:association_has_object><rcv_disease_curi>  .
            write_spo(
                monarch_assoc, 'OBAN:association_has_object', rcv_disease_curi)
            # <rcv_disease_curi><rdfs:label><rcv_disease_label>  .
            write_spo(rcv_disease_curi, 'rdfs:label', rcv_disease_label)
            # <monarch_assoc><SEPIO:0000007><:_evidence_id>  .
            write_spo(
                monarch_assoc,
                resolve('has_supporting_evidence', LTT), _evidence_id)
            # <monarch_assoc><SEPIO:0000015><:_assertion_id>  .
            write_spo(
                monarch_assoc, resolve('is_asserted_in', LTT), _assertion_id)

            # <:_evidence_id><rdf:type><ECO:0000000> .
            write_spo(_evidence_id, 'rdf:type', resolve('evidence', LTT))

            # <:_assertion_id><rdf:type><SEPIO:0000001> .
            write_spo(_assertion_id, 'rdf:type', resolve('assertion', LTT))
            # <:_assertion_id><rdfs:label><'assertion'>  .
            write_spo(
                _assertion_id, 'rdfs:label', 'ClinVarAssertion_' + scv_id)

            # <:_assertion_id><SEPIO_0000111><:_evidence_id>
            write_spo(
                _assertion_id,
                resolve('is_assertion_supported_by_evidence', LTT),
                _evidence_id)

            # <:_assertion_id><dc:identifier><scv_acc + '.' + scv_accver>
            write_spo(
                _assertion_id, 'dc:identifier', scv_acc + '.' + scv_accver)
            # <:_assertion_id><SEPIO:0000018><ClinVarSubmitters:scv_orgid>  .
            write_spo(
                _assertion_id, resolve('created_by', LTT),
                'ClinVarSubmitters:' + scv_orgid)
            # <ClinVarSubmitters:scv_orgid><rdf:type><foaf:organization>  .
            write_spo(
                'ClinVarSubmitters:' + scv_orgid,
                'rdf:type',
                'foaf:organization')
            # <ClinVarSubmitters:scv_orgid><rdfs:label><scv_submitter>  .
            write_spo(
                'ClinVarSubmitters:' + scv_orgid,
                'rdfs:label',
                scv_submitter)
            ################################################################
            ClinicalSignificance = SCV_Assertion.find('./ClinicalSignificance')
            if ClinicalSignificance is not None:
                scv_eval_date = str(
                    ClinicalSignificance.get('DateLastEvaluated'))

            # bummer. cannot specify xpath parent '..' targeting above .find()
            for SCV_AttributeSet in SCV_Assertion.findall('./AttributeSet'):
                # /SCV/AttributeSet/Attribute[@Type="AssertionMethod"]
                SCV_Attribute = SCV_AttributeSet.find(
                    './Attribute[@Type="AssertionMethod"]')
                if SCV_Attribute is not None:
                    SCV_Citation = SCV_AttributeSet.find(
                        './Citation')

                    # <:_assertion_id><SEPIO:0000021><scv_eval_date>  .
                    if scv_eval_date != "None":
                        write_spo(
                            _assertion_id,
                            resolve('date_created', LTT), scv_eval_date)

                    scv_assert_method = SCV_Attribute.text
                    #  need to be mapped to a <sepio:100...n> curie ????
                    # if scv_assert_method in TT:
                    # scv_assert_id = resolve(scv_assert_method, LTT)
                    # _assertion_method_id = '_:' + monarch_id + \
                    #    '_assertionmethod_' + digest_id(scv_assert_method)
                    #
                    # changing to not include context till we have IRI

                    # blank node, would be be nice if these were only made once
                    _assertion_method_id = '_:' + digest_id(
                        scv_assert_method + '_assertionmethod')
                    write_spo(
                        _assertion_method_id, 'rdfs:label',
                        scv_assert_method + '_assertionmethod')

                    #       TRIPLES   specified_by
                    # <:_assertion_id><SEPIO:0000041><_assertion_method_id>
                    write_spo(
                        _assertion_id,
                        resolve('is_specified_by', LTT), _assertion_method_id)

                    # <_assertion_method_id><rdf:type><SEPIO:0000037>
                    write_spo(
                        _assertion_method_id,
                        'rdf:type', resolve('assertion method', LTT))

                    # <_assertion_method_id><rdfs:label><scv_assert_method>
                    write_spo(
                        _assertion_method_id, 'rdfs:label', scv_assert_method)

                    # <_assertion_method_id><ERO:0000480><scv_citation_url>
                    if SCV_Citation is not None:
                        SCV_Citation_URL = SCV_Citation.find('./URL')
                        if SCV_Citation_URL is not None:
                            write_spo(
                                _assertion_method_id,
                                resolve('has_url', LTT),
                                SCV_Citation_URL.text)

            # scv_type = ClinVarAccession.get('Type')  # assert == 'SCV' ?
            # RecordStatus                             # assert =='current' ?

            # SCV_ReviewStatus = ClinicalSignificance.find('./ReviewStatus')
            # if SCV_ReviewStatus is not None:
            #    scv_review = SCV_ReviewStatus.text

            # SCV/ClinicalSignificance/Citation/ID
            # see also:
            # SCV/ObservedIn/ObservedData/Citation/'ID[@Source="PubMed"]
            for SCV_Citation in \
                    ClinicalSignificance.findall(
                        './Citation/ID[@Source="PubMed"]'):
                scv_citation_id = SCV_Citation.text
                #           TRIPLES
                # has_part -> has_supporting_reference
                # <:_evidence_id><SEPIO:0000124><PMID:scv_citation_id>  .
                write_spo(
                    _evidence_id,
                    resolve('has_supporting_reference', LTT),
                    'PMID:' + scv_citation_id)
                # <:monarch_assoc><dc:source><PMID:scv_citation_id>
                write_spo(
                    monarch_assoc,
                    'dc:source',
                    'PMID:' + scv_citation_id)

                # <PMID:scv_citation_id><rdf:type><IAO:0000013>
                write_spo(
                    'PMID:' + scv_citation_id,
                    'rdf:type',
                    resolve('journal article', LTT))

                # <PMID:scv_citation_id><SEPIO:0000123><literal>

            scv_significance = scv_geno = None
            SCV_Description = ClinicalSignificance.find('./Description')
            if SCV_Description is not None:
                scv_significance = SCV_Description.text
                scv_geno = resolve(scv_significance, LTT)
                if scv_geno is not None \
                        and scv_geno != resolve('uncertain significance', LTT):
                    # we have the association's (SCV) pathnogicty call
                    # and its significance is explicit
                    ##########################################################
                    # 2016 july.
                    # We do not want any of the proceeding triples
                    # unless we get here (no implicit "uncertain significance")
                    # TRIPLES
                    # <monarch_assoc>
                    #   <OBAN:association_has_predicate>
                    #       <scv_geno>
                    write_spo(
                        monarch_assoc,
                        'OBAN:association_has_predicate',
                        scv_geno)
                    # <rcv_variant_id><scv_geno><rcv_disease_db:rcv_disease_id>
                    write_spo(rcv_variant_id, scv_geno, rcv_disease_curi)
                    # <monarch_assoc><OIO:hasdbxref><ClinVar:rcv_acc>  .
                    write_spo(
                        monarch_assoc, 'OIO:hasdbxref', 'ClinVar:' + rcv_acc)

                    # store association's significance to compare w/sibs
                    pathocalls[monarch_assoc] = scv_geno
                else:
                    del rcvtriples[:]
                    continue
            # if we have deleted the triples buffer then
            # there is no point in continueing  (I don't think)
            if len(rcvtriples) == 0:
                continue
            # scv_assert_type = SCV_Assertion.find('./Assertion').get('Type')
            # check scv_assert_type == 'variation to disease'?
            # /SCV/ObservedIn/ObservedData/Citation/'ID[@Source="PubMed"]
            for SCV_ObsIn in SCV_Assertion.findall('./ObservedIn'):
                # /SCV/ObservedIn/Sample
                # /SCV/ObservedIn/Method
                for SCV_ObsData in SCV_ObsIn.findall('./ObservedData'):
                    for SCV_Citation in SCV_ObsData.findall('./Citation'):

                        for scv_citation_id in \
                                SCV_Citation.findall('./ID[@Source="PubMed"]'):
                            # has_supporting_reference
                            # see also: SCV/ClinicalSignificance/Citation/ID
                            # <_evidence_id><SEPIO:0000124><PMID:scv_citation_id>
                            write_spo(
                                _evidence_id,
                                resolve('has_supporting_reference', LTT),
                                'PMID:' + scv_citation_id.text)
                            # <PMID:scv_citation_id><rdf:type><IAO:0000013>
                            write_spo(
                                'PMID:' + scv_citation_id.text,
                                'rdf:type',
                                resolve('journal article', LTT))

                            # <:monarch_assoc><dc:source><PMID:scv_citation_id>
                            write_spo(
                                monarch_assoc,
                                'dc:source',
                                'PMID:' + scv_citation_id.text)
                        for scv_pub_comment in \
                                SCV_Citation.findall(
                                    './Attribute[@Type="Description"]'):
                            # <PMID:scv_citation_id><rdf:comment><scv_pub_comment>
                            write_spo(
                                'PMID:' + scv_citation_id.text,
                                'rdf:comment',
                                scv_pub_comment)
                    # for SCV_Citation in SCV_ObsData.findall('./Citation'):
                    for SCV_Description in \
                            SCV_ObsData.findall(
                                'Attribute[@Type="Description"]'):
                        # <_evidence_id> <dc:description> "description"
                        if SCV_Description.text != 'not provided':
                            write_spo(
                                _evidence_id,
                                'dc:description',
                                SCV_Description.text)

                # /SCV/ObservedIn/TraitSet
                # /SCV/ObservedIn/Citation
                # /SCV/ObservedIn/Co-occurrenceSet
                # /SCV/ObservedIn/Comment
                # /SCV/ObservedIn/XRef

                # /SCV/Sample/Origin
                # /SCV/Sample/Species@TaxonomyId="9606" is a constant
                # scv_affectedstatus = \
                #    SCV_ObsIn.find('./Sample').find('./AffectedStatus').text

                # /SCV/ObservedIn/Method/NamePlatform
                # /SCV/ObservedIn/Method/TypePlatform
                # /SCV/ObservedIn/Method/Description
                # /SCV/ObservedIn/Method/SourceType
                # /SCV/ObservedIn/Method/MethodType
                # /SCV/ObservedIn/Method/MethodType
                for SCV_OIMT in SCV_ObsIn.findall('./Method/MethodType'):
                    if SCV_OIMT.text != 'not provided':
                        scv_evidence_type = resolve(SCV_OIMT.text, LTT)
                        # blank node
                        _provenance_id = '_:' + digest_id(
                            _evidence_id + scv_evidence_type)

                        write_spo(
                            _provenance_id, 'rdfs:label',
                            _evidence_id + scv_evidence_type)

                        # TRIPLES
                        # has_provenance -> has_supporting_study
                        # <_evidence_id><SEPIO:0000011><_provenence_id>
                        write_spo(
                            _evidence_id,
                            resolve('has_supporting_activity', LTT),
                            _provenance_id)

                        # <_:provenance_id><rdf:type><scv_evidence_type>
                        write_spo(
                            _provenance_id, 'rdf:type', scv_evidence_type)

                        # <_:provenance_id><rdfs:label><SCV_OIMT.text>
                        write_spo(
                            _provenance_id, 'rdfs:label', SCV_OIMT.text)
            # End of a SCV (a.k.a. MONARCH association)
        # End of the ClinVarSet.
        # Output triples that only are known after processing sibbling records
        scv_link(pathocalls, rcvtriples)
        # put this RCV's triples in the SET of all triples in this data release
        releasetriple.update(set(rcvtriples))
        del rcvtriples[:]
        ClinVarSet.clear()

    ###############################################################
    # first in is last out
    if ReleaseSet is not None and ReleaseSet.get('Type') != 'full':
        LOG.warning('Not a full release')
    rs_dated = ReleaseSet.get('Dated')  # "2016-03-01 (date_last_seen)
    releasetriple.add(
         make_spo(
            'MonarchData:' + ARGS.output, 'owl:versionInfo', rs_dated))
    # not finalized
    # releasetriple.add(
    #     make_spo(
    #        'MonarchData:' + ARGS.output, owl:versionIRI,
    #        'MonarchArchive:' RELEASEDATE + '/ttl/' + ARGS.output'))

    # write all remaining triples out
    print('\n'.join(list(releasetriple)), file=OUTTMP)
if rjct_cnt > 0:
    LOG.warning(
        'The %i out of %i records not included are written back to \n%s',
        rjct_cnt, tot_cnt, str(REJECT))
OUTTMP.close()
REJECT.close()
os.replace(OUTFILE, OUTPUT)
