#! /usr/bin/python3

'''
    First pass at converting ClinVar XML into
    RDF triples to be ingested by SciGraph.
    These triples comform to the core of the
    SEPIO Evidence & Provenance model 2016 Apr

'''
import os
import re
import sys
import gzip
import hashlib
import logging
import argparse
import xml.etree.ElementTree as ET
# import Requests


LOG = logging.getLogger(__name__)

# from dipper import curie_map  # hangs on to stale data?


# http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/sample_xml/RCV000077146.xml
# I'm running in another dir so you will have to also,
# be sure to have the xml and the mapping file there too
# FILENAME = 'BRCA_ClinVarSet.xml.gz'
# FILENAME = 'ClinVarFullRelease_00-latest.xml.gz'

# The name of the ingest we are doing
IPATH = re.split(r'/', os.path.realpath(__file__))
(INAME, DOTPY) = re.split(r'\.', IPATH[-1].lower())
RPATH = '/' + '/'.join(IPATH[1:-3])
FILES = {'f1': 'ClinVarFullRelease_00-latest.xml.gz'}

# regular expression to limit what is found in the CURIE identifier
# it is ascii centric and may(will) not pass valid utf8 curies
CURIERE = re.compile(r'^.*:[A-Za-z0-9_][A-Za-z0-9_.]*[A-Za-z0-9_]$')

ENIGMA = \
    'https://submit.ncbi.nlm.nih.gov/ft/byid/hxnfuuxx/enigma_rules_2015-03-26.pdf'

# handle arguments for IO
ARGPARSER = argparse.ArgumentParser()

# INPUT
ARGPARSER.add_argument(
    '-f', '--filename', default=FILES['f1'],
    help="path to '" + FILES['f1'] + "'")

ARGPARSER.add_argument(
    '-i', '--inputdir', default=RPATH + '/raw/' + INAME,
    help="path to '" + FILES['f1'] + "'")

ARGPARSER.add_argument(
    '-t', "--transtab",
    default=RPATH + '/translationtable/' + INAME + '.tt',
    help="'pOtatoe'\t'potAtoe'")

# OUTPUT '/dev/stdout' would be my first choice
ARGPARSER.add_argument(
    '-d', "--destination", default=RPATH + '/out',
    help='directory to write into')

ARGPARSER.add_argument(
    '-o', "--output", default=INAME + '.nt',
    help='file name to write to')

# TODO validate IO arguments
ARGS = ARGPARSER.parse_args()

FILENAME = ARGS.inputdir + '/' + ARGS.filename

# CURIEMAP = curie_map.get()
# this still fails miserably and returns a copy
# other than the one in this tree that I am updating. i.e.
# print(CURIEMAP['SEPIO']) -> key error
# after it is added to the .yaml

# hardcoding this while loading from curie_map.yaml is wonky
CURIEMAP = {
    'dc': 'http://purl.org/dc/elements/1.1/',
    'rdf':  'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
    'foaf': 'http://xmlns.com/foaf/0.1/',
    '_':    'https://monarchinitiave.org/.well-known/genid/',
    'BFO':  'http://purl.obolibrary.org/obo/BFO_',
    'ECO':  'http://purl.obolibrary.org/obo/ECO_',
    'ERO':  'http://purl.obolibrary.org/obo/ERO_',
    'GENO': 'http://purl.obolibrary.org/obo/GENO_',
    'GO':   'http://purl.obolibrary.org/obo/GO_',
    'RO':   'http://purl.obolibrary.org/obo/RO_',
    'MP':   'http://purl.obolibrary.org/obo/MP_',
    'OBAN': 'http://purl.org/oban/',
    'OMIM': 'http://purl.obolibrary.org/obo/OMIM_',
    'owl':  'http://www.w3.org/2002/07/owl#',
    'OBO':  'http://purl.obolibrary.org/obo/',
    'OIO':  'http://www.geneontology.org/formats/oboInOwl#',
    'Orphanet': 'http://www.orpha.net/ORDO/Orphanet_',
    'MONARCH':  'http://monarchinitiative.org/MONARCH_',
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


def make_spo(sub, prd, obj):
    '''
    Decorates the three given strings as a line of ntriples

    '''
    # To expand we could use utils/CuriUtil.get_uri(curi)
    # sub & prd are allways uri (unless prd is 'a')
    # should fail loudly if curie does not exist
    if prd == 'a':
        prd = 'rdf:type'

    (subcuri, subid) = re.split(r':', sub)
    (prdcuri, prdid) = re.split(r':', prd)
    objt = ''

    # obj  is a curie or literal [string|number]
    match = re.match(CURIERE, obj)
    objcuri = None
    if match is not None:
        try:
            (objcuri, objid) = re.split(r':', obj,)
        except ValueError:
            match = None
    if match is not None and objcuri in CURIEMAP:
        objt = '<' + CURIEMAP[objcuri] + objid + '>'
    elif obj.isnumeric():
        objt = obj
    else:
        objt = '"' + obj.strip('"') + '"'

    return '<' + CURIEMAP[subcuri] + subid + '> ' + \
           '<' + CURIEMAP[prdcuri] + prdid + '> ' + objt + ' .'


def scv_link(scv_sig):
    '''
    Creates links between SCV based on their pathonnicty significancce

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

    sig = {
        'GENO:0000840': 1,  # pathogenic
        'GENO:0000841': 2,  # likely pathogenic
        'GENO:0000844': 4,  # likely benign
        'GENO:0000843': 8}  # benign

    lnk = {
        0: 'SEPIO:0000098',
        1: 'SEPIO:0000099',
        2: 'SEPIO:0000101',
        3: 'SEPIO:0000101',
        4: 'SEPIO:0000099',
        6: 'SEPIO:0000101',
        7: 'SEPIO:0000100'}
    keys = sorted(scv_sig.keys())
    for scv_a in keys:
        scv_av = scv_sig.pop(scv_a)
        for scv_b in scv_sig.keys():
            link = lnk[abs(sig[scv_av] - sig[scv_sig[scv_b]])]
            print(make_spo(scv_a, link, scv_b))
            print(make_spo(scv_b, link, scv_a))


################################################################
# CONSTANTS once at the beginning (would be better not at all).

#                   TRIPLES
# <SEPIO:0000007><rdfs:label><'has_supporting_evidence'>  .
print(make_spo('SEPIO:0000007', 'rdfs:label', 'has_supporting_evidence'))
# <SEPIO:0000011><rdfs:label><'has_provenance'>  .
print(make_spo('SEPIO:0000011', 'rdfs:label', 'has_provenance'))
# <SEPIO:000001><rdfs:label><'has_agent'>  .
print(make_spo('SEPIO:0000001', 'rdfs:label', 'has_agent'))
# <SEPIO:0000095><rdfs:label><'before_date'>  .
print(make_spo('SEPIO:0000095', 'rdfs:label', 'before_date'))
# <SEPIO:0000041><rdfs:label><'specified_by'>  .
print(make_spo('SEPIO:0000041', 'rdfs:label', 'specified_by'))

# tec check if these or others exist
print(
    make_spo(
        'OBAN:association_has_subject', 'rdf:type', 'owl:ObjectProperty'))
print(
    make_spo(
        'OBAN:association_has_predicate', 'rdf:type', 'owl:ObjectProperty'))
print(
    make_spo('OBAN:association_has_object', 'rdf:type', 'owl:ObjectProperty'))
print(
    make_spo(
        'OBAN:association_has_object_property',
        'rdf:type', 'owl:ObjectProperty'))
print(make_spo('OBO:RO_0002085', 'rdf:type', 'owl:ObjectProperty'))
print(make_spo('OBO:RO_0002162', 'rdf:type', 'owl:ObjectProperty'))
print(make_spo('OBO:RO_0003303', 'rdf:type', 'owl:ObjectProperty'))
print(make_spo('OBO:GENO_0000418', 'rdf:type', 'owl:ObjectProperty'))

# larval stage term mapping file
# will want namespace I expect.
# strips comments and blank lines
# ... OR convert to YAML see: curi_map.yaml
TT = {}
with open(ARGS.transtab) as f:
    for line in f:
        line = line.partition('#')[0].strip()
        if line != "":
            (key, val) = re.split(r'\t+', line, 2)
            TT[key.strip()] = val.strip()

#######################################################
# main loop over xml
with gzip.open(FILENAME, 'rt') as fh:
    TREE = ET.parse(fh)
    ReleaseSet = TREE.getroot()
    if ReleaseSet.get('Type') != 'full':
        LOG.warning('Not a full release')
        sys.exit(-1)

    rs_dated = ReleaseSet.get('Dated')  # "2016-03-01 (date_last_seen)

    for ClinVarSet in ReleaseSet.findall('ClinVarSet[RecordStatus]'):
        if ClinVarSet.find('RecordStatus').text != 'current':
            LOG.warning(
                ClinVarSet.get('ID') + " is not current as of " + rs_dated)
            continue  # or break?

        # collect svc significance calls within a rcv
        pathocalls = {}

        # There is only one RCV per ClinVarSet
        rcv_variant_id = rcv_variant_type = rcv_variant_label = None
        rcv_disease_db = rcv_disease_id = rcv_disease_label = None
        rcv_disease_curi = None

        RCVAssertion = ClinVarSet.find('ReferenceClinVarAssertion')
        rcv_created = RCVAssertion.get('DateCreated')
        rcv_updated = RCVAssertion.get('DateLastUpdated')
        rcv_id = RCVAssertion.get('ID')
        # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/ClinVarAccession/@Acc
        # 162,466  2016-Mar
        rcv_acc = RCVAssertion.find('ClinVarAccession').get('Acc')

        # I do not expect we care as we shouldn't keep the RCV.
        if RCVAssertion.find('RecordStatus').text != 'current':
            LOG.warning(
                rcv_acc + " <is not current on> " + rs_dated)
            continue

        # Child elements
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

        RCV_MeasureSet = RCVAssertion.find('MeasureSet')
        # Note: it is a "set" but have only seen a half dozen with two,
        # all of type:  copy number gain  SO:0001742
        rcv_variant_id = RCV_MeasureSet.get('ID')

        for RCV_Measure in \
                RCV_MeasureSet.findall('Measure'):
            rcv_variant_type = TT.get(RCV_Measure.get('Type'))
            if rcv_variant_type is None:
                LOG.warning(
                    rcv_acc + " UNKNOWN VARIANT TYPE " +
                    RCV_Measure.get('Type').text)
                continue

            RCV_VariantName = RCV_Measure.find(
                'Name/ElementValue[@Type="Preferred"]')
            if RCV_VariantName is not None:
                rcv_variant_label = RCV_VariantName.text
            else:
                LOG.warning(
                    rcv_acc + " VARIANT MISSING LABEL")

        # /RCV/MeasureSet/Measure/Name/ElementValue/[@Type="Preferred"]
        #######################################################################
        # the Object is the Disease, here is called a "trait"
        # reluctantly starting with the RCV disease
        # not the SCV traits as submitted due to time constraints

        for RCV_TraitSet in RCVAssertion.findall('TraitSet'):

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
                'Trait[@Type="Disease"]/Name/ElementValue[@Type="Preferred"]')

            if RCV_TraitName is not None:
                rcv_disease_label = RCV_TraitName.text
                # print(rcv_acc + ' ' + rcv_disease_label)
            else:
                LOG.warning(rcv_acc + " MISSING DISEASE NAME")

            # Prioritize OMIM
            for RCV_Trait in RCV_TraitSet.findall('Trait[@Type="Disease"]'):
                if rcv_disease_db is not None:
                    break
                for RCV_TraitXRef in RCV_Trait.findall('XRef[@DB="OMIM"]'):
                    rcv_disease_db = RCV_TraitXRef.get('DB')
                    rcv_disease_id = RCV_TraitXRef.get('ID')
                    break

            # Accept Orphanet if no OMIM
            if rcv_disease_db is None or rcv_disease_id is None:
                for RCV_Trait in \
                        RCV_TraitSet.findall('Trait[@Type="Disease"]'):
                    if rcv_disease_db is not None:
                        break
                    for RCV_TraitXRef in RCV_Trait.findall(
                            'XRef[@DB="Orphanet"]'):
                        rcv_disease_db = RCV_TraitXRef.get('DB')
                        rcv_disease_id = RCV_TraitXRef.get('ID')
                        break

            # Otherwise go with MedGen
            if rcv_disease_db is None or rcv_disease_id is None:
                for RCV_Trait in \
                        RCV_TraitSet.findall('Trait[@Type="Disease"]'):
                    if rcv_disease_db is not None:
                        break
                    for RCV_TraitXRef in RCV_Trait.findall(
                            'XRef[@DB="MedGen"]'):
                        rcv_disease_db = RCV_TraitXRef.get('DB')
                        rcv_disease_id = RCV_TraitXRef.get('ID')
                        break

            # See if there are any leftovers. Possibilities include:
            # EFO, Gene, Human Phenotype Ontology
            if rcv_disease_db is None:
                for RCV_Trait in\
                        RCV_TraitSet.findall('Trait[@Type="Disease"]'):
                    for RCV_TraitXRef in RCV_Trait.findall('XRef'):
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
            LOG.warning(rcv_acc + " RCV IS WONKY, BYEBYE")
            continue

        rcv_disease_curi = rcv_disease_db + ':' + rcv_disease_id

        #######################################################################
        # Descend into each SCV grouped with the current RCV
        #######################################################################

        # keep a collection of a SCV's associations and patho significance call
        # when this RCV's set is complete, interlink based on patho call

        pathocalls = {}

        for SCV_Assertion in ClinVarSet.findall('ClinVarAssertion'):

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
            monarch_id = hashlib.md5(
                (rcv_id + scv_id).encode('utf-8')).hexdigest()[1:17]
            monarch_assoc = 'MONARCH:' + monarch_id

            ClinVarAccession = SCV_Assertion.find('ClinVarAccession')
            scv_acc = ClinVarAccession.get('Acc')
            scv_accver = ClinVarAccession.get('Version')
            scv_orgid = ClinVarAccession.get('OrgID')
            scv_updated = ClinVarAccession.get('DateUpdated')

            # blank node identifiers
            _evidence_id = '_:' + monarch_id + '_evidence'
            _assertion_id = '_:' + monarch_id + '_assertion'

            #                   TRIPLES
            # <monarch_assoc><rdf:type><OBAN:association>  .
            print(make_spo(monarch_assoc, 'rdf:type', 'OBAN:association'))
            # <monarch_assoc>
            #   <OBAN:association_has_subject>
            #       <ClinVarVariant:rcv_variant_id>
            print(make_spo(
                monarch_assoc,
                'OBAN:association_has_subject',
                'ClinVarVariant:' + rcv_variant_id))
            # <ClinVarVariant:rcv_variant_id><rdfs:label><rcv_variant_label>  .
            print(make_spo(
                'ClinVarVariant:' + rcv_variant_id,
                'rdfs:label', rcv_variant_label))
            # <ClinVarVariant:rcv_variant_id><rdf:type><rcv_variant_type>  .
            print(make_spo(
                'ClinVarVariant:' + rcv_variant_id,
                'rdf:type',
                rcv_variant_type))

            # <ClinVarVariant:rcv_variant_id><GENO:0000418>
            # <ClinVarVariant:rcv_variant_id><rdf:type><owl:Class> TODO ???

            # <monarch_assoc><OBAN:association_has_object><rcv_disease_curi>  .
            print(
                make_spo(
                    monarch_assoc,
                    'OBAN:association_has_object', rcv_disease_curi))
            # <rcv_disease_curi><rdfs:label><rcv_disease_label>  .
            print(make_spo(rcv_disease_curi, 'rdfs:label', rcv_disease_label))
            # <monarch_assoc><SEPIO:0000007><:_evidence_id>  .
            print(
                make_spo(monarch_assoc, 'SEPIO:0000007', _evidence_id))
            # <monarch_assoc><SEPIO:0000015><:_assertion_id>  is asserted in
            print(
                make_spo(
                    monarch_assoc, 'SEPIO:0000015', _assertion_id))

            # <:_evidence_id><rdf:type><SEPIO:0000000> .
            print(make_spo(_evidence_id, 'rdf:type', 'ECO:0000000'))
            # <:_evidence_id><rdfs:label><'evidence line'> .  # nope
            # print(make_spo(_evidence_id, 'rdfs:label', 'evidence line'))
            # <:_assertion_id><rdf:type><SEPIO:0000001> .
            print(make_spo(_assertion_id, 'rdf:type', 'SEPIO:0000001'))
            # <:_assertion_id><rdfs:label><'assertion'>  .
            print(
                make_spo(
                    _assertion_id, 'rdfs:label', scv_id))
            # <:_assertion_id><SEPIO_0000111><:_evidence_id>  is_ass_supprt_by
            print(
                make_spo(
                    _assertion_id, 'SEPIO:0000111', _evidence_id))

            # <:_assertion_id><dc:identifier><scv_acc + '.' + scv_accver>
            print(make_spo(
                _assertion_id,
                'dc:identifier',
                scv_acc + '.' + scv_accver))
            # <:_assertion_id><SEPIO:0000018><ClinVarSubmitters:scv_orgid>  .
            print(make_spo(
                _assertion_id,
                'SEPIO:0000018',
                'ClinVarSubmitters:' + scv_orgid))
            # <ClinVarSubmitters:scv_orgid><rdf:type><foaf:organization>  .
            print(make_spo(
                'ClinVarSubmitters:' + scv_orgid,
                'rdf:type',
                'foaf:organization'))

            # /SCV/AttributeSet/Attribute[@Type="AssertionMethod"]
            SCV_Attribute = \
                SCV_Assertion.find(
                    'AttributeSet/Attribute[@Type="AssertionMethod"]')

            ClinicalSignificance = SCV_Assertion.find('ClinicalSignificance')
            scv_eval_date = ClinicalSignificance.get('DateLastEvaluated')

            # changed from SEPIO:0000105 -> SEPIO:0000105 -> SEPIO:0000021
            # <:_assertion_id><SEPIO:0000021><scv_eval_date>  .
            print(
                make_spo(_assertion_id, 'SEPIO:0000021', str(scv_eval_date)))
            if SCV_Attribute is not None:
                scv_assert_method = SCV_Attribute.text
                # this string needs to be mapped to a <sepio:100...n> curie
                if scv_assert_method in TT:
                    scv_assert_id = TT[scv_assert_method]
                    # TRIPLES   specified_by
                    # <:_assertion_id><SEPIO:0000041><scv_att_id>
                    print(make_spo(
                        _assertion_id, 'SEPIO:0000041', scv_assert_id))
                    # <scv_assert_id> <class?> <SEPIO:0000037>

                    # <scv_att_id><rdf:type><SEPIO:1000001>
                    print(make_spo(
                        scv_assert_id, 'rdf:type', 'SEPIO:1000001'))
                    # <scv_att_id><rdfs:label><scv_assert_method>
                    print(make_spo(
                        scv_assert_id, 'rdfs:label', scv_assert_method))
                    # has_url
                    # <scv_att_id><ERO:0000480><ENIGMA>
                    print(make_spo(
                        scv_assert_id, 'ERO:0000480', ENIGMA))

            # scv_type = ClinVarAccession.get('Type')  # assert == 'SCV' ?
            # RecordStatus                             # assert =='current' ?

            # SCV_ReviewStatus = ClinicalSignificance.find('ReviewStatus')
            # if SCV_ReviewStatus is not None:
            #    scv_review = SCV_ReviewStatus.text

            for SCV_Citation in \
                    ClinicalSignificance.findall('Citation/ID[@Source="PubMed"]'):
                scv_citation_id = SCV_Citation.text
                # TRIPLES
                # <:_evidence_id><BFO:0000051><PMID:scv_citation_id>  .
                print(make_spo(
                    _evidence_id,
                    'BFO:0000051',
                    'PMID:' + scv_citation_id))
                # <PMID:scv_citation_id><rdfs:label><'literature-based study'>
                print(make_spo(
                    'PMID:' + scv_citation_id,
                    'rdfs:label',
                    'literature-based study'))

            scv_significance = scv_geno = None
            SCV_Description = ClinicalSignificance.find('Description')
            if SCV_Description is not None:
                scv_significance = SCV_Description.text
                if scv_significance in TT:
                    scv_geno = TT[scv_significance]
                if scv_geno is not None:
                    # we have the association's (SCV) pathnogicty call
                    # TRIPLES
                    # <monarch_assoc><OBAN:association_has_predicate><scv_geno>
                    print(make_spo(
                        monarch_assoc,
                        'OBAN:association_has_predicate',
                        scv_geno))
                    # <rcv_variant_id><scv_geno><rcv_disease_db:rcv_disease_id>
                    print(make_spo(
                        'ClinVarVariant:' + rcv_variant_id, scv_geno,
                        rcv_disease_curi))
                    # <monarch_assoc><OIO:hasdbxref><ClinVar:rcv_acc>  .
                    print(make_spo(
                        monarch_assoc, 'OIO:hasdbxref', 'ClinVar:' + rcv_acc))

                    # store association's significance to compare w/sibs
                    if scv_geno in (
                            'GENO:0000840', 'GENO:0000841',
                            'GENO:0000844', 'GENO:0000843'):
                        pathocalls[monarch_assoc] = scv_geno

            # scv_assert_type = SCV_Assertion.find('Assertion').get('Type')
            # check scv_assert_type == 'variation to disease'?

            for SCV_ObsIn in SCV_Assertion.findall('ObservedIn'):

                # /SCV/ObservedIn/Sample
                # /SCV/ObservedIn/Method
                # /SCV/ObservedIn/ObservedData
                # /SCV/ObservedIn/TraitSet
                # /SCV/ObservedIn/Citation
                # /SCV/ObservedIn/Co-occurrenceSet
                # /SCV/ObservedIn/Comment
                # /SCV/ObservedIn/XRef

                # /SCV/Sample/Origin
                # /SCV/Sample/Species@TaxonomyId="9606" is a constant
                # scv_affectedstatus = \
                #    SCV_ObsIn.find('Sample').find('AffectedStatus').text

                # /SCV/ObservedIn/Method/NamePlatform
                # /SCV/ObservedIn/Method/TypePlatform
                # /SCV/ObservedIn/Method/Description
                # /SCV/ObservedIn/Method/SourceType
                # /SCV/ObservedIn/Method/MethodType
                # /SCV/ObservedIn/Method/MethodType
                for SCV_OIMT in SCV_ObsIn.findall('Method/MethodType'):
                    if SCV_OIMT.text != 'not provided':
                        scv_evidence_type = TT[SCV_OIMT.text]
                        # blank node
                        _provenance_id = \
                            '_:' + \
                            hashlib.md5(
                                (_evidence_id + scv_evidence_type).
                                encode('utf-8')).hexdigest()[1:17]
                        # TRIPLES
                        # has_provenance
                        # <_evidence_id><SEPIO:0000011><_provenence_id>
                        print(make_spo(
                            _evidence_id, 'SEPIO:0000011', _provenance_id))

                        # <_:provenance_id><rdf:type><scv_evidence_type>
                        print(make_spo(
                            _provenance_id, 'rdf:type', scv_evidence_type))

                        # <_:provenance_id><rdfs:label><SCV_OIMT.text>
                        print(make_spo(
                            _provenance_id, 'rdfs:label',
                            _evidence_id + ' ' + SCV_OIMT.text))

                # /SCV/ObservedIn/ObservedData/Attribute@Type
                # /SCV/ObservedIn/ObservedData/Attribute@integerValue

                # Trait being taken from RCV instead (,for now)
                # for SCV_TraitSet in SCV_ObsIn.findall('TraitSet'):

                    # /SCV/ObservedIn/TraitSet/Comment
                    # /SCV/ObservedIn/TraitSet/Trait
                    # for SCV_Trait in SCV_TraitSet.findall('Trait'):
                    #   # /SCV/ObservedIn/TraitSet/Trait/Name
                    #   # /SCV/ObservedIn/TraitSet/Trait/Symbol
                    #   # /SCV/ObservedIn/TraitSet/Trait/TraitRelationship
                    #   # /SCV/ObservedIn/TraitSet/Trait/XRef

                    #   # init
                    #   scv_traitname = scv_xref = scv_db = None

                    #   scv_finding = SCV_TraitSet.find('Trait').get('Type')
                    #   SCV_TraitName = SCV_Trait.find('Name')
                    #   if SCV_TraitName is None:
                    #       # continue
                    #        break
                    #    SCV_TraitNameEV = SCV_TraitName.find('ElementValue')
                    #    if SCV_TraitNameEV is None:
                    #        # continue
                    #        break
                    #        # if 'Preferred' == SCV_TraitName.get('Type');
                    #    scv_traitname = SCV_TraitNameEV.text

                    #    SCV_XRef = SCV_Trait.find('XRef')
                    #    if SCV_XRef:
                    #        scv_xref = SCV_XRef.get('ID')
                    #        scv_db = SCV_XRef.get('DB')
                    #    # write out

            for SCV_MeasureSet in SCV_Assertion.findall('MeasureSet'):
                for SCV_Measure in SCV_MeasureSet:
                    for SCV_AttributeSet in\
                            SCV_Measure.findall('AttributeSet'):
                        for SCV_Attribute in\
                                SCV_AttributeSet.findall('Attribute'):
                            scv_atttype = SCV_Attribute.get('Type')
                            scv_att = SCV_Attribute.text

                    for SCV_MeasureRelationship in \
                            SCV_Measure.findall('MeasureRelationship'):
                        if SCV_MeasureRelationship.get('Type') == \
                                'variant in gene':
                            # 136,850  all only 'variant in gene'
                            # we can look for what gene
                            SCV_Symbol = SCV_MeasureRelationship.find(
                                'Symbol/ElementName[@Type="Preferred"]')
                            if SCV_Symbol is not None:
                                scv_gene_symbol = SCV_Symbol.text
                        # XRef[@DB="Gene"]/@ID
                        SCV_NCBI = SCV_Measure.find('XRef[@DB="Gene"]')
                        if SCV_NCBI is not None:
                            scv_ncbigene_id = 'NCBIGene:' + SCV_NCBI.get('ID')
                            # TRIPLES
                            # <rcv_variant_id><GENO:0000418><scv_ncbigene_id>
                            print(make_spo(
                                rcv_variant_id,
                                'GENO:0000418',
                                scv_ncbigene_id))
                            # <scv_ncbigene_id><rdfs:label><scv_gene_symbol>
                            print(make_spo(
                                scv_ncbigene_id,
                                'rdfs:label',
                                scv_gene_symbol))

            # for SCV_Citation in SCV_Assertion.findall('Citation'):
            #    # init
            #    scv_citesource = scv_citeid = None
            #    scv_citeabbrev = SCV_Citation.get('Abbrev')
            #    scv_citetype = SCV_Citation.get('Type')
            #    scv_citeabbrev = SCV_Citation.get('Abbrev')
            #    SCV_CiteId = SCV_Citation.find('ID')
            #    if SCV_CiteId:
            #        scv_citesource = SCV_CiteId.get('Source')
            #        scv_citeid = SCV_CiteId.text
        # print("write out any scv links for ", pathocalls)
        scv_link(pathocalls)
