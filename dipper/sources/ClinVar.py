#! /usr/bin/env python3

"""
    Converts ClinVar XML into
    RDF triples to be ingested by SciGraph.
    These triples conform to the core of the
    SEPIO Evidence & Provenance model

    We also use the clinvar curated gene to disease
    mappings to discern the functional consequence of
    a variant on a gene in cases where this is ambiguous.
    For example, some variants are located in two
    genes overlapping on different strands, and may
    only have a functional consequence on one gene.
    This is suboptimal and we should look for a source
    that directly provides this.

    creating a test set.
        get a full dataset   default ClinVarFullRelease_00-latest.xml.gz
        get the mapping file default gene_condition_source_id
        get a list of RCV    default CV_test_RCV.txt
        put the input files the raw directory
        write the test set back to the raw directory
    ./scripts/ClinVarXML_Subset.sh | gzip > raw/clinvarxml_alpha/ClinVarTestSet.xml.gz


    parsing a test set  (Skolemizing blank nodes  i.e. for Protege)
    dipper/sources/ClinVar.py -f ClinVarTestSet.xml.gz -o ClinVarTestSet_`datestamp`.nt

    For while we are still required to redundantly conflate the owl properties
    in with the data files.

    python3 ./scripts/add-properties2turtle.py --input ./out/ClinVarTestSet_`datestamp`.nt --output ./out/ClinVarTestSet_`datestamp`.nt --format nt

"""

import os
import re
import gzip
import csv
import hashlib
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import List, Dict
import yaml
from dipper.models.ClinVarRecord import ClinVarRecord, Gene,\
    Variant, Allele, Condition, Genotype
from dipper import curie_map

LOG = logging.getLogger(__name__)

# The name of the ingest we are doing
IPATH = re.split(r'/', os.path.realpath(__file__))
(INAME, DOTPY) = re.split(r'\.', IPATH[-1].lower())
RPATH = '/' + '/'.join(IPATH[1:-3])

GLOBAL_TT_PATH = RPATH + '/translationtable/GLOBAL_TERMS.yaml'
LOCAL_TT_PATH = RPATH + '/translationtable/' + INAME + '.yaml'

CV_FTP = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar'

# Global translation table
# Translate labels found in ontologies
# to the terms they are for
GLOBALTT = {}
with open(GLOBAL_TT_PATH) as fh:
    GLOBALTT = yaml.safe_load(fh)

# Local translation table
# Translate external strings found in datasets
# to specific labels found in ontologies
LOCALTT = {}
with open(LOCAL_TT_PATH) as fh:
    LOCALTT = yaml.safe_load(fh)

CURIEMAP = curie_map.get()
CURIEMAP['_'] = 'https://monarchinitiative.org/.well-known/genid/'

# regular expression to limit what is found in the CURIE identifier
# it is ascii centric and may(will) not pass some valid utf8 curies
CURIERE = re.compile(r'^.*:[A-Za-z0-9_][A-Za-z0-9_.]*[A-Za-z0-9_]*$')


def make_spo(sub, prd, obj):
    """
    Decorates the three given strings as a line of ntriples

    """
    # To establish string as a curie and expand,
    # we use a global curie_map(.yaml)
    # sub are always uri  (unless a bnode)
    # prd are always uri (unless prd is 'a')
    # should fail loudly if curie does not exist
    if prd == 'a':
        prd = 'rdf:type'

    try:
        (subcuri, subid) = re.split(r':', sub)
    except Exception:
        LOG.error("not a Subject Curie  '%s'", sub)
        raise ValueError

    try:
        (prdcuri, prdid) = re.split(r':', prd)
    except Exception:
        LOG.error("not a Predicate Curie  '%s'", prd)
        raise ValueError
    objt = ''
    subjt = ''

    # object is a curie or bnode or literal [string|number] NOT None.
    assert (obj is not None), '"None" object for subject ' + sub + ' & pred ' + prd

    objcuri = None
    match = re.match(CURIERE, obj)
    if match is not None:
        try:
            (objcuri, objid) = re.split(r':', obj)
        except ValueError:
            match = None
    if match is not None and objcuri in CURIEMAP:
        objt = CURIEMAP[objcuri] + objid.strip()
        # allow unexpanded bnodes in object
        if objcuri != '_' or CURIEMAP[objcuri] != '_:b':
            objt = '<' + objt + '>'
    elif obj.isdigit():
        objt = '"' + obj + '"^^<http://www.w3.org/2001/XMLSchema#integer>'
    elif obj.isnumeric():
        objt = '"' + obj + '"^^<http://www.w3.org/2001/XMLSchema#double>'
    else:
        # Literals may not contain the characters ", LF, CR '\'
        # except in their escaped forms. internal quotes as well.
        obj = obj.strip('"').replace('\\', '\\\\').replace('"', '\'')
        obj = obj.replace('\n', '\\n').replace('\r', '\\r')
        objt = '"' + obj + '"'

    # allow unexpanded bnodes in subject
    if subcuri is not None and subcuri in CURIEMAP and \
            prdcuri is not None and prdcuri in CURIEMAP:
        subjt = CURIEMAP[subcuri] + subid.strip()
        if subcuri != '_' or CURIEMAP[subcuri] != '_:b':
            subjt = '<' + subjt + '>'
    else:
        raise ValueError(
            "Cant work with: <{}> {} , <{}> {}, {}".format(
                subcuri, subid, prdcuri, prdid, objt))

    return subjt + ' <' + CURIEMAP[prdcuri] + prdid.strip() + '> ' + objt + ' .'


def write_spo(sub, prd, obj, triples):
    """
        write triples to a buffer incase we decide to drop them
    """
    triples.append(make_spo(sub, prd, obj))


def scv_link(scv_sig, rcv_trip):
    '''
    Creates links between SCV based on their pathonicty/significance calls

    # GENO:0000840 - GENO:0000840 --> is_equilavent_to SEPIO:0000098
    # GENO:0000841 - GENO:0000841 --> is_equilavent_to SEPIO:0000098
    # GENO:0000843 - GENO:0000843 --> is_equilavent_to SEPIO:0000098
    # GENO:0000844 - GENO:0000844 --> is_equilavent_to SEPIO:0000098
    # GENO:0000840 - GENO:0000844 --> contradicts SEPIO:0000101
    # GENO:0000841 - GENO:0000844 --> contradicts SEPIO:0000101
    # GENO:0000841 - GENO:0000843 --> contradicts SEPIO:0000101
    # GENO:0000840 - GENO:0000841 --> is_consistent_with SEPIO:0000099
    # GENO:0000843 - GENO:0000844 --> is_consistent_with SEPIO:0000099
    # GENO:0000840 - GENO:0000843 --> strongly_contradicts SEPIO:0000100
    '''

    sig = {  # 'arbitrary scoring scheme increments as powers of two'
        'GENO:0000840': 1,   # pathogenic
        'GENO:0000841': 2,   # likely pathogenic
        'GENO:0000844': 4,   # likely benign
        'GENO:0000843': 8,   # benign
        'GENO:0000845': 16,  # uncertain significance
    }

    lnk = {  # specific result from diff in 'arbitrary scoring scheme'
        0: 'SEPIO:0000098',  # is_equilavent_to
        1: 'SEPIO:0000099',  # is_consistent_with
        2: 'SEPIO:0000101',  # contradicts
        3: 'SEPIO:0000101',  # contradicts
        4: 'SEPIO:0000099',  # is_consistent_with
        6: 'SEPIO:0000101',  # contradicts
        7: 'SEPIO:0000100',  # strongly_contradicts
        8: 'SEPIO:0000126',   # is_inconsistent_with
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


def digest_id(wordage):
    """
    return a deterministic digest of input
    the 'b' is an experiment forcing the first char to be non numeric
    but valid hex; which is in no way required for RDF
    but may help when using the identifier in other contexts
    which do not allow identifiers to begin with a digit

    :param wordage  the string to hash
    :returns 20 hex char digest
    """
    return 'b' + hashlib.sha1(wordage.encode('utf-8')).hexdigest()[1:20]


def process_measure_set(measure_set, rcv_acc) -> Variant:
    """
    Given a MeasureSet, create a Variant object
    :param measure_set: XML object
    :param rcv_acc: str rcv accession
    :return: Variant object
    """
    rcv_variant_id = measure_set.get('ID')
    measure_set_type = measure_set.get('Type')

    # Create Variant object
    rcv_variant_id = 'ClinVarVariant:' + rcv_variant_id
    variant = Variant(id=rcv_variant_id)

    if measure_set_type in ["Haplotype", "Phase unknown", "Distinct chromosomes"]:
        variant.variant_type = measure_set_type
    elif measure_set_type == "Variant":
        # We will attempt to infer the type
        pass
    else:
        raise ValueError(
            rcv_acc + " UNKNOWN VARIANT SUPERTYPE / TYPE \n" + measure_set_type)

    for rcv_measure in measure_set.findall('./Measure'):

        allele_name = rcv_measure.find('./Name/ElementValue[@Type="Preferred"]')
        rcv_allele_label = None
        if allele_name is not None:
            rcv_allele_label = allele_name.text
        # else:
        #    LOG.warning(
        #       rcv_acc + " VARIANT MISSING LABEL")

        allele_type = rcv_measure.get('Type').strip()

        # Create Variant object
        rcv_allele_id = 'ClinVarVariant:' + rcv_measure.get('ID')
        allele = Allele(
            id=rcv_allele_id,
            label=rcv_allele_label,
            variant_type=allele_type
        )

        # this xpath works but is not supported by ElementTree.
        # ./AttributeSet/Attribute[starts-with(@Type, "HGVS")]
        for synonym in rcv_measure.findall('./AttributeSet/Attribute[@Type]'):
            if synonym.get('Type') is not None and \
                    synonym.text is not None and \
                    re.match(r'^HGVS', synonym.get('Type')):
                allele.synonyms.append(synonym.text)

        # XRef[@DB="dbSNP"]/@ID
        for dbsnp in rcv_measure.findall('./XRef[@DB="dbSNP"]'):
            allele.dbsnps.append('dbSNP:' + dbsnp.get('ID'))
            allele.synonyms.append('rs' + dbsnp.get('ID'))

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

        rcv_allele_rels = rcv_measure.findall('./MeasureRelationship')

        if rcv_allele_rels is None:  # try letting them all through
            LOG.info(ET.tostring(rcv_measure).decode('utf-8'))
        else:
            for measure in rcv_allele_rels:
                allele_rel_type = measure.get('Type').strip()
                # if rcv_variant_relationship_type is not None:
                #    LOG.warning(
                #        rcv_acc +
                #        ' rcv_variant_relationship_type ' +
                #        rcv_variant_relationship_type)

                # XRef[@DB="Gene"]/@ID
                ncbigene_id = None
                allele_gene = measure.find('./XRef[@DB="Gene"]')
                if allele_gene is not None:
                    ncbigene_id = allele_gene.get('ID')

                allele.genes.append(Gene(
                    id=ncbigene_id,
                    association_to_allele=allele_rel_type
                ))
        variant.alleles.append(allele)

    # If a variant only has one allele
    # Infer variant type from allele type
    # and allele ID from the variant ID
    if len(variant.alleles) == 1:
        variant.alleles[0].id = variant.id
        variant.variant_type = variant.alleles[0].variant_type

    if variant.variant_type is None:
        raise ValueError("{} Unable to infer type "
                         "from {}".format(rcv_acc, measure_set_type))

    return variant


def resolve(label):
    '''
    composite mapping
    given f(x) and g(x)    here:  GLOBALTT & LOCALTT respectivly
    in order of preference
    return g(f(x))|f(x)|g(x) | x
    TODO consider returning x on fall through
    : return label's mapping

    '''
    term_id = label
    if label is not None and label in LOCALTT:
        term_id = LOCALTT[label]
        if term_id in GLOBALTT:
            term_id = GLOBALTT[term_id]
        else:
            LOG.warning(
                'Local translation but do not have a global term_id for %s', label)
    elif label is not None and label in GLOBALTT:
        term_id = GLOBALTT[label]
    else:
        LOG.error('Do not have any mapping for label: %s', label)
    return term_id


def allele_to_triples(allele, triples) -> None:
    """
    Process allele info such as dbsnp ids and synonyms
    :param allele: Allele
    :param triples: List, Buffer to store the triples
    :return: None
    """

    write_spo(allele.id, 'rdf:type', resolve(allele.variant_type), triples)
    write_spo(allele.id, GLOBALTT['in taxon'], GLOBALTT['Homo sapiens'], triples)
    if allele.label is not None:
        write_spo(allele.id, 'rdfs:label', allele.label, triples)

    # <ClinVarVariant:rcv_variant_id><OWL:hasDbXref><dbSNP:rs>
    #
    # Note that making clinvar variants and dbSNPs equivalent
    # causes clique merge bugs, so best to leave them as xrefs
    # Example: https://www.ncbi.nlm.nih.gov/clinvar/variation/31915/
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/21303/
    for dbsnp_id in allele.dbsnps:
        # sameAs or hasdbxref?
        write_spo(
            allele.id,
            GLOBALTT['database_cross_reference'],
            dbsnp_id,
            triples)

    for syn in allele.synonyms:
        write_spo(allele.id, GLOBALTT['has_exact_synonym'], syn, triples)


def record_to_triples(rcv: ClinVarRecord, triples: List, g2p_map: Dict) -> None:
    """
    Given a ClinVarRecord, adds triples to the triples list

    :param rcv: ClinVarRecord
    :param triples: List, Buffer to store the triples
    :param g2p_map: Gene to phenotype dict
    :return: None
    """
    # For all genotypes variants we add a type, label, and has_taxon human
    write_spo(rcv.genovar.id, 'rdf:type', resolve(rcv.genovar.variant_type), triples)
    write_spo(rcv.genovar.id, GLOBALTT['in taxon'], GLOBALTT['Homo sapiens'], triples)
    if rcv.genovar.label is not None:
        write_spo(rcv.genovar.id, 'rdfs:label', rcv.genovar.label, triples)

    gene_allele = []  # List of two tuples (gene, association_to_allele)

    # Check the type of genovar
    if isinstance(rcv.genovar, Variant):
        if len(rcv.genovar.alleles) > 1:
            for allele in rcv.genovar.alleles:
                write_spo(
                    rcv.genovar.id, GLOBALTT['has_variant_part'], allele.id, triples)

        for allele in rcv.genovar.alleles:
            allele_to_triples(allele, triples)
            for gene in allele.genes:
                gene_allele.append((gene.id, gene.association_to_allele))

        # Hack to determine what relationship to make between a gene and variant.
        # First look at the rcv variant gene relationship type to get the correct
        # curie, but override has_affected_feature in cases where a gene to disease
        # association has not been curated

        # TODO refactor this, the intention is to avoid
        # cases where a variant is mapped to two genes on different strands
        # and we want to connect the correct one
        # see https://github.com/monarch-initiative/monarch-app/issues/1591
        # https://github.com/monarch-initiative/dipper/issues/593
        if len([val[1] for val in gene_allele
                if LOCALTT[val[1]] == 'has_affected_feature']) == len(gene_allele):
            for gene, allele_rel in gene_allele:
                is_affected = True
                if not rcv.significance == GLOBALTT['pathogenic_for_condition'] \
                        and not rcv.significance == \
                        GLOBALTT['likely_pathogenic_for_condition']:
                    is_affected = False
                else:
                    for condition in rcv.conditions:
                        if condition.medgen_id is None \
                                or gene not in g2p_map \
                                or condition.medgen_id not in g2p_map[gene]:
                            is_affected = False
                            break
                if is_affected:
                    write_spo(
                        rcv.genovar.id,
                        resolve(allele_rel),
                        'NCBIGene:' + gene,
                        triples)
                else:
                    write_spo(
                        rcv.genovar.id,
                        GLOBALTT['part_of'],
                        'NCBIGene:' + gene,
                        triples)

        else:
            for allele in rcv.genovar.alleles:
                for gene in allele.genes:
                    write_spo(
                        allele.id,
                        GLOBALTT['part_of'],
                        'NCBIGene:' + gene.id,
                        triples)

    elif isinstance(rcv.genovar, Genotype):
        for variant in rcv.genovar.variants:
            write_spo(
                rcv.genovar.id, GLOBALTT['has_variant_part'], variant.id, triples)

            for allele in variant.alleles:
                allele_to_triples(allele, triples)
                for gene in allele.genes:
                    gene_allele.append((gene.id, gene.association_to_allele))
                    write_spo(
                        allele.id,
                        resolve(gene.association_to_allele),
                        'NCBIGene:' + gene.id,
                        triples)

        # Zygosity if we can infer it from the type
        if rcv.genovar.variant_type == "CompoundHeterozygote":
            write_spo(
                rcv.genovar.id,
                GLOBALTT['has_zygosity'],
                GLOBALTT['compound heterozygous'],
                triples)

        # If all variants are within the same single gene,
        # the genotype affects the gene
        if len([val[1] for val in gene_allele
                if val[1] in ['within single gene', 'variant in gene']
                ]) == len(gene_allele) \
                and len({val[0] for val in gene_allele}) == 1:
            write_spo(
                rcv.genovar.id,
                GLOBALTT['has_affected_feature'],
                'NCBIGene:' + gene_allele[0][0],
                triples)
    else:
        raise ValueError("Invalid type for genovar in rcv {}".format(rcv.id))


def write_review_status_scores():
    """
    Make triples that attach a "star" score to each of ClinVar's review statuses.
    (Stars are basically a 0-4 rating of the review status.)

    Per https://www.ncbi.nlm.nih.gov/clinvar/docs/details/
    Table 1. The review status and assignment of stars( with changes made mid-2015)
    Number of gold stars Description and review statuses

    NO STARS:
    <ReviewStatus> "no assertion criteria provided"
    <ReviewStatus> "no assertion provided"
    No submitter provided an interpretation with assertion criteria (no assertion
    criteria provided), or no interpretation was provided (no assertion provided)

    ONE STAR:
    <ReviewStatus> "criteria provided, single submitter"
    <ReviewStatus> "criteria provided, conflicting interpretations"
    One submitter provided an interpretation with assertion criteria (criteria
    provided, single submitter) or multiple submitters provided assertion criteria
    but there are conflicting interpretations in which case the independent values
    are enumerated for clinical significance (criteria provided, conflicting
    interpretations)

    TWO STARS:
    <ReviewStatus> "criteria provided, multiple submitters, no conflicts"
    Two or more submitters providing assertion criteria provided the same
    interpretation (criteria provided, multiple submitters, no conflicts)

    THREE STARS:
    <ReviewStatus> "reviewed by expert panel"
    reviewed by expert panel

    FOUR STARS:
    <ReviewStatus> "practice guideline"
    practice guideline
    A group wishing to be recognized as an expert panel must first apply to ClinGen
    by completing the form that can be downloaded from our ftp site.

    :param None
    :return: list of triples that attach a "star" score to each of ClinVar's review
    statuses

    """
    triples = []
    status_and_scores = {
        "no assertion criteria provided": '0',
        "no assertion provided": '0',
        "criteria provided, single submitter": '1',
        "criteria provided, conflicting interpretations": '1',
        "criteria provided, multiple submitters, no conflicts": '2',
        "reviewed by expert panel": '3',
        "practice guideline": '4',
    }
    for status, score in status_and_scores.items():
        triples.append(
            make_spo(
                GLOBALTT[status],
                GLOBALTT['has specified numeric value'],
                score))
    return triples


def parse():
    """
    Main function for parsing a clinvar XML release and outputting triples
    """

    files = {
        'f1': {
            'file': 'ClinVarFullRelease_00-latest.xml.gz',
            'url': CV_FTP + '/xml/ClinVarFullRelease_00-latest.xml.gz'
        },
        'f2': {
            'file': 'gene_condition_source_id',
            'url': CV_FTP + 'gene_condition_source_id'
        }
    }

    # handle arguments for IO
    argparser = argparse.ArgumentParser()

    # INPUT
    argparser.add_argument(
        '-f', '--filename', default=files['f1']['file'],
        help="input filename. default: '" + files['f1']['file'] + "'")

    argparser.add_argument(
        '-m', '--mapfile', default=files['f2']['file'],
        help="input g2d mapping file. default: '" + files['f2']['file'] + "'")

    argparser.add_argument(
        '-i', '--inputdir', default=RPATH + '/raw/' + INAME,
        help="path to input file. default: '" + RPATH + '/raw/' + INAME + "'")

    argparser.add_argument(
        '-l', "--localtt",
        help="'spud'\t'potato'   default: " +
        RPATH + '/translationtable/' + INAME + '.yaml')

    argparser.add_argument(
        '-g', "--globaltt",
        help="'potato'\t'PREFIX:p123'   default: " +
        RPATH + '/translationtable/GLOBAL_TERM.yaml')

    # output '/dev/stdout' would be my first choice
    argparser.add_argument(
        '-d', "--destination", default=RPATH + '/out',
        help='directory to write into. default: "' + RPATH + '/out"')

    argparser.add_argument(
        '-o', "--output", default=INAME + '.nt',
        help='file name to write to. default: ' + INAME + '.nt')

    argparser.add_argument(
        '-s', '--skolemize', default=True,
        help='default: True. False keeps plain blank nodes  "_:xxx"')

    args = argparser.parse_args()

    basename = re.sub(r'\.xml.gz$', '', args.filename)
    filename = args.inputdir + '/' + args.filename
    mapfile = args.inputdir + '/' + args.mapfile

    # be sure I/O paths exist
    try:
        os.makedirs(args.inputdir)
    except FileExistsError:
        pass  # no problem

    try:
        os.makedirs(args.destination)
    except FileExistsError:
        pass  # no problem

    # check input exists

    # avoid clobbering existing output until we are finished
    outfile = args.destination + '/TMP_' + args.output + '_PART'
    try:
        os.remove(outfile)
    except FileNotFoundError:
        # no problem
        LOG.info("fresh start for %s", outfile)

    outtmp = open(outfile, 'a')
    output = args.destination + '/' + args.output

    # catch and release input for future study
    reject = args.inputdir + '/' + basename + '_reject.xml'
    # ignore = args.inputdir + '/' + basename + '_ignore.txt'  # unused
    try:
        os.remove(reject)
    except FileNotFoundError:
        # no problem
        LOG.info("fresh start for %s", reject)
    reject = open(reject, 'w')

    # default to /dev/stdout if anything amiss

    # Buffer to store the triples below a MONARCH_association
    # before we decide to whether to keep or not"
    rcvtriples = []

    # Buffer to store non redundant triples between RCV sets
    releasetriple = set()

    # make triples to relate each review status to Clinvar's "score" - 0 to 4 stars
    # releasetriple.update(set(write_review_status_scores()))

    g2pmap = {}
    # this needs to be read first
    with open(mapfile, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")

        next(reader)  # header
        for row in reader:
            if row[0] in g2pmap:
                g2pmap[row[0]].append(row[3])
            else:
                g2pmap[row[0]] = [row[3]]

    # Override default global translation table
    if args.globaltt:
        with open(args.globaltt) as globaltt_fh:
            global GLOBALTT
            GLOBALTT = yaml.safe_load(globaltt_fh)

    # Overide default local translation table
    if args.localtt:
        with open(args.localtt) as localtt_fh:
            global LOCALTT
            LOCALTT = yaml.safe_load(localtt_fh)

    # Overide the given Skolem IRI for our blank nodes
    # with an unresovable alternative.
    if args.skolemize is False:
        global CURIEMAP
        CURIEMAP['_'] = '_:'

    # Seed releasetriple to avoid union with the empty set
    # <MonarchData: + args.output> <a> <owl:Ontology>
    releasetriple.add(make_spo('MonarchData:' + args.output, 'a', 'owl:Ontology'))

    rjct_cnt = tot_cnt = 0

    status_and_scores = {
        "no assertion criteria provided": '0',
        "no assertion provided": '0',
        "criteria provided, single submitter": '1',
        "criteria provided, conflicting interpretations": '1',
        "criteria provided, multiple submitters, no conflicts": '2',
        "reviewed by expert panel": '3',
        "practice guideline": '4',
    }

    #######################################################
    # main loop over xml
    # taken in chunks composed of ClinVarSet stanzas
    with gzip.open(filename, 'rt') as clinvar_fh:
        TREE = ET.iterparse(clinvar_fh)  # w/o specifing events it defaults to 'end'
        for event, element in TREE:
            if element.tag != 'ClinVarSet':
                ReleaseSet = element
                continue
            else:
                ClinVarSet = element
                tot_cnt += 1

            if ClinVarSet.find('RecordStatus').text != 'current':
                LOG.warning(
                    "%s is not current", ClinVarSet.get('ID'))

            RCVAssertion = ClinVarSet.find('./ReferenceClinVarAssertion')
            # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/ClinVarAccession/@Acc
            # 162,466  2016-Mar
            rcv_acc = RCVAssertion.find('./ClinVarAccession').get('Acc')

            # I do not expect we care as we shouldn't keep the RCV.
            if RCVAssertion.find('./RecordStatus').text != 'current':
                LOG.warning(
                    "%s <is not current on>", rcv_acc)  # + rs_dated)

            ClinicalSignificance = RCVAssertion.find(
                './ClinicalSignificance/Description').text
            significance = resolve(ClinicalSignificance)

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

            RCV_ClinicalSignificance = RCVAssertion.find('./ClinicalSignificance')
            if RCV_ClinicalSignificance is not None:
                RCV_ReviewStatus = RCV_ClinicalSignificance.find('./ReviewStatus')
                if RCV_ReviewStatus is not None:
                    rcv_review = RCV_ReviewStatus.text.strip()

            #######################################################################
            # Our Genotype/Subject is a sequence alteration / Variant
            # which apparently was Measured

            # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/MeasureSet/@ID
            # 162,466  2016-Mar
            # 366,566  2017-Mar

            # are now >4 types
            # <GenotypeSet ID="424700" Type="CompoundHeterozygote">
            # <MeasureSet  ID="242681" Type="Variant">
            # <MeasureSet  ID="123456" Type="Haplotype">
            # <Measure     ID="46900"  Type="single nucleotide variant">

            # As of 04/2019
            # Measure is no longer a direct child of ReferenceClinVarAssertion
            # Unless a MeasureSet Type="Variant", both the MeasureSet ID and Measure IDs
            # will be resolvable, eg:
            # https://www.ncbi.nlm.nih.gov/clinvar/variation/431733/
            # https://www.ncbi.nlm.nih.gov/clinvar/variation/425238/

            # If MeasureSet Type == Variant, make the ID the child ID
            # Genotypes can have >1 MeasureSets (Variants)
            # MeasureSets can have >1 Measures (Alleles)
            # Measures (Alleles) can have >1 gene

            RCV_MeasureSet = RCVAssertion.find('./MeasureSet')
            # Note: it is a "set" but have only seen a half dozen with two,
            # all of type:  copy number gain  SO:0001742

            genovar = None  # Union[Genotype, Variant, None]

            if RCV_MeasureSet is None:
                #  201705 introduced GenotypeSet a CompoundHeterozygote
                #  with multiple variants
                RCV_GenotypeSet = RCVAssertion.find('./GenotypeSet')
                genovar = Genotype(
                    id="ClinVarVariant:" + RCV_GenotypeSet.get('ID'),
                    label=RCV_GenotypeSet.find(
                        './Name/ElementValue[@Type="Preferred"]').text,
                    variant_type=RCV_GenotypeSet.get('Type')
                )
                for RCV_MeasureSet in RCV_GenotypeSet.findall('./MeasureSet'):
                    genovar.variants.append(
                        process_measure_set(RCV_MeasureSet, rcv_acc))
            else:
                genovar = process_measure_set(RCV_MeasureSet, rcv_acc)

            # Create ClinVarRecord object
            rcv = ClinVarRecord(
                id=RCVAssertion.get('ID'),
                accession=rcv_acc,
                created=RCVAssertion.get('DateCreated'),
                updated=RCVAssertion.get('DateLastUpdated'),
                genovar=genovar,
                significance=significance
            )

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

                for RCV_Trait in RCV_TraitSet.findall('./Trait[@Type="Disease"]'):
                    has_medgen_id = False
                    rcv_disease_db = None
                    rcv_disease_id = None
                    medgen_id = None
                    disease_label = None

                    RCV_TraitName = RCV_Trait.find(
                        './Name/ElementValue[@Type="Preferred"]')

                    if RCV_TraitName is not None:
                        disease_label = RCV_TraitName.text
                    # else:
                    #    LOG.warning(rcv_acc + " MISSING DISEASE NAME")

                    for RCV_TraitXRef in RCV_Trait.findall('./XRef[@DB="OMIM"]'):
                        rcv_disease_db = RCV_TraitXRef.get('DB')
                        rcv_disease_id = RCV_TraitXRef.get('ID')
                        if rcv_disease_id.startswith('PS'):
                            rcv_disease_db = 'OMIMPS'
                        break

                    # Accept Orphanet if no OMIM
                    if rcv_disease_db is None or rcv_disease_id is None:
                        if rcv_disease_db is not None:
                            break
                        for RCV_TraitXRef in RCV_Trait.findall(
                                './XRef[@DB="Orphanet"]'):
                            rcv_disease_db = 'ORPHA'  # RCV_TraitXRef.get('DB')
                            rcv_disease_id = RCV_TraitXRef.get('ID')
                            break

                    # Always get medgen for g2p mapping file
                    for RCV_TraitXRef in RCV_Trait.findall('./XRef[@DB="MedGen"]'):
                        has_medgen_id = True
                        medgen_id = RCV_TraitXRef.get('ID')
                        break

                    if rcv_disease_db is None and has_medgen_id:
                        # use UMLS prefix instead of MedGen
                        # see https://github.com/monarch-initiative/dipper/issues/874
                        rcv_disease_db = 'UMLS'  # RCV_TraitXRef.get('DB')
                    if rcv_disease_id is None and has_medgen_id:
                        rcv_disease_id = medgen_id

                    # See if there are any leftovers. Possibilities include:
                    # EFO, Gene, Human Phenotype Ontology
                    if rcv_disease_db is None:
                        for RCV_TraitXRef in RCV_Trait.findall('./XRef'):
                            LOG.warning(
                                rcv_acc + " UNKNOWN DISEASE DB:\t" +
                                RCV_TraitXRef.get('DB') + ":" + RCV_TraitXRef.get('ID'))
                            # 82372 MedGen
                            #    58 EFO
                            #     1 Human Phenotype Ontology
                            break

                    rcv.conditions.append(Condition(
                        id=rcv_disease_id,
                        label=disease_label,
                        database=rcv_disease_db,
                        medgen_id=medgen_id
                    ))

            # Check that we have enough info from the RCV
            # to justify parsing the related SCVs
            # check that no members of rcv.genovar are none
            # and that at least one condition has an id and db
            if [1 for member in vars(rcv.genovar) if member is None] \
                    or not [
                        1 for condition in rcv.conditions
                        if condition.id is not None and
                        condition.database is not None]:
                LOG.info('%s is under specified. SKIPPING', rcv_acc)
                rjct_cnt += 1
                # Write this Clinvar set out so we can know what we are missing
                print(
                    # minidom.parseString(
                    #    ET.tostring(
                    #        ClinVarSet)).toprettyxml(
                    #           indent="   "), file=reject)
                    #  too slow. doubles time
                    ET.tostring(ClinVarSet).decode('utf-8'), file=reject)
                ClinVarSet.clear()
                continue

            # start anew
            del rcvtriples[:]

            # At this point we should have a ClinVarRecord object with all
            # necessary data.  Next convert it to triples
            record_to_triples(rcv, rcvtriples, g2pmap)

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

                for condition in rcv.conditions:

                    if condition.database is None:
                        continue

                    rcv_disease_curie = condition.database + ':' + condition.id

                    scv_id = SCV_Assertion.get('ID')
                    monarch_id = digest_id(rcv.id + scv_id + condition.id)
                    monarch_assoc = 'MONARCH:' + monarch_id

                    # if we parsed a review status up above, attach this review status
                    # to this association to allow filtering of RCV by review status
                    if rcv_review is not None:
                        write_spo(
                            monarch_assoc,
                                  GLOBALTT['confidence_score'],
                            status_and_scores[rcv_review],
                            rcvtriples)

                    ClinVarAccession = SCV_Assertion.find('./ClinVarAccession')
                    scv_acc = ClinVarAccession.get('Acc')
                    scv_accver = ClinVarAccession.get('Version')
                    scv_orgid = ClinVarAccession.get('OrgID')
                    # scv_updated = ClinVarAccession.get('DateUpdated')  # not used
                    SCV_SubmissionID = SCV_Assertion.find('./ClinVarSubmissionID')
                    if SCV_SubmissionID is not None:
                        scv_submitter = SCV_SubmissionID.get('submitter')

                    # blank node identifiers
                    _evidence_id = '_:' + digest_id(monarch_id + '_evidence')
                    write_spo(
                        _evidence_id, 'rdfs:label', monarch_id + '_evidence',
                        rcvtriples)

                    _assertion_id = '_:' + digest_id(monarch_id + '_assertion')
                    write_spo(
                        _assertion_id, 'rdfs:label', monarch_id + '_assertion',
                        rcvtriples)

                    #                   TRIPLES
                    # <monarch_assoc><rdf:type><OBAN:association>  .
                    write_spo(monarch_assoc, 'rdf:type', 'OBAN:association', rcvtriples)
                    # <monarch_assoc>
                    #   <OBAN:association_has_subject>
                    #       <ClinVarVariant:rcv_variant_id>
                    write_spo(
                        monarch_assoc, 'OBAN:association_has_subject', rcv.genovar.id,
                        rcvtriples)
                    # <ClinVarVariant:rcv_variant_id><rdfs:label><rcv.variant.label>  .

                    # <monarch_assoc><OBAN:association_has_object><rcv_disease_curi>  .
                    write_spo(
                        monarch_assoc, 'OBAN:association_has_object', rcv_disease_curie,
                        rcvtriples)
                    # <rcv_disease_curi><rdfs:label><rcv_disease_label>  .
                    # medgen might not have a disease label
                    if condition.label is not None:
                        write_spo(
                            rcv_disease_curie, 'rdfs:label', condition.label,
                            rcvtriples)

                    # <monarch_assoc><SEPIO:0000007><:_evidence_id>  .
                    write_spo(
                        monarch_assoc,
                        GLOBALTT['has_supporting_evidence_line'],
                        _evidence_id,
                        rcvtriples)
                    # <monarch_assoc><SEPIO:0000015><:_assertion_id>  .
                    write_spo(
                        monarch_assoc,
                        GLOBALTT['proposition_asserted_in'],
                        _assertion_id,
                        rcvtriples)

                    # <:_evidence_id><rdf:type><ECO:0000000> .
                    write_spo(
                        _evidence_id, 'rdf:type', GLOBALTT['evidence'], rcvtriples)

                    # <:_assertion_id><rdf:type><SEPIO:0000001> .
                    write_spo(
                        _assertion_id, 'rdf:type', GLOBALTT['assertion'], rcvtriples)
                    # <:_assertion_id><rdfs:label><'assertion'>  .
                    write_spo(
                        _assertion_id, 'rdfs:label', 'ClinVarAssertion_' + scv_id,
                        rcvtriples)

                    # <:_assertion_id><SEPIO_0000111><:_evidence_id>
                    write_spo(
                        _assertion_id,
                        GLOBALTT['is_assertion_supported_by_evidence'], _evidence_id,
                        rcvtriples)

                    # <:_assertion_id><dc:identifier><scv_acc + '.' + scv_accver>
                    write_spo(
                        _assertion_id, 'dc:identifier', scv_acc + '.' + scv_accver,
                        rcvtriples)
                    # <:_assertion_id><SEPIO:0000018><ClinVarSubmitters:scv_orgid>  .
                    write_spo(
                        _assertion_id,
                        GLOBALTT['created_by'],
                        'ClinVarSubmitters:' + scv_orgid,
                        rcvtriples)
                    # <ClinVarSubmitters:scv_orgid><rdf:type><foaf:organization>  .
                    write_spo(
                        'ClinVarSubmitters:' + scv_orgid,
                        'rdf:type',
                        'foaf:organization',
                        rcvtriples)
                    # <ClinVarSubmitters:scv_orgid><rdfs:label><scv_submitter>  .
                    write_spo(
                        'ClinVarSubmitters:' + scv_orgid, 'rdfs:label', scv_submitter,
                        rcvtriples)
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
                            SCV_Citation = SCV_AttributeSet.find('./Citation')

                            # <:_assertion_id><SEPIO:0000021><scv_eval_date>  .
                            if scv_eval_date != "None":
                                write_spo(
                                    _assertion_id,
                                    GLOBALTT['Date Created'],
                                    scv_eval_date,
                                    rcvtriples)

                            scv_assert_method = SCV_Attribute.text
                            #  need to be mapped to a <sepio:100...n> curie ????
                            # if scv_assert_method in TT:
                            # scv_assert_id = resolve(scv_assert_method)
                            # _assertion_method_id = '_:' + monarch_id + \
                            #    '_assertionmethod_' + digest_id(scv_assert_method)
                            #
                            # changing to not include context till we have IRI

                            # blank node, would be be nice if these were only made once
                            _assertion_method_id = '_:' + digest_id(
                                scv_assert_method + '_assertionmethod')
                            write_spo(
                                _assertion_method_id, 'rdfs:label',
                                scv_assert_method + '_assertionmethod',
                                rcvtriples)

                            #       TRIPLES   specified_by
                            # <:_assertion_id><SEPIO:0000041><_assertion_method_id>
                            write_spo(
                                _assertion_id, GLOBALTT['is_specified_by'],
                                _assertion_method_id,
                                rcvtriples)

                            # <_assertion_method_id><rdf:type><SEPIO:0000037>
                            write_spo(
                                _assertion_method_id,
                                'rdf:type',
                                GLOBALTT['assertion method'],
                                rcvtriples)

                            # <_assertion_method_id><rdfs:label><scv_assert_method>
                            write_spo(
                                _assertion_method_id, 'rdfs:label', scv_assert_method,
                                rcvtriples)

                            # <_assertion_method_id><ERO:0000480><scv_citation_url>
                            if SCV_Citation is not None:
                                SCV_Citation_URL = SCV_Citation.find('./URL')
                                if SCV_Citation_URL is not None:
                                    write_spo(
                                        _assertion_method_id, GLOBALTT['has_url'],
                                        SCV_Citation_URL.text, rcvtriples)

                    # scv_type = ClinVarAccession.get('Type')  # assert == 'SCV' ?
                    # RecordStatus                             # assert =='current' ?

                    # SCV_ReviewStatus = ClinicalSignificance.find('./ReviewStatus')
                    # if SCV_ReviewStatus is not None:
                    #    scv_review = SCV_ReviewStatus.text

                    # SCV/ClinicalSignificance/Citation/ID
                    # see also:
                    # SCV/ObservedIn/ObservedData/Citation/'ID[@Source="PubMed"]
                    for SCV_Citation in ClinicalSignificance.findall(
                            './Citation/ID[@Source="PubMed"]'):
                        scv_citation_id = SCV_Citation.text
                        #           TRIPLES
                        # has_part -> evidence_has_supporting_reference
                        # <:_evidence_id><SEPIO:0000124><PMID:scv_citation_id>  .
                        write_spo(
                            _evidence_id,
                            GLOBALTT['evidence_has_supporting_reference'],
                            'PMID:' + scv_citation_id,
                            rcvtriples)
                        # <:monarch_assoc><dc:source><PMID:scv_citation_id>
                        write_spo(
                            monarch_assoc, 'dc:source', 'PMID:' + scv_citation_id,
                            rcvtriples)

                        # <PMID:scv_citation_id><rdf:type><IAO:0000013>
                        write_spo(
                            'PMID:' + scv_citation_id,
                            'rdf:type', GLOBALTT['journal article'], rcvtriples)

                        # <PMID:scv_citation_id><SEPIO:0000123><literal>

                    scv_significance = scv_geno = None
                    SCV_Description = ClinicalSignificance.find('./Description')
                    if SCV_Description is not None:
                        scv_significance = SCV_Description.text.strip()
                        scv_geno = resolve(scv_significance)
                        unkwn = 'has_uncertain_significance_for_condition'
                        if scv_geno is not None and \
                                LOCALTT[scv_significance] != unkwn and \
                                scv_significance != 'protective':
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
                                scv_geno,
                                rcvtriples)
                            # <rcv_variant_id><scv_geno><rcv_disease_db:rcv_disease_id>
                            write_spo(
                                genovar.id, scv_geno, rcv_disease_curie, rcvtriples)
                            # <monarch_assoc><oboInOwl:hasdbxref><ClinVar:rcv_acc>  .
                            write_spo(
                                monarch_assoc,
                                GLOBALTT['database_cross_reference'],
                                'ClinVar:' + rcv_acc,
                                rcvtriples)

                            # store association's significance to compare w/sibs
                            pathocalls[monarch_assoc] = scv_geno
                        else:
                            del rcvtriples[:]
                            continue
                    # if we have deleted the triples buffer then
                    # there is no point in continueing  (I don't think)
                    if not rcvtriples:
                        continue
                    # scv_assert_type = SCV_Assertion.find('./Assertion').get('Type')
                    # check scv_assert_type == 'variation to disease'?
                    # /SCV/ObservedIn/ObservedData/Citation/'ID[@Source="PubMed"]
                    for SCV_ObsIn in SCV_Assertion.findall('./ObservedIn'):
                        # /SCV/ObservedIn/Sample
                        # /SCV/ObservedIn/Method
                        for SCV_ObsData in SCV_ObsIn.findall('./ObservedData'):
                            for SCV_Citation in SCV_ObsData.findall('./Citation'):

                                for scv_citation_id in SCV_Citation.findall(
                                        './ID[@Source="PubMed"]'):
                                    # evidence_has_supporting_reference
                                    # see also: SCV/ClinicalSignificance/Citation/ID
                                    # <_evidence_id><SEPIO:0000124><PMID:scv_citation_id>
                                    write_spo(
                                        _evidence_id,
                                        GLOBALTT['evidence_has_supporting_reference'],
                                        'PMID:' + scv_citation_id.text, rcvtriples)
                                    # <PMID:scv_citation_id><rdf:type><IAO:0000013>
                                    write_spo(
                                        'PMID:' + scv_citation_id.text,
                                        'rdf:type', GLOBALTT['journal article'],
                                        rcvtriples)

                                    # <:monarch_assoc><dc:source><PMID:scv_citation_id>
                                    write_spo(
                                        monarch_assoc,
                                        'dc:source',
                                        'PMID:' + scv_citation_id.text, rcvtriples)
                                for scv_pub_comment in SCV_Citation.findall(
                                        './Attribute[@Type="Description"]'):
                                    # <PMID:scv_citation_id><rdf:comment><scv_pub_comment>
                                    write_spo(
                                        'PMID:' + scv_citation_id.text,
                                        'rdf:comment', scv_pub_comment, rcvtriples)
                            # for SCV_Citation in SCV_ObsData.findall('./Citation'):
                            for SCV_Description in SCV_ObsData.findall(
                                    'Attribute[@Type="Description"]'):
                                # <_evidence_id> <dc:description> "description"
                                if SCV_Description.text != 'not provided':
                                    write_spo(
                                        _evidence_id,
                                        'dc:description',
                                        SCV_Description.text,
                                        rcvtriples)

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
                                scv_evidence_type = resolve(SCV_OIMT.text.strip())
                                if scv_evidence_type is None:
                                    LOG.warning(
                                        'No mapping for scv_evidence_type: %s',
                                        SCV_OIMT.text)
                                    continue
                                # blank node
                                _provenance_id = '_:' + digest_id(
                                    _evidence_id + scv_evidence_type)

                                write_spo(
                                    _provenance_id, 'rdfs:label',
                                    _evidence_id + scv_evidence_type, rcvtriples)

                                # TRIPLES
                                # has_provenance -> has_supporting_study
                                # <_evidence_id><SEPIO:0000011><_provenence_id>
                                write_spo(
                                    _evidence_id,
                                    GLOBALTT['has_supporting_activity'],
                                    _provenance_id,
                                    rcvtriples)

                                # <_:provenance_id><rdf:type><scv_evidence_type>
                                write_spo(
                                    _provenance_id, 'rdf:type', scv_evidence_type,
                                    rcvtriples)

                                # <_:provenance_id><rdfs:label><SCV_OIMT.text>
                                write_spo(
                                    _provenance_id, 'rdfs:label', SCV_OIMT.text,
                                    rcvtriples)
                    # End of a SCV (a.k.a. MONARCH association)
            # End of the ClinVarSet.
            # output triples that only are known after processing sibbling records
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
            make_spo('MonarchData:' + args.output, 'owl:versionInfo', rs_dated))
        # not finalized
        # releasetriple.add(
        #     make_spo(
        #        'MonarchData:' + args.output, owl:versionIRI,
        #        'MonarchArchive:' RELEASEDATE + '/ttl/' + args.output'))

        # write all remaining triples out
        print('\n'.join(list(releasetriple)), file=outtmp)
    if rjct_cnt > 0:
        LOG.warning(
            'The %i out of %i records not included are written back to \n%s',
            rjct_cnt, tot_cnt, str(reject))
    outtmp.close()
    reject.close()
    os.replace(outfile, output)

    # If the intermediate file is there it is because of a problem to fix elsewhere
    # try:
    #    os.remove(outfile)
    # except FileNotFoundError:


if __name__ == "__main__":
    parse()
