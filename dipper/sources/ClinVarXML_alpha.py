#! /usr/bin/python3
import hashlib
import re
import sys
import gzip
import xml.etree.ElementTree as ET

# from dipper import curie_map  # not there yet
# from dipper import CurieUtils

# http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/sample_xml/RCV000077146.xml

# I'm running in another dir so you will have to also,
# be sure to have the xml and the mapping file there too

# FILENAME = 'BRCA_ClinVarSet.xml.gz'
#
FILENAME = 'ClinVarFullRelease_00-latest.xml.gz'


def write_spo(s, p, o):
    # to expand we could use utils/CuriUtil.get_uri(curi)
    match = re.match(r'.*:[!-~]+', o)
    if match is not None:
        # TODO: really need a more rigorous literal/curi check
        o = '<' + o + '>'
    else:
        o = '"' + o + '"'
    print('<' + s + '> <' + p + '> ' + o + '  .')
    return

# Edge labels done once at the beginning
# TRIPLES
# <SEPIO:0000007><rdfs:label><'has_supporting_evidence'>  .
write_spo('SEPIO:0000007', 'rdfs:label', 'has_supporting_evidence')
# <SEPIO:0000011><rdfs:label><'has_provenance'>  .
write_spo('SEPIO:0000011', 'rdfs:label', 'has_provenance')
# <SEPIO:000001><rdfs:label><'has_agent'>  .
write_spo('SEPIO:0000001', 'rdfs:label', 'has_agent')
# <SEPIO:0000095><rdfs:label><'before_date'>  .
write_spo('SEPIO:0000095', 'rdfs:label', 'before_date')
# <SEPIO:0000041><rdfs:label><'specified_by'>  .
write_spo('SEPIO:0000041', 'rdfs:label', 'specified_by')

# larval stage term mapping file
# will want namespace I expect.
# strips comments and blank lines
# ... OR convert to YAML see: curi_map.yaml
onto_map = {}
with open("clinvar_alpha_word_ontology.txt") as f:
    for line in f:
        line = line.partition('#')[0].strip()
        if line != "":
            (key, val) = re.split(r'\t+', line, 2)
            onto_map[key.strip()] = val.strip()

#######################################################
# main loop over xml
with gzip.open(FILENAME, 'rt') as fh:
    tree = ET.parse(fh)
    ReleaseSet = tree.getroot()
    if ReleaseSet.get('Type') != 'full':
        print("Not a full release", file=sys.stderr)
        sys.exit(-1)

    rs_dated = ReleaseSet.get('Dated')  # "2016-03-01 (date_last_seen)

    for ClinVarSet in ReleaseSet.findall('ClinVarSet[RecordStatus]'):
        if ClinVarSet.find('RecordStatus').text != 'current':
            print(
                ClinVarSet.get('ID') + " <is not current as of> " +
                rs_dated + '  .', file=sys.stderr)
            continue  # or break?

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
            print(
                rcv_acc + " <is not current on> " + rs_dated, file=sys.stderr)

        # Child elements
        #
        # /*/*/ReferenceClinVarAssertion/Assertion
        # /*/*/ReferenceClinVarAssertion/AttributeSet
        # /*/*/ReferenceClinVarAssertion/Citation
        # /*/*/ReferenceClinVarAssertion/ClinVarAccession
        # /*/*/ReferenceClinVarAssertion/ClinicalSignificance
        # /*/*/ReferenceClinVarAssertion/MeasureSet
        # /*/*/ReferenceClinVarAssertion/ObservedIn
        # /*/*/ReferenceClinVarAssertion/RecordStatus
        # /*/*/ReferenceClinVarAssertion/TraitSet

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
            rcv_variant_type = onto_map.get(RCV_Measure.get('Type'))
            if rcv_variant_type is None:
                print(
                    rcv_acc + " UNKNOWN VARIANT TYPE " +
                    RCV_Measure.get('Type').text, file=sys.stderr)
                continue

            RCV_VariantName = RCV_Measure.find(
                'Name/ElementValue[@Type="Preferred"]')
            if RCV_VariantName is not None:
                rcv_variant_label = RCV_VariantName.text
            else:
                print(
                    rcv_acc + " VARIANT MISSING LABEL", file=sys.stderr)

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
                print(rcv_acc + " MISSING DISEASE NAME ", file=sys.stderr)

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

            # See if there are any leftovers. Possibilities include:
            # EFO, Gene, Human Phenotype Ontology, MedGen
            if rcv_disease_db is None:
                for RCV_Trait in\
                        RCV_TraitSet.findall('Trait[@Type="Disease"]'):
                    for RCV_TraitXRef in RCV_Trait.findall('XRef'):
                        # print(
                        #    rcv_acc + " UNKNOWN DISEASE DB:\t" +
                        #    RCV_TraitXRef.get('DB') + ":" +
                        #    RCV_TraitXRef.get('ID'), file=sys.stderr)
                        # 82372 MedGen
                        #    58 EFO
                        #     1 Human Phenotype Ontology
                        break

        # Check that we have enough info from the RCV
        # to justify parsing the related SCVs
        if rcv_disease_db is None or rcv_disease_id is None or \
                rcv_disease_label is None or rcv_variant_id is None or \
                rcv_variant_type is None or rcv_variant_label is None:
            print(rcv_acc + " ERROR IS WONKY BYEBYE", file=sys.stderr)
            continue

        rcv_disease_curi = rcv_disease_db + ':' + rcv_disease_id

        #######################################################################
        # Descend into each SCV grouped with the current RCV
        #######################################################################

        for SCV_Assertion in ClinVarSet.findall('ClinVarAssertion'):

            # /*/*/ClinVarAssertion/AdditionalSubmitters
            # /*/*/ClinVarAssertion/Assertion
            # /*/*/ClinVarAssertion/AttributeSet
            # /*/*//ClinVarAssertion/Citation
            # /*/*/ClinVarAssertion/ClinVarAccession
            # /*/*/ClinVarAssertion/ClinVarSubmissionID
            # /*/*/ClinVarAssertion/ClinicalSignificance
            # /*/*/ClinVarAssertion/Comment
            # /*/*/ClinVarAssertion/CustomAssertionScore
            # /*/*/ClinVarAssertion/ExternalID
            # /*/*/ClinVarAssertion/MeasureSet
            # /*/*/ClinVarAssertion/ObservedIn
            # /*/*/ClinVarAssertion/RecordStatus
            # /*/*/ClinVarAssertion/StudyDescription
            # /*/*/ClinVarAssertion/StudyName
            # /*/*/ClinVarAssertion/TraitSet
            # scv_assertcount += 1

            # init
            # scv_review = scv_significance = None

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
            _evidence_id = '_' + monarch_assoc + '|evidence'
            _provenance_id = '_' + monarch_assoc + '|provenance'
            _assertion_id = _provenance_id + '|assertion'

            # TRIPLES
            # <monarch_assoc><rdf:type><OBAN:association>  .
            write_spo(monarch_assoc, 'rdf:type', 'OBAN:association')
            # <monarch_assoc>
            #   <OBAN:association_has_subject>
            #       <ClinVarVariant:rcv_variant_id>
            write_spo(
                monarch_assoc,
                'OBAN:association_has_subject',
                'ClinVarVariant:rcv_variant_id')
            # <ClinVarVariant:rcv_variant_id><rdfs:label><rcv_variant_label>  .
            write_spo(
                'ClinVarVariant:' + rcv_variant_id,
                'rdfs:label', rcv_variant_label)
            # <ClinVarVariant:rcv_variant_id><rdf:type><rcv_variant_type>  .
            write_spo(
                'ClinVarVariant:' + rcv_variant_id,
                'rdf:type',
                rcv_variant_type)
            # <monarch_assoc><OBAN:association_has_object><rcv_disease_curi>  .
            write_spo(
                monarch_assoc, 'OBAN:association_has_object', rcv_disease_curi)
            # <rcv_disease_curi><rdfs:label><rcv_disease_label>  .
            write_spo(rcv_disease_curi, 'rdfs:label', rcv_disease_label)
            # <monarch_assoc><SEPIO:0000007><:_evidence_id>  .
            write_spo(monarch_assoc, 'SEPIO:0000007', ':' + _evidence_id)
            # <monarch_assoc><SEPIO:0000011><:_provenance_id>  .
            write_spo(monarch_assoc, 'SEPIO:0000011', ':' + _provenance_id)
            # <:_evidence_id><rdf:type><SEPIO:0000000> .
            write_spo(':' + _evidence_id, 'rdf:type', 'SEPIO:0000000')
            # <:_evidence_id><rdfs:label><'evidence line'> .
            write_spo(':' + _evidence_id, 'rdfs:label', 'evidence line')
            # <:_provenance_id><rdf:type><SEPIO:0000003> .
            write_spo(':' + _provenance_id, 'rdf:type', 'SEPIO:0000003')
            # <:_provenance_id><rdfs:label><'assertion process'>  .
            write_spo(':' + _provenance_id, 'rdfs:label', 'assertion process')
            # <:_provenance_id><ECO:9000001><:_evidence_id>  .
            write_spo(':' + _provenance_id, 'ECO:9000001', ':' + _evidence_id)
            # <:_provenance_id><had_output><:_assertion_id> .
            write_spo(':' + _provenance_id, 'had_output', ':' + _assertion_id)
            # <:_assertion_id><rdf:type><SEPIO:0000001> .
            write_spo(':' + _assertion_id, 'rdf:type', 'SEPIO:0000001')
            # <:_assertion_id><rdfs:label><'assertion'>  .
            write_spo(':' + _assertion_id, 'rdfs:label', 'assertion')
            # <:_assertion_id><dc:identifier><scv_acc + '.' + scv_accver>
            write_spo(
                ':' + _assertion_id,
                'dc:identifier',
                scv_acc + '.' + scv_accver)
            # <:_provenance_id><SEPIO:0000017><ClinVar:submitters/scv_orgid>  .
            write_spo(
                ':' + _provenance_id,
                'SEPIO:0000017',
                'ClinVar:submitters/' + scv_orgid)
            # <ClinVar:submitters/scv_orgid><rdf:type><FOAF:organization>  .
            write_spo(
                'ClinVar:submitters/' + scv_orgid,
                'rdf:type',
                'FOAF:organization')
            # <:_provenance_id><SEPIO:0000095><scv_updated>  .
            write_spo(':' + _provenance_id, 'SEPIO:0000095', scv_updated)

            # /SCV/AttributeSet/Attribute[@Type="AssertionMethod"]
            SCV_Attribute = \
                SCV_Assertion.find(
                    'AttributeSet/Attribute[@Type="AssertionMethod"]')
            if SCV_Attribute is not None:
                scv_att_method = SCV_Attribute.text
                # TRIPLES
                # specified_by
                # <:_provenance_id><SEPIO:0000041><scv_att_method>
                write_spo(':' + _provenance_id, 'EPIO:0000041', scv_att_method)
                # TODO: not done!!!
                # SEPIO:1000001
                # 'ENIGMA BRCA1/2 Classification Criteria (2015)'
                # 'variant classification guideline'
                # SEPIO:0000037 'variant classification guideline' #class label

            # scv_type = ClinVarAccession.get('Type')  # assert == 'SCV' ?
            # RecordStatus                             # assert =='current' ?

            ClinicalSignificance = SCV_Assertion.find('ClinicalSignificance')
            # scv_eval_date = ClinicalSignificance.get('DateLastEvaluated')
            # SCV_ReviewStatus = ClinicalSignificance.find('ReviewStatus')
            # if SCV_ReviewStatus is not None:
            #    scv_review = SCV_ReviewStatus.text
            SCV_Description = ClinicalSignificance.find('Description')

            SCV_Citation = \
                ClinicalSignificance.find('Citation/ID[@Source="PubMed"]')
            if SCV_Citation is not None:
                scv_citation_id = SCV_Citation.text
                # TRIPLES
                # <:_evidence_id><SEPIO:0000084><PMID:scv_citation_id>  .
                write_spo(
                    ':' + _evidence_id,
                    'SEPIO:0000084',
                    'PMID:' + scv_citation_id)
                # <PMID:scv_citation_id><rdfs:label><'literature-based study'>
                write_spo(
                    'PMID:' + scv_citation_id,
                    'rdfs:label',
                    'literature-based study')

            scv_significance = scv_geno = None
            if SCV_Description:
                scv_significance = SCV_Description.text
                scv_geno = onto_map[scv_significance]
                if scv_geno is not None:
                    # we have the association's pathnogicty call
                    # TRIPLES
                    # <monarch_assoc><OBAN:association_has_predicate><scv_geno>
                    write_spo(
                        monarch_assoc,
                        'OBAN:association_has_predicate',
                        scv_geno)
                    # <rcv_variant_id><scv_geno><rcv_disease_db:rcv_disease_id>
                    write_spo(rcv_variant_id, scv_geno, rcv_disease_curi)
                    # <monarch_assoc><OIO:hasdbxref><ClinVar:rcv_acc>  .
                    write_spo(
                        monarch_assoc, 'OIO:hasdbxref', 'ClinVar:' + rcv_acc)

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

                # Sample/Origin
                # Sample/Species@TaxonomyId="9606" is a constant
                # scv_affectedstatus = \
                #    SCV_ObsIn.find('Sample').find('AffectedStatus').text

                # /SCV/ObservedIn/Method/NamePlatform
                # /SCV/ObservedIn/Method/TypePlatform
                # /SCV/ObservedIn/Method/Description
                # /SCV/ObservedIn/Method/SourceType
                # /SCV/ObservedIn/Method/MethodType
                # /SCV/ObservedIn/Method/MethodType
                for SCV_OIMT in SCV_ObsIn.findall('Method/MethodType'):
                    if 'not provided' != SCV_OIMT.text:
                        scv_evidence_type = onto_map[SCV_OIMT.text]
                        # TODO need 'not provided' mapping? prolly not.

                        # TRIPLES
                        # has_supporting_process
                        # <:_evidence_id><SEPIO:0000085><scv_evidence_type>
                        write_spo(
                            ':' + _evidence_id,
                            'SEPIO:0000085',
                            scv_evidence_type)

                # /SCV/ObservedIn/ObservedData/Attribute@Type
                # /SCV/ObservedIn/ObservedData/Attribute@integerValue

                # Trait being taken from RCV instead (,for now)
                # for SCV_TraitSet in SCV_ObsIn.findall('TraitSet'):

                    # /*/*/*/*/TraitSet/Comment
                    # /*/*/*/*/TraitSet/Trait
                    # for SCV_Trait in SCV_TraitSet.findall('Trait'):
                    #   # /*/*/*/*/*/Trait/Name
                    #   # /*/*/*/*/*/Trait/Symbol
                    #   # /*/*/*/*/*/Trait/TraitRelationship
                    #   # /*/*/*/*/*/Trait/XRef

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

                            SCV_NCBI = SCV_Measure.find('XRef[@DB="Gene"]')
                            if SCV_NCBI is not None:
                                scv_ncbigene_id = '\
                                    NCBIGene:' + SCV_NCBI.get('ID')

                                # TRIPLES
                                # <rcv_variant_id><GENO:0000418><scv_ncbigene_id>
                                write_spo(
                                    rcv_variant_id,
                                    'GENO:0000418',
                                    scv_ncbigene_id)
                                # <scv_ncbigene_id><rdfs:label><scv_gene_symbol>
                                write_spo(
                                    scv_ncbigene_id,
                                    'rdfs:labe',
                                    scv_gene_symbol)

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
