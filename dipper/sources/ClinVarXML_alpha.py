#! /usr/bin/python3
import hashlib
import re
import sys
import gzip
import xml.etree.ElementTree as ET

# from dipper import curie_map  # not there yet

# http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/sample_xml/RCV000077146.xml

# I'm running in another dir so you will have to also
# have the xml and the mapping file there too

# FILENAME = 'BRCA_ClinVarSet.xml.gz'
#
FILENAME = 'ClinVarFullRelease_00-latest.xml.gz'


# scv_assertcount = scv_measurecount = scv_traitcount = scv_citecount = 0
# rs_cvset = 0

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
        rcv_variant_id = rcv_variant_type = None
        rcv_disease_db = rcv_disease_id = None
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
            pass

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
        # all of type:  copy number gain            SO:0001742t.
        rcv_variant_id = rcv_variant_type = rcv_variant_label = None
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

        # RCV_TraitSet = RCVAssertion.find('TraitSet')

        for RCV_TraitSet in RCVAssertion.findall('TraitSet'):

            # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/
            # TraitSet/Trait[@Type="Disease"]/@ID
            # 144,327   2016-Mar

            # /ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/
            # TraitSet/Trait[@Type="Disease"]/XRef/@DB
            #     29 Human Phenotype Ontology
            #     82 EFO
            #    659 Gene
            #  53218 Orphanet
            #  57356 OMIM
            # 142532 MedGen

            RCV_TraitName = \
                    RCV_TraitSet.find(
                        'Trait[@Type="Disease"]/Name/ElementValue[@Type="Preferred"]')

            if RCV_TraitName is not None:
                rcv_disease_label = RCV_TraitName.text
                print("rcv_disease_label: ", rcv_disease_label)
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
            if rcv_disease_db is None:
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
                        print(
                            rcv_acc + "\tUNKNOWN DISEASE DB:\t" +
                            RCV_TraitXRef.get('DB') + "\t" +
                            RCV_TraitXRef.get('ID'), file=sys.stderr)
                        # 82372 MedGen
                        #    58 EFO
                        #     1 Human Phenotype Ontology

        # Check that we have enough info from the RCV
        # to justify parsing the related SCVs
        if rcv_disease_db is None or rcv_disease_id is None or \
                rcv_variant_id is None or rcv_variant_type is None:
            continue

        rcv_disease_curi = rcv_disease_db + '.' + rcv_disease_id
        # else:
        #    print(rcv_acc + "\t" + rcv_variant_id  + "\t" + rcv_variant_type)
        #    print(rcv_acc + "\t" + rcv_disease_db + "\t" + rcv_disease_id)

        #######################################################################
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
                rcv_id.encode('utf-8') + scv_id.encode('utf-8')).hexdigest()[1:17]

            monarch_assoc = 'MONARCH:' + monarch_id

            # blank node identifiers
            evidence_id = r'_' + monarch_assoc + '->evidence'
            provenance_id = r'_' + monarch_assoc + '->provenance'
            assertion_id = provenance_id + '->assertion'

            # TODO: Edge labels done once at the beginning
            # <SEPIO:0000007><rdfs:label><'has_supporting_evidence'>  .
            # <SEPIO:0000011><rdfs:label><'has_provenance'>  .
            # <SEPIO:0000017><rdfs:label><'has_agent'>  .
            # <SEPIO:0000095><rdfs:label><'before_date'>  .
            # <SEPIO:0000041><rdfs:label><'specified_by'>  .

            # TRIPLES
            # <monarch_assoc><rdf:type><OBAN:association>  .
            # <monarch_assoc><association_has_subject><ClinVarVariant:rcv_variant_id>  .
            # <ClinVarVariant:rcv_variant_id><rdfs:label><rcv_variant_label>  .
            # <ClinVarVariant:rcv_variant_id><rdf:type><rcv_variant_type>  .
            # <monarch_assoc><association_has_object><rcv_disease_db:rcv_disease_id>  .
            # <rcv_disease_db:rcv_disease_id><rdfs:label><rcv_disease_label>  .

            # <monarch_assoc><SEPIO:0000007><_:evidence_id>  .
            # <monarch_assoc><SEPIO:0000011><_:provenance_id>  .

            # <_:evidence_id><rdf:type><SEPIO:0000000> .
            # <_:evidence_id><rdfs:label><'evidence line'> .
            # <_:provenance_id><rdf:type><SEPIO:0000003> .
            # <_:provenance_id><rdfs:label><'assertion process'>  .
            # <_:provenance_id><ECO:9000001><_:evidence_id>  .
            # <_:provenance_id><had_output><_:assertion_id> .
            # <_:assertion_id><rdf:type><SEPIO:0000001> .
            # <_:assertion_id><rdfs:label><'assertion'>  .

            # scv_name = SCV_Assertion.get('SubmissionName')
            # ClinVarSubmissionID
            ClinVarAccession = SCV_Assertion.find('ClinVarAccession')
            scv_acc = ClinVarAccession.get('Acc')
            scv_accver = int(ClinVarAccession.get('Version'))
            scv_orgid = ClinVarAccession.get('OrgID')
            scv_updated = ClinVarAccession.get('DateUpdated')

            # TRIPLES
            # TODO CURI for
            # CVS: = 'http://www.ncbi.nlm.nih.gov/clinvar/submitters/'
            #
            # <_:assertion_id><dc:identifier><scv_acc + '.' + scv_accver>
            # <_:provenance_id><SEPIO:0000017><CVS:scv_orgid>  .
            # <CVS:scv_orgid><rdf:type><FOAF:organization>  .
            # <_:provenance_id><SEPIO:0000095><scv_updated>  .

            # /SCV/AttributeSet/Attribute[@Type="AssertionMethod"]
            SCV_Attribute = \
                SCV_Assertion.find(
                    'AttributeSet/Attribute[@Type="AssertionMethod"]')
            if SCV_Attribute is not None:
                scv_att_method = SCV_Attribute.text
                # TRIPLES
                # specified_by
                # <_:provenance_id><SEPIO:0000041><scv_att_method>
                # not done!!!
                # SEPIO:1000001
                # 'ENIGMA BRCA1/2 Classification Criteria (2015)'
                # 'variant classification guideline'
                # SEPIO:0000037 'variant classification guideline' #class label

            # scv_type = ClinVarAccession.get('Type')  # assert == 'SCV' ?
            # AdditionalSubmitters
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
                # <_:evidence_id><SEPIO:0000084><PMID:scv_citation_id>
                # <PMID:scv_citation_id><rdfs:label><'literature-based study'>

            scv_significance = scv_geno = None
            if SCV_Description:
                scv_significance = SCV_Description.text
                scv_geno = onto_map[scv_significance]
                # if scv_geno is not None:
                #   we have the association's pathnogicty call
                #   TRIPLES
                #   CLINVAR: = 'http://www.ncbi.nlm.nih.gov/clinvar/'
                #   <monarch_assoc><OBAN:association_has_predicate><scv_geno>
                #   <rcv_variant_id><scv_geno><rcv_disease_db:rcv_disease_id>
                #   <monarch_assoc><OIO:hasdbxref><CLINVAR:rcv_acc>

            # scv_assert_type = SCV_Assertion.find('Assertion').get('Type')
            # check scv_assert_type == 'variation to disease'?

            for SCV_ObsIn in SCV_Assertion.findall('ObservedIn'):

                # /*/*/*/ObservedIn/Sample
                # /*/*/*/ObservedIn/Method
                # /*/*/*/ObservedIn/ObservedData
                # /*/*/*/ObservedIn/TraitSet
                # /*/*/*/ObservedIn/Citation
                # /*/*/*/ObservedIn/Co-occurrenceSet
                # /*/*/*/ObservedIn/Comment
                # /*/*/*/ObservedIn/XRef

                # Sample/Origin
                # Sample/Species@TaxonomyId="9606" is a constant
                # scv_affectedstatus = \
                #    SCV_ObsIn.find('Sample').find('AffectedStatus').text

                # Method/NamePlatform
                # Method/TypePlatform
                # Method/Description
                # Method/SourceType
                # Method/MethodType
                # SCV/ObservedIn/Method/MethodType
                for SCV_OIMT in SCV_ObsIn.findall('Method/MethodType'):
                    if 'not provided' != SCV_OIMT.text:
                        scv_evidence_type = onto_map[SCV_OIMT.text]
                        # TODO need 'not provided' mapping

                    # TRIPLES
                    # has_supporting_process
                    # <_:evidence_id><SEPIO:0000085><scv_evidence_type>

                # ObservedData/Attribute@Type
                # ObservedData/Attribute@integerValue

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
                                scv_ncbigene_id = 'NCBIGene:' + SCV_NCBI.get('ID')

                            # TRIPLES
                            # <rcv_variant_id><GENO:0000418><scv_ncbigene_id>
                            # <scv_ncbigene_id><rdfs:label><scv_gene_symbol>

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
