"""
ClinVarUtils: basic 
"""

import xml.etree.ElementTree as ET


def isCurrent(element):
    """Determine if the indicated clinvar set is current"""
    rr = element.find("RecordStatus")
    if rr == None:
        return False
    else:
        return(rr.text == "current")

def textIfPresent(element, field):
    """Return the text associated with a field under the element, or
    None if the field is not present"""
    ff = element.find(field)
    if ff == None:
        return None
    else:
        return(ff.text)

class genomicCoordinates:
    """Contains the genomic information on the variant"""

    def __init__(self, element, useNone=False, debug=False):
        if debug:
            print("Parsing genomic coordinates")
        if useNone:
            self.element = None
            self.chrom = None
            self.start = None
            self.stop = None
            self.length = None
            self.referenceAllele = None
            self.alternateAllele = None
        else:
            self.element = element
            self.chrom = element.get("Chr")
            self.start = element.get("start")
            self.stop = element.get("stop")
            self.length = element.get("variantLength")
            self.referenceAllele = element.get("referenceAllele")
            self.alternateAllele = element.get("alternateAllele")



class variant:
    """The Measure set.  We are interested in the variants specifically, 
    but measure sets can be other things as well, such as haplotypes"""
    
    def __init__(self, element, name, id, debug=False):
        self.element = element
        self.id = id
        if debug:
            print("Parsing variant", self.id)
        self.name = name
        self.attribute = dict()
        attrs = element.find("AttributeSet")
        if attrs != None:
            for attrib in attrs.findall("Attribute"):
                self.attribute[attrib.get("Type")] = attrib.text
        self.coordinates = dict()
        for item in element.findall("SequenceLocation"):
            assembly = item.get("Assembly")
            genomic = genomicCoordinates(item, debug=debug)
            self.coordinates[assembly] = genomic
        self.geneSymbol = None
        measureRelationship = element.find("MeasureRelationship")
        if measureRelationship != None:
            symbol = measureRelationship.find("Symbol")
            if symbol != None:
                self.geneSymbol = textIfPresent(symbol, "ElementValue")
            

class referenceAssertion:
    """For gathering the reference assertion"""
    
    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ReferenceClinVarAssertion", self.id)
        self.reviewStatus = None
        self.clinicalSignificance = None
        cs = element.find("ClinicalSignificance")
        if cs != None:
            self.reviewStatus = textIfPresent(cs, "ReviewStatus")
            self.clinicalSignificance = textIfPresent(cs, "Description")
        obs = element.find("ObservedIn")
        if obs == None:
            self.origin = None
            self.ethnicity = None
            self.geographicOrigin = None
            self.age = None
            self.gender = None
            self.familyData = None
            self.method = None
        else:
            sample = obs.find("Sample")
            if sample != None:
                self.origin = textIfPresent(sample, "Origin")
                self.ethnicity = textIfPresent(sample, "Ethnicity")
                self.geographicOrigin = textIfPresent(sample, "GeographicOrigin")
                self.age = textIfPresent(sample, "Age")
                self.gender = textIfPresent(sample, "Gender")
                self.familyData = textIfPresent(sample, "FamilyData")
            method = obs.find("Method")
            if method != None:
                self.method = textIfPresent(method, "MethodType")
        self.variant = None
        measureSet = element.find("MeasureSet")
        #if measureSet.get("Type") == "Variant":
        if debug:
            if len(measureSet.findall("Measure")) > 1:
                print(self.id, "has multiple measures")
        if len(measureSet.findall("Measure")) == 1:
            name = measureSet.find("Name")
            if name == None:
                variantName = none
            else:
                variantName = name.find("ElementValue").text
            self.variant = variant(measureSet.find("Measure"), variantName, 
                                   measureSet.get("ID"), debug=debug)
                

class clinVarAssertion:
    """Class for representing one submission (i.e. one annotation of a 
    submitted variant"""
    
    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ClinVarAssertion", self.id)
        cvsd = element.find("ClinVarSubmissionID")
        if cvsd == None:
            self.submitter = None
            self.dateSubmitted = None
        else:
            self.submitter = cvsd.get("submitter", default=None)
            self.dateSubmitted = cvsd.get("submitterDate")
        cva = element.find("ClinVarAccession")
        if cva == None:
            self.accession = None
        else:
            self.accession = cva.get("Acc", default=None)
        self.origin = None
        self.method = None
        oi = element.find("ObservedIn")
        if oi != None:
            sample = oi.find("Sample")
            if sample != None:
                self.origin = textIfPresent(sample, "Origin")
            method = oi.find("Method")
            if method != None:
                self.method = textIfPresent(method, "MethodType")
        self.clinicalSignificance = None
        self.reviewStatus = None
        self.dateLastUpdated = None
        cs = element.find("ClinicalSignificance")
        if cs != None:
            self.dateLastUpdated = cs.get("DateLastEvaluated")
            self.clinicalSignificance = textIfPresent(cs, "Description")
            self.reviewStatus = textIfPresent(cs, "ReviewStatus")
                 

class clinVarSet:
    """Container class for a ClinVarSet record, which is a set of submissions
    that were submitted to ClinVar together.  In the ClinVar terminology, 
    each ClinVarSet is one aggregate record ("RCV Accession"), which contains
    one or more submissions ("SCV Accessions").
    """

    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ClinVarSet ID", self.id)
        rcva = element.find("ReferenceClinVarAssertion")
        if isCurrent(rcva):
            self.referenceAssertion = referenceAssertion(rcva, debug=debug)
        self.otherAssertions = dict()
        for item in element.findall("ClinVarAssertion"):
            if isCurrent(item):
                cva = clinVarAssertion(item)
                accession = cva.accession
                self.otherAssertions[accession] = cva

            

        
