#! /bin/bash

# extract ClinVarSets containing specific SCV
# not at all efficent, but does not need to be.
# give a list of SCV identifiers and a dataset to find them in
#
# Usage:
# ClinVarXML_Subset_SCV.sh	 [<scvlist> <cvxml>] > ClinVarTestSet.xml
# e.g.
# ./scripts/ClinVarXML_Subset.sh | gzip > raw/clinvarxml_alpha/ClinVarTestSet.xml.gz

RPTH='raw/clinvarxml_alpha'
# Defaults if not given
TEST=${1:-"${RPTH}/CV_test_RCV.txt"}
CXML=${2:-"${RPTH}/ClinVarFullRelease_00-latest.xml.gz"}

CVSET='ReleaseSet/ClinVarSet'
SCV='ClinVarAssertion/ClinVarAccession/@Acc'

### header
if [ 'gz' == "${CXML##*.}" ] ; then
    zcat "${CXML}" | head -3
else
    head -3 "${CXML}"
fi

### body
for scv in `< "${TEST}"`; do
    xmlstarlet sel -t -c "${CVSET}[${SCV}=\"${scv}\"]" -n "${CXML}"
done

### footer
if [ 'gz' == "${CXML##*.}" ] ; then
    zcat "${CXML}" | tail -2
else
    tail -2  "${CXML}"
fi
