#! /bin/bash

# extract ClinVarSets containing specific RCV
# not at all efficent, but does not need to be.
# give a list of RCV identifiers and a dataset to find them in
#
# Usage:
# ClinVarXML_Subset.sh	 [<rcvlist> <cvxml>] > ClinVarTestSet.xml
# e.g.
# ./scripts/ClinVarXML_Subset.sh | gzip > raw/clinvarxml_alpha/ClinVarTestSet.xml.gz

RPTH='raw/clinvarxml_alpha'
# Defaults if not given
TEST=${1:-"${RPTH}/CV_test_RCV.txt"}
CXML=${2:-"${RPTH}/ClinVarFullRelease_00-latest.xml.gz"}

CVSET='ReleaseSet/ClinVarSet'
RCV='ReferenceClinVarAssertion/ClinVarAccession/@Acc'

### header
if [ 'gz' == "${CXML##*.}" ] ; then
    zcat "${CXML}" | head -3
else
    head -3 "${CXML}"
fi

### body
for rcv in `< "${TEST}"`; do
    xmlstarlet sel -t -c "${CVSET}[${RCV}=\"${rcv}\"]" -n "${CXML}"
done

### footer
if [ 'gz' == "${CXML##*.}" ] ; then
    zcat "${CXML}" | tail -2
else
    tail -2  "${CXML}"
fi
