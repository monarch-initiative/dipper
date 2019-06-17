#!/usr/bin/env bash

# Generates test set files needed for dipper tests in tests/resources/clinvar/input
# Requires xmlstarlet installed and on your path (and a decent amount of memory, >32gb)
#
# Usage:
# clinvar.sh
# mv RCV*.xml.gz /path/to/tests/resources/clinvar/input

# Target RCVs
rcvs=(
    'RCV000112698'
    'RCV000162061'
    'RCV000175394'
    'RCV000416376'
    'RCV000498447'
    'RCV000763295'
    'RCV000087646'
)

#  ClinVar XML file
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

for rcv in "${rcvs[@]}"
do
    echo $rcv > ./tmp_rcv_list.txt
    ../ClinVarXML_Subset.sh tmp_rcv_list.txt ClinVarFullRelease_00-latest.xml.gz >${rcv}.xml
    gzip ${rcv}.xml

done

rm ./tmp_rcv_list.txt

# Get gene_condition_source_id
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id
grep "BRCA1\|FBN1\|FBN2\|HADHA\|POLR3B\|ASPM\|NPHP4" gene_condition_source_id > gene_condition_test_set.tsv
rm gene_condition_source_id
