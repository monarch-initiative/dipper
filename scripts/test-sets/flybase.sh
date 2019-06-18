#!/usr/bin/env bash

# Generates test set files needed for dipper tests in tests/resources/flybase/input
# Requires dipper for fetching,
# cd /path/to/local/dipper
# python install setup.py
#
# Usage:
# flybase.sh /path/to/venv
# mv ./raw/flybase/FBal* ../../tests/resources/flybase/input/
#
# Fetches all flybase files and creates a raw/flybase directory relative
# to where the script is run, then greps out the lines we need for testing
# specific alleles
#
# Alleles:
# http://flybase.org/reports/FBal0195705, standard allele, fly gene
# http://flybase.org/reports/FBal0256668:
#      humana allele, human transgene, 1 phenotype, 1 disease model
#
# Genes:
# FBgn0033159
# FBgn0250787

# Potential alleles of interest:
# http://flybase.org/reports/FBal0190789,FBgn0084473 transgene, only phenotype manifests in

test_set_1=(
    'FBal0195705'
    'FBgn0033159'
)

test_set_2=(
    'FBal0256668'
    'FBgn0250787'
)

counter=1
max_test_sets=2

PTH=`pwd`
# Default looks back two directories for a venv
VENV=${1:-"${PTH}/../../venv"}

source ${VENV}/bin/activate

dipper-etl.py --sources flybase --fetch_only

cd raw/flybase

# Strip all fields from fbrf_pmid_pmcid_doi_fb except 1,2
zcat fbrf_pmid_pmcid_doi_fb.tsv.gz | head -4 > fbrf_pmid_pmcid_doi_fb.tsv
zcat fbrf_pmid_pmcid_doi_fb.tsv.gz | tail -n+4 | cut -f1,2 >> fbrf_pmid_pmcid_doi_fb.tsv
rm fbrf_pmid_pmcid_doi_fb.tsv.gz
gzip fbrf_pmid_pmcid_doi_fb.tsv

while [[ ${counter} -le ${max_test_sets} ]]
do
    test_set="test_set_${counter}[@]"
    test_set_var=(${!test_set})
    allele_id=${test_set_var[0]}
    gene_id=${test_set_var[1]}

    mkdir ${allele_id} && cd ${allele_id}
    cp ../fbrf_pmid_pmcid_doi_fb.tsv.gz ./
    cp ../species.ab.gz ./

    # Allele to gene file
    zcat ../fbal_to_fbgn_fb.tsv.gz | head -2 > fbal_to_fbgn_fb.tsv
    zgrep ${allele_id} ../fbal_to_fbgn_fb.tsv.gz >> fbal_to_fbgn_fb.tsv
    gzip fbal_to_fbgn_fb.tsv

    # Gene xref file
    head -1 ../gene_xref.tsv > gene_xref.tsv
    grep ${gene_id} ../gene_xref.tsv >> gene_xref.tsv

    # Allele phenotype file
    head -1 ../allele_phenotype.tsv > allele_phenotype.tsv
    grep "^${allele_id}" ../allele_phenotype.tsv >> allele_phenotype.tsv

    # Model of disease
    zcat ../disease_model_annotations.tsv.gz | head -5 > disease_model_annotations.tsv
    zgrep -P "\t${allele_id}\t" ../disease_model_annotations.tsv.gz >> disease_model_annotations.tsv
    gzip disease_model_annotations.tsv

    cd ..
    counter=$[$counter+1]
done
