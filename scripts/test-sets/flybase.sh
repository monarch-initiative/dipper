#!/usr/bin/env bash
#
# Generates test set files needed for dipper tests in tests/resources/flybase/input
# If no raw dir is supplied, requires dipper on venv for fetching,
# source venv/bin/activate
# cd /path/to/local/dipper
# python install setup.py
#
# Usage:
# flybase.sh /path/to/raw/ /path/to/venv
# mv ./FBal* /path/to/tests/resources/flybase/input/
#
# Fetches all flybase files and creates a raw/flybase directory relative
# to where the script is run, then greps out the lines we need for testing
# specific alleles
#
# Alleles:
# http://flybase.org/reports/FBal0195705, standard allele, fly gene
# http://flybase.org/reports/FBal0263199, autism model
# http://flybase.org/reports/FBal0256668:
#      human allele, human transgene, 1 phenotype, 1 disease model
#
set -e

test_set_1=(
    'FBal0195705'
    'FBgn0033159'
)

# Transgene, should be filtered out of ingest
test_set_2=(
    'FBal0256668'
    'FBgn0250787'
)

test_set_3=(
    'FBal0263199'
    'FBgn0028734'
)

counter=1
max_test_sets=3

PTH=`pwd`
RAWDIR=${1:-"${PTH}/raw/flybase"}
# Default looks back two directories
VENV=${2:-"${PTH}/../../venv"}

if [[ -z "$1" ]]; then
    source ${VENV}/bin/activate
    dipper-etl.py --sources flybase --fetch_only
fi

# Strip all fields from fbrf_pmid_pmcid_doi_fb except 1,2
zcat $RAWDIR/fbrf_pmid_pmcid_doi_fb.tsv.gz | head -4 > fbrf_pmid_pmcid_doi_fb.tmp
zcat $RAWDIR/fbrf_pmid_pmcid_doi_fb.tsv.gz | tail -n+4 | cut -f1,2 >> fbrf_pmid_pmcid_doi_fb.tmp
gzip fbrf_pmid_pmcid_doi_fb.tmp

while [[ ${counter} -le ${max_test_sets} ]]
do
    test_set="test_set_${counter}[@]"
    test_set_var=(${!test_set})
    allele_id=${test_set_var[0]}
    gene_id=${test_set_var[1]}

    outdir=$allele_id

    mkdir ${outdir}
    cp fbrf_pmid_pmcid_doi_fb.tmp.gz ${outdir}/fbrf_pmid_pmcid_doi_fb.tsv.gz
    cp $RAWDIR/species.ab.gz ${outdir}

    # Allele to gene file
    outfile=${outdir}/fbal_to_fbgn_fb.tsv
    zcat $RAWDIR/fbal_to_fbgn_fb.tsv.gz | head -2 > ${outfile}
    zgrep ${allele_id} $RAWDIR/fbal_to_fbgn_fb.tsv.gz >> ${outfile}
    gzip ${outfile}

    # Gene xref file
    outfile=${outdir}/gene_xref.tsv
    head -1 $RAWDIR/gene_xref.tsv > ${outfile}
    grep ${gene_id} $RAWDIR/gene_xref.tsv >> ${outfile}

    # Allele phenotype file
    outfile=${outdir}/allele_phenotype.tsv
    head -1 $RAWDIR/allele_phenotype.tsv > ${outfile}
    grep "^${allele_id}" $RAWDIR/allele_phenotype.tsv >> ${outfile}

    # Model of disease
    outfile=${outdir}/disease_model_annotations.tsv
    zcat $RAWDIR/disease_model_annotations.tsv.gz | head -5 > ${outfile}
    zgrep -P "\t${allele_id}\t" $RAWDIR/disease_model_annotations.tsv.gz >> ${outfile}
    gzip ${outfile}

    counter=$[$counter+1]
done

rm fbrf_pmid_pmcid_doi_fb.tmp.gz
