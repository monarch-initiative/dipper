#!/usr/bin/env bash
#
# Generates test set files needed for dipper tests in tests/resources/gwascatalog/input
#
# Usage:
# gwas.sh /path/to/gwas-catalog-associations_ontology-annotated.tsv
# mv ./gwascatalog-*.tsv /path/to/tests/resources/gwascatalog/input/
#
# Tested Variants:
#
# rs1329573-?; rs7020413-?; rs3824344-?; rs3758171-?: haplotype block
# kgp8851185-? kgp variant, 100k genomes
# rs1491921-C intergenic variant
#
set -e

# Target Variant IDs
variants=(
    'rs1329573-?; rs7020413-?; rs3824344-?; rs3758171-?'
    'kgp8851185'
    'rs1491921-C'
)

#  Gwas full file
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv

for variant in "${variants[@]}"
do
    filename=$(echo ${variant} | sed 's/-.*//')
    head -1 gwas-catalog-associations_ontology-annotated.tsv > gwascatalog-${filename}.tsv
    grep "${variant}" gwas-catalog-associations_ontology-annotated.tsv >> gwascatalog-${filename}.tsv
done
