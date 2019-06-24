#!/usr/bin/env bash
#
# Usage:
# wormbase.sh /path/to/raw/ /path/to/venv
# mv ./WB* /path/to/tests/resources/flybase/input/
#
# Original test ids from WormBase.py:
#    test_ids = {
#        'gene': [
#            'WBGene00001414', 'WBGene00004967', 'WBGene00003916', 'WBGene00004397',
#            'WBGene00001531'],
#        'allele': [
#            'WBVar00087800', 'WBVar00087742', 'WBVar00144481', 'WBVar00248869',
#            'WBVar00250630'],
#        'strain': ['BA794', 'RK1', 'HE1006'],
#        'pub': []
#    }
#
#
set -e

test_genes=(
    'WBGene00001414'
    'WBGene00004967'
    'WBGene00003916'
)

PTH=`pwd`

# Default looks back two directories
RAWDIR=${1:-"${PTH}/raw/wormbase"}
VENV=${2:-"${PTH}/../../venv"}

if [[ -z "$1" ]]; then
    source ${VENV}/bin/activate
    dipper-etl.py --sources wormbase --fetch_only
fi


for gene in "${test_genes[@]}"
do
  mkdir ${gene}

  for file in ${RAWDIR}/*
  do
      filename=$(basename $file)
      zgrep $gene ${file} > ${gene}/$filename || echo -n > ${gene}/$filename
      file_type=$(file $file | grep 'gzip compressed data' || echo '')
      if [[ ! -z "$file_type" ]]; then
          new_name=${filename::-3}
          mv ${gene}/$filename ${gene}/$new_name
          gzip ${gene}/${new_name}
      fi
  done

  rm ${gene}/CHECKSUMS

done
