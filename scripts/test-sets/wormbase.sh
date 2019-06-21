#!/usr/bin/env bash
#
# Original test ids from py source:
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

set -e

test_genes=(
    'WBGene00001414'
    'WBGene00004967'
    'WBGene00003916'
)


PTH=`pwd`
RAWDIR=${1:-"${PTH}/raw/wormbase"}
# Default looks back two directories
VENV=${2:-"${PTH}/../../venv"}

if [[ -z "$1" ]]; then
    source ${VENV}/bin/activate
    dipper-etl.py --sources wormbase --fetch_only
fi


for gene in "${test_genes[@]}"
do
  echo "TODO"
done
