#!/usr/bin/env python3

"""
Script to generate entrez gene IDs from HGNC gene symbols in the format
NCBIGene:1234

Queries mygene and monarch and joins the results, noting conflicts
likely from dupe/ambiguous ids in both dbs
"""

import mygene
import argparse
import requests
import time

MONARCH_URL = 'https://beta.monarchinitiative.org/search/'
requests.adapters.DEFAULT_RETRIES = 10

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', type=str, help='Path to input file')
parser.add_argument('--output', '-o', type=str, help='Path to output file')

args = parser.parse_args()
mg = mygene.MyGeneInfo()
symbol_list = []
output_file = open(args.output, 'w')

# Get symbols from mapping file

with open(args.input, 'rt') as tsvfile:
    for row in tsvfile:
        gene_symbol = row.rstrip('\n')
        symbol_list.append(gene_symbol)

print("Querying mygene")
mygene_results = mg.querymany(symbol_list, scopes='symbol', fields='entrezgene', species='human', returnall=True)

mygene_output = mygene_results['out']
mygene_dict = {}
monarch_dict = {}

for result in mygene_output:
    if 'entrezgene' in result:
        mygene_dict[result['query']] = "NCBIGene:{0}".format(result['entrezgene'])

# Now try with monarch
for gene in symbol_list:
    url = MONARCH_URL + gene + '.json'
    try:
        monarch_request = requests.get(url)
        results = monarch_request.json()
        for res in results['results']:
            if res['taxon'] == 'Human' \
                    and 'gene' in res['categories'] \
                    and (gene in res['labels'] or gene in res['synonyms']):
                monarch_dict[gene] = res['id']
                break
    except requests.ConnectionError:
        time.sleep(20)
        print("error fetching {0}".format(url))

for gene in symbol_list:
    if gene in mygene_dict and gene not in monarch_dict:
        line = "{0}\t{1}".format(gene, mygene_dict[gene])
    elif gene not in mygene_dict and gene in monarch_dict:
        line = "{0}\t{1}".format(gene, monarch_dict[gene])
    elif gene in mygene_dict and gene in monarch_dict:
        # Check if they agree
        if monarch_dict[gene] == mygene_dict[gene]:
            line = "{0}\t{1}".format(gene, monarch_dict[gene])
        else:
            print("monarch and mygene disagree on {0},"
                  " mygene says {1}, monarch says {2}".format(
                      gene, mygene_dict[gene], monarch_dict[gene]))
            line = "{0}\t{1}".format(gene, monarch_dict[gene])
    else:
        line = "{0}\t".format(gene)

    output_file.write(line+"\n")

output_file.close()