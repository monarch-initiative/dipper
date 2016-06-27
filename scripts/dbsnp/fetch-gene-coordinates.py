#!/usr/bin/env python3
import argparse
import requests
import xml.etree.ElementTree as ET
import csv
import re

"""
Input example, tab delimited file with gene symbols and NCBI curie formatted ids
see for example, the output of scripts/fetch-gene-ids:

AANAT   NCBIGene:15
ABCA1   NCBIGene:19
ABCA2   NCBIGene:20

Output
NCBIGene:15     74449432        74466198        plus    GRCh37.p13
NCBIGene:19     107543282       107690526       minus   GRCh37.p13
NCBIGene:20     139901685       139923373       minus   GRCh37.p13

Note this has only been tested with GRCh37.p13
"""


def main():

    EFETCH_BASE = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')
    parser.add_argument('--build', '-b', type=str, default="GRCh37.p13",
                        help='Location of output file')

    args = parser.parse_args()

    gene_list = []

    # Open file
    with open(args.input, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        for row in reader:
            gene = row[1].rstrip('\n')
            gene = re.sub(r'NCBIGene:', '', gene)
            if gene != '':
                gene_list.append(gene)

    output_file = open(args.output, 'w')

    # Split into chunks of 100, courtesy http://stackoverflow.com/a/1751478
    n = max(1, 20)
    chunked_list = [gene_list[i:i + 20] for i in range(0, len(gene_list), 20)]

    for gene_chunk in chunked_list:
        print("Fetching chunk of ids: {0}".format(gene_chunk))
        params = {
            'id': gene_chunk,
            'rettype': 'xml',
            'db': 'gene'
        }
        eutils_request = requests.get(EFETCH_BASE, params)
        root = ET.fromstring(eutils_request.content)
        for entrezgene in root.findall("Entrezgene"):
            gene_id = entrezgene.find(".//Entrezgene_track-info/Gene-track/Gene-track_geneid").text
            build_elt = \
                entrezgene.find(".//Entrezgene_comments/Gene-commentary/Gene-commentary_comment/"
                                "Gene-commentary/Gene-commentary_comment/"
                                "Gene-commentary[Gene-commentary_heading='" + args.build + "']")

            if build_elt is not None:

                try:
                    interval_start = \
                        build_elt.find(".//Gene-commentary_comment/Gene-commentary/"
                                       "Gene-commentary_seqs/Seq-loc/Seq-loc_int/"
                                       "Seq-interval/Seq-interval_from").text
                except AttributeError:
                    interval_start = ''

                try:
                    interval_end = \
                        build_elt.find(".//Gene-commentary_comment/Gene-commentary/"
                                       "Gene-commentary_seqs/Seq-loc/Seq-loc_int/"
                                       "Seq-interval/Seq-interval_to").text
                except AttributeError:
                    interval_end = ''

                try:
                    strand = \
                        build_elt.find(".//Gene-commentary_comment/Gene-commentary/"
                                       "Gene-commentary_seqs/Seq-loc/Seq-loc_int/"
                                       "Seq-interval/Seq-interval_strand/Na-strand").attrib['value']
                except AttributeError:
                    strand = ''

            else:
                interval_start = ''
                interval_end = ''
                strand = ''

            output_file.write("NCBIGene:{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                gene_id, interval_start, interval_end, strand, args.build))

            # sanity check
            if gene_id not in gene_chunk:
                print("error, gene id not found in input")
                exit(1)

    output_file.close()

if __name__ == "__main__":
    main()
