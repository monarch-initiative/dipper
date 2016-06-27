#!/usr/bin/env python3
import argparse
import requests
import xml.etree.ElementTree as ET


"""
Input example, whitespace separated file of chromosome, location, rsid
1 100343254 141043166
1 103355059 55821405
1 103471457 565657269
1 103471457 71752747
1 103471457 763506645


Output: tab delimited file with additional columns variant type, and observed alleles
1       100343254       141043166       snp     A/C/G
1       103355059       55821405        snp     G/T
1       103471457       565657269       in-del  -/CATCAT
1       103471457       71752747        in-del  -/CAT
1       103471457       763506645       in-del  -/CATCATCAT
"""


def main():

    EFETCH_BASE = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    NS = '{http://www.ncbi.nlm.nih.gov/SNP/docsum}'

    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')

    args = parser.parse_args()

    # Open file
    with open(args.input) as file:
        rs_ids = file.read().splitlines()

    output_file = open(args.output, 'w')

    # Split into chunks of 100, courtesy http://stackoverflow.com/a/1751478
    n = max(1, 100)
    chunked_list = [rs_ids[i:i + 100] for i in range(0, len(rs_ids), 100)]

    for snp_chunk in chunked_list:
        rs_list = [rs.split(' ')[2] for rs in snp_chunk]

        params = {
            'id': rs_list,
            'report': 'xml',
            'db': 'snp'
        }
        eutils_request = requests.get(EFETCH_BASE, params)
        root = ET.fromstring(eutils_request.content)
        for snp in snp_chunk:
            info = snp.split(' ')
            chromosome = info[0]
            coordinate = info[1]
            snp_id = info[2]

            rs = root.find(".//"+NS+"Rs[@rsId='"+snp_id+"']")
            try:
                snp_class = rs.attrib['snpClass']
                sequence = rs.find(NS+'Sequence')
                allele_diff = sequence.find(NS+'Observed').text
            except AttributeError:
                snp_class = ''
                allele_diff = ''
            output_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                chromosome, coordinate, snp_id, snp_class, allele_diff))

    output_file.close()

if __name__ == "__main__":
    main()
