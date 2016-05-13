#!/usr/bin/env python3
import argparse
import json
import urllib.request
import io
import csv
import gzip


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Location of output file')

    args = parser.parse_args()
    output_fh = open(args.output, 'w')
    procedure_list = json.load(open(args.input, 'r'))

    param_set = set()
    tmp_param_map = {}
    param_map = {}

    for mp in procedure_list:
        for code in mp:
            tmp_param_map[code] = mp[code]

    ftp = 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/ALL_genotype_phenotype.csv.gz'
    with urllib.request.urlopen(ftp) as response:
        buf = io.BytesIO(response.read())
        f = gzip.GzipFile(fileobj=buf)
        data = f.read().decode()

    reader = csv.reader(data.split('\n'), delimiter=',', quotechar='\"')
    for row in reader:
        if not row or row[0] == 'marker_accession_id':
            continue
        else:
            param_set.add(row[17])

    for param in param_set:
        param_map[param] = tmp_param_map[param]

    json.dump(param_map, output_fh)

    output_fh.close()


if __name__ == "__main__":
    main()
