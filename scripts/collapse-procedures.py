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

    proc_set = set()
    tmp_proc_map = {}
    proc_map = {}

    for mp in procedure_list:
        for code in mp:
            tmp_proc_map[code] = mp[code]

    ftp = 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/ALL_genotype_phenotype.csv.gz'
    with urllib.request.urlopen(ftp) as response:
        buf = io.BytesIO(response.read())
        f = gzip.GzipFile(fileobj=buf)
        data = f.read()

    reader = csv.reader(data.split('\n'), delimiter=',', quotechar='\"')
    for row in reader:
        if not row or row[0] == 'marker_accession_id':
            continue
        else:
            proc_set.add(row[17])

    for procedure in proc_set:
        proc_map[procedure] = tmp_proc_map[procedure]

    json.dump(proc_map, output_fh)

    output_fh.close()


if __name__ == "__main__":
    main()
