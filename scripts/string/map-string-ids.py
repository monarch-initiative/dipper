""" Map StringDB Protein IDs (Ensembl) to Mod Db Genes

Finds identifier mappings between Ensembl protein identifiers
found in StringDB and gene IDs in model organism dbs

"""
import argparse
import pymysql
from pymysql.cursors import Cursor
from pymysql import Connection
from dipper.sources.Ensembl import Ensembl
import yaml
import logging
from typing import List, Union
from intermine.webservice import Service
import requests
from pathlib import Path
import gzip
import pandas as pd


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main():

    """
    Zebrafish:
        1. Map ENSP to ZFIN Ids using Intermine
        2. Map deprecated ENSP IDs to ensembl genes
           by querying the ensembl database then use
           intermine to resolve to gene IDs
    Mouse: Map deprecated ENSP IDs to ensembl genes
           by querying the ensembl database then use
           intermine to resolve to MGI IDs
    Fly: ENSP IDs appear as xrefs on translation IDs
    Worm: Use UniProt Mapping file provided by String
    """

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--config', '-c', required=True, help='JSON configuration file')
    parser.add_argument('--out', '-o', required=False, help='output directory', default="./")
    parser.add_argument('--use_cache', '-cached', action="store_true",
                        required=False, help='use cached files', default=False)
    args = parser.parse_args()

    # Hardcoded dir for raw files
    out_path = Path(args.out)
    raw_dir = out_path / "out"
    raw_dir.mkdir(parents=True, exist_ok=True)

    # Hardcoded unmapped file

    VERSION = 'v10.5'
    STRING_BASE = "http://string-db.org/download/" \
                  "protein.links.detailed.{}".format(VERSION)

    config_file = open(args.config, 'r')
    config = yaml.load(config_file)
    config_file.close()

    out_unmapped_file = out_path / "unmapped_ids.tsv"
    unmapped_file = out_unmapped_file.open("w")

    # Connect to ensembl
    connection = connect_to_database(host=config['database']['host'],
                                     username=config['database']['username'],
                                     port=config['database']['port'])

    cursor = connection.cursor()

    # Process MGI eqs #
    ####################
    taxon = config['taxa_specific']['mouse']['tax_id']

    # IO
    dump_file = raw_dir / '{}.protein.links.detailed.{}.txt.gz' \
        .format(taxon, VERSION)

    mouse_map_file = out_path / config['taxa_specific']['mouse']['output_file']
    mouse_file = mouse_map_file.open('w')

    path = '{}/{}.protein.links.detailed.{}.txt.gz' \
        .format(STRING_BASE, taxon, VERSION)
    if not args.use_cache:
        download_file(path, dump_file)

    ensembl = Ensembl("rdf_graph", True)
    p2gene_map = ensembl.fetch_protein_gene_map(taxon)

    fh = gzip.open(str(dump_file), 'rb')
    df = pd.read_csv(fh, sep='\s+')
    fh.close()
    proteins = pd.unique(df[['protein1', 'protein2']].values.ravel())
    logger.info("Processing {} proteins".format(len(proteins)))
    for protein in proteins:
        prot = protein.replace('{}.'.format(str(taxon)), '')
        try:
            ens_gene = p2gene_map[prot]
            ens_curie = "ENSEMBL:{}".format(ens_gene)
            mouse_file.write("{}\t{}\n".format(prot, ens_curie))
            continue
        except KeyError:
            pass

        ens_gene = get_deprecated_protein_gene_rel(
            cursor, prot, config['taxa_specific']['mouse']['ensembl'],
            config)
        intermine_resp = query_mousemine(
            config['taxa_specific']['mouse']['intermine'], ens_gene)
        if intermine_resp.is_successful:
            mouse_file.write("{}\t{}\n".format(prot, intermine_resp.gene_id))
        else:
            unmapped_file.write("{}\t{}\t{}\n".format(prot, ens_gene, taxon))

    mouse_file.close()

    # Process Fly eqs #
    ####################
    taxon = config['taxa_specific']['fly']['tax_id']

    # IO
    dump_file = raw_dir / '{}.protein.links.detailed.{}.txt.gz' \
        .format(taxon, VERSION)

    fly_map_file = out_path / config['taxa_specific']['fly']['output_file']
    fly_file = fly_map_file.open('w')

    path = '{}/{}.protein.links.detailed.{}.txt.gz' \
        .format(STRING_BASE, taxon, VERSION)
    if not args.use_cache:
        download_file(path, dump_file)

    ensembl = Ensembl("rdf_graph", True)
    p2gene_map = ensembl.fetch_protein_gene_map(taxon)

    fh = gzip.open(str(dump_file), 'rb')
    df = pd.read_csv(fh, sep='\s+')
    fh.close()
    proteins = pd.unique(df[['protein1', 'protein2']].values.ravel())
    logger.info("Processing {} proteins".format(len(proteins)))
    for protein in proteins:
        prot = protein.replace('{}.'.format(str(taxon)), '')
        try:
            ens_gene = p2gene_map[prot]
            ens_curie = "ENSEMBL:{}".format(ens_gene)
            fly_file.write("{}\t{}\n".format(prot, ens_curie))
            continue
        except KeyError:
            pass

        ens_gene = get_xref_protein_gene_rel(
            cursor, prot, config['taxa_specific']['fly']['ensembl'],
            config, taxon)

        if ens_gene is not None:
            fly_file.write("{}\t{}\n".format(prot, "ENSEMBL:{}".format(ens_gene)))
        else:
            unmapped_file.write("{}\t{}\t{}\n".format(prot, '', taxon))

    fly_file.close()

    # Process Worm eqs #
    ####################
    taxon = config['taxa_specific']['worm']['tax_id']

    # IO
    dump_file = raw_dir / '{}.protein.links.detailed.{}.txt.gz' \
        .format(taxon, VERSION)

    uniprot_file = raw_dir / config['taxa_specific']['worm']['uniprot_file']

    worm_map_file = out_path / config['taxa_specific']['worm']['output_file']
    worm_file = worm_map_file.open('w')

    path = '{}/{}.protein.links.detailed.{}.txt.gz' \
        .format(STRING_BASE, taxon, VERSION)
    if not args.use_cache:
        download_file(path, dump_file)
        download_file(config['taxa_specific']['worm']['uniprot_mappings'],
                      uniprot_file)

    ensembl = Ensembl("rdf_graph", True)
    p2gene_map = ensembl.fetch_protein_gene_map(taxon)
    uni2gene_map = ensembl.fetch_uniprot_gene_map(taxon)

    fh = gzip.open(str(uniprot_file), 'rb')
    df = pd.read_csv(fh, sep='\s+')
    fh.close()
    string_uniprot_map = {}
    for index, row in df.iterrows():
        uniprot_ac = row['uniprot_ac|uniprot_id'].split('|')[0]
        string_uniprot_map[row['string_id']] = uniprot_ac

    fh = gzip.open(str(dump_file), 'rb')
    df = pd.read_csv(fh, sep='\s+')
    fh.close()
    proteins = pd.unique(df[['protein1', 'protein2']].values.ravel())
    logger.info("Processing {} proteins".format(len(proteins)))
    for protein in proteins:
        prot = protein.replace('{}.'.format(str(taxon)), '')
        try:
            ens_gene = p2gene_map[prot]
            ens_curie = "ENSEMBL:{}".format(ens_gene)
            worm_file.write("{}\t{}\n".format(prot, ens_curie))
            continue
        except KeyError:
            pass

        try:
            uniprot_ac = string_uniprot_map[prot]
            ens_gene = uni2gene_map[uniprot_ac]
            ens_curie = "ENSEMBL:{}".format(ens_gene)
            worm_file.write("{}\t{}\n".format(prot, ens_curie))
            continue
        except KeyError:
            pass

        unmapped_file.write("{}\t{}\t{}\n".format(prot, '', taxon))

    worm_file.close()

    # Process ZFIN eqs #
    ####################
    taxon = config['taxa_specific']['zebrafish']['tax_id']

    # IO
    dump_file = raw_dir / '{}.protein.links.detailed.{}.txt.gz' \
        .format(taxon, VERSION)

    zfin_map_file = out_path / config['taxa_specific']['zebrafish']['output_file']
    zfin_file = zfin_map_file.open('w')

    path = '{}/{}.protein.links.detailed.{}.txt.gz' \
        .format(STRING_BASE, taxon, VERSION)
    if not args.use_cache:
        download_file(path, dump_file)

    ensembl = Ensembl("rdf_graph", True)
    p2gene_map = ensembl.fetch_protein_gene_map(taxon)

    # in 3.6 gzip accepts Paths
    fh = gzip.open(str(dump_file), 'rb')
    df = pd.read_csv(fh, sep='\s+')
    fh.close()
    proteins = pd.unique(df[['protein1', 'protein2']].values.ravel())
    logger.info("Processing {} proteins".format(len(proteins)))
    for protein in proteins:
        prot = protein.replace('{}.'.format(str(taxon)), '')
        try:
            ens_gene = p2gene_map[prot]
            ens_curie = "ENSEMBL:{}".format(ens_gene)
            zfin_file.write("{}\t{}\n".format(prot, ens_curie))
            continue
        except KeyError:
            pass

        intermine_resp = query_fishmine(
            config['taxa_specific']['zebrafish']['intermine'], prot)
        if intermine_resp.is_successful:
            zfin_file.write("{}\t{}\n".format(prot, intermine_resp.gene_id))
            continue

        ens_gene = get_deprecated_protein_gene_rel(
            cursor, prot, config['taxa_specific']['zebrafish']['ensembl'],
            config)
        intermine_resp = query_fishmine(
            config['taxa_specific']['zebrafish']['intermine'], ens_gene)
        if intermine_resp.is_successful:
            zfin_file.write("{}\t{}\n".format(prot, intermine_resp.gene_id))
            continue

        intermine_resp = query_fishmine(
            config['taxa_specific']['zebrafish']['intermine'],
            ens_gene, "Pseudogene")
        if intermine_resp.is_successful:
            zfin_file.write("{}\t{}\n".format(prot, intermine_resp.gene_id))
        else:
            unmapped_file.write("{}\t{}\t{}\n".format(prot, ens_gene, taxon))

    zfin_file.close()

    unmapped_file.close()
    connection.close()

    logger.info("ID Map Finished")


class IntermineResult(object):
    def __init__(self, is_successful: bool, gene_id: str=None,  msg: str=None):
        self.is_successful = is_successful
        self.gene_id = gene_id
        self.msg = msg


def query_fishmine(intermine_url: str, protein_id: str, query: str="Gene") -> IntermineResult:
    service = Service(intermine_url)
    query = service.new_query(query)
    query.add_view("primaryIdentifier")
    query.add_constraint("primaryIdentifier", "CONTAINS", "ZDB*", code="A")
    query.add_constraint("crossReferences.identifier", "=", "{}".format(protein_id), code="B")
    result_list = ["ZFIN:{}".format(val['primaryIdentifier']) for val in query.rows()]
    return intermine_response_factory(result_list, protein_id)


def query_mousemine(intermine_url: str, gene_id: str) -> IntermineResult:
    """
    :param intermine_url: intermine server, eg
                          http://www.mousemine.org/mousemine/service
    :param gene_id: gene ID, eg ENSMUSG00000063180
    :return: Intermine_Result object
    """
    service = Service(intermine_url)
    query = service.new_query("SequenceFeature")
    query.add_view("primaryIdentifier")
    query.add_constraint("SequenceFeature", "LOOKUP", "{}".format(gene_id), code="A")
    query.add_constraint("organism.shortName", "=", "M. musculus", code="B")
    result_list = ["{}".format(val['primaryIdentifier']) for val in query.rows()]
    return intermine_response_factory(result_list, gene_id)


def intermine_response_factory(result_list: List[str], input_id: str) -> IntermineResult:
    if len(result_list) == 1:
        gene = result_list[0]
        is_successful = True
        intermine_res = IntermineResult(is_successful, gene)
    elif len(result_list) > 1:
        msg = "Intermine: Ambiguous gene mapping"
        logger.warn(msg + " for {}: {}".format(input_id, result_list))
        is_successful = False
        intermine_res = IntermineResult(is_successful, msg=msg)
    else:
        msg = "Intermine: No mapping found"
        logger.warn(msg + " for {}".format(input_id))
        is_successful = False
        intermine_res = IntermineResult(is_successful, msg=msg)

    return intermine_res


def select_database(cursor: Cursor, database: str) -> None:
    query = "USE {};".format(database)
    cursor.execute(query)


def connect_to_database(host: str,
                        username: str=None,
                        password: str=None,
                        database: str=None,
                        port: Union[str, int]=None) -> Connection:
    logger.info("Connecting to database %s on %s", database, host)

    connection = pymysql.connect(host=host, user=username,
                                 passwd=password, db=database, port=port)
    logger.info("Connected to {}:{}".format(host, database))
    return connection


def download_file(remote_path: str, local_path:Path) -> None:
    req = requests.get(remote_path)
    logger.info("Downloading remote file: "
                "{} to {}".format(remote_path, local_path))
    if req.status_code == 200:
        fh = local_path.open('wb')
        for chunk in req:
            fh.write(chunk)
    else:
        raise IOError("Unable to fetch file {}".format(remote_path))


def get_xref_protein_gene_rel(cursor: Cursor, protein: str,
                              database: str, config: dict, taxon: str) -> str:
    select_database(cursor, database)
    find_gene_query = \
        config['queries']['protein2xref'].format("{}.{}".format(taxon, protein))
    cursor.execute(find_gene_query)
    query_res = cursor.fetchall()

    if len(query_res) > 1:
        logger.warn("Ambiguous protein to gene "
                    "mapping for {} {}".format(protein, query_res))
        gene = None
    elif len(query_res) == 0:
        logger.warn("No mapping for {}".format(protein))
        gene = None
    else:
        gene = query_res[0][1]
    return gene


def get_deprecated_protein_gene_rel(cursor: Cursor, protein: str,
                                    database: str, config: dict) -> str:
    select_database(cursor, database)
    find_db_query = config['queries']['id_mapping'].format(protein)
    cursor.execute(find_db_query)
    old_db = cursor.fetchall()
    try:
        select_database(cursor, old_db[0][0])
    except pymysql.err.InternalError:
        return ''
    except IndexError:
        return ''

    find_gene_query = config['queries']['protein2gene'].format(protein)
    try:
        cursor.execute(find_gene_query)
    except pymysql.err.ProgrammingError:
        find_gene_query = config['queries']['p2g_two'].format(protein)
        cursor.execute(find_gene_query)
    query_res = cursor.fetchall()
    if len(query_res) > 1:
        raise ValueError("Ambiguous protein to gene "
                         "mapping for {} {}".format(protein, query_res))
    elif len(query_res) == 0:
        raise ValueError("No mapping for protein {}".format(protein))
    return query_res[0][0]

if __name__ == "__main__":
    main()
