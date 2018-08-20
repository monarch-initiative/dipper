import logging
import urllib
import csv
import http
import xml.etree.ElementTree as etree

from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.Genotype import Genotype
from dipper import config


logger = logging.getLogger(__name__)


class Ensembl(Source):
    """
    This is the processing module for Ensembl.

    It only includes methods to acquire the equivalences between NCBIGene and
    ENSG ids using ENSEMBL's Biomart services.

    """

    files = {
        '9606': {'file': 'ensembl_9606.txt'},  # human
        '7955': {'file': 'ensembl_7955.txt'},  # fish
        '10090': {'file': 'ensembl_10090.txt'},  # mouse
        '28377': {'file': 'ensembl_28377.txt'},  #
        '3702': {'file': 'ensembl_3702.txt'},  #
        '9913': {'file': 'ensembl_9913.txt'},  #
        '6239': {'file': 'ensembl_6239.txt'},  #
        '9615': {'file': 'ensembl_9615.txt'},  #
        '9031': {'file': 'ensembl_9031.txt'},  #
        '44689': {'file': 'ensembl_44689.txt'},  #
        '7227': {'file': 'ensembl_7227.txt'},  #
        '9796': {'file': 'ensembl_9796.txt'},  #
        '9544': {'file': 'ensembl_9544.txt'},  #
        '13616': {'file': 'ensembl_13616.txt'},  #
        '9258': {'file': 'ensembl_9258.txt'},  #
        '9823': {'file': 'ensembl_9823.txt'},  #
        '10116': {'file': 'ensembl_10116.txt'},  #
        '4896': {'file': 'ensembl_4896.txt'},  #
        '31033': {'file': 'ensembl_31033.txt'},  #
        '8364': {'file': 'ensembl_8364.txt'},  #
        '4932': {'file': 'ensembl_4932.txt'},  #
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None, gene_ids=None):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'ensembl',
            ingest_title='ENSEMBL',
            ingest_url='http://uswest.ensembl.org'
            # license_url=None,
            # data_rights=None,
            # file_handle=None
        )

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]

        self.gene_ids = []
        if 'test_ids' not in config.get_config() \
                or 'gene' not in config.get_config()['test_ids']:
            logger.warning("not configured with gene test ids.")
        else:
            self.gene_ids = config.get_config()['test_ids']['gene']

        logger.setLevel(logging.INFO)

        return

    def fetch(self, is_dl_forced=False):

        for t in self.tax_ids:
            logger.info("Fetching genes for %s", str(t))
            loc_file = '/'.join((self.rawdir, 'ensembl_'+str(t)+'.txt'))
            # todo move into util?
            params = urllib.parse.urlencode(
                {'query': self._build_biomart_gene_query(str(t))})
            conn = http.client.HTTPConnection('uswest.ensembl.org')
            conn.request("GET", '/biomart/martservice?'+params)
            resp = conn.getresponse()
            with open(loc_file, 'wb') as f:
                f.write(resp.read())

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        if self.testOnly:
            self.testMode = True

        logger.info("Parsing files...")

        for t in self.tax_ids:
            self._process_genes(str(t), limit)

        logger.info("Done parsing files.")

        return

    def fetch_protein_list(self, taxon_id):
        """
        Fetch a list of proteins for a species in biomart
        :param taxid:
        :return: list
        """
        protein_list = list()
        params = urllib.parse.urlencode(
            {'query': self._build_biomart_gene_query(str(taxon_id))})
        conn = http.client.HTTPConnection('www.ensembl.org')
        conn.request("GET", '/biomart/martservice?' + params)
        response = conn.getresponse()
        for line in response:
            line = line.decode('utf-8').rstrip()
            row = line.split('\t')
            if len(row) < 6:
                logger.warning("Data error for query on %d", taxon_id)
                continue
            (ensembl_gene_id, external_gene_name,
             description, gene_biotype, entrezgene,
             peptide_id) = row[0:6]
            protein_list.append(peptide_id)
        conn.close()
        return protein_list

    def fetch_protein_gene_map(self, taxon_id):
        """
        Fetch a list of proteins for a species in biomart
        :param taxid:
        :return: dict
        """
        protein_dict = dict()
        params = urllib.parse.urlencode(
            {'query': self._build_biomart_gene_query(str(taxon_id))})
        conn = http.client.HTTPConnection('www.ensembl.org')
        conn.request("GET", '/biomart/martservice?' + params)
        response = conn.getresponse()
        for line in response:
            line = line.decode('utf-8').rstrip()
            row = line.split('\t')
            if len(row) < 6:
                logger.warning("Data error for query on %d", taxon_id)
                continue
            (ensembl_gene_id, external_gene_name,
             description, gene_biotype, entrezgene,
             peptide_id) = row[0:6]
            protein_dict[str(peptide_id)] = ensembl_gene_id
        conn.close()
        return protein_dict

    def fetch_uniprot_gene_map(self, taxon_id):
        """
        Fetch a dict of uniprot-gene  for a species in biomart
        :param taxid:
        :return: dict
        """
        protein_dict = dict()
        params = urllib.parse.urlencode(
            {'query': self._build_biomart_gene_query(str(taxon_id))})
        conn = http.client.HTTPConnection('www.ensembl.org')
        conn.request("GET", '/biomart/martservice?' + params)
        response = conn.getresponse()
        for line in response:
            line = line.decode('utf-8').rstrip()
            row = line.split('\t')
            if len(row) < 7:
                continue
            (ensembl_gene_id, external_gene_name,
             description, gene_biotype, entrezgene,
             peptide_id, uniprot_swissprot) = row[0:7]
            protein_dict[str(uniprot_swissprot)] = ensembl_gene_id
        conn.close()
        return protein_dict

    def _build_biomart_gene_query(self, taxid):
        """
        Building url to fetch equivalent identifiers via Biomart Restful API.
        Documentation at
        http://uswest.ensembl.org/info/data/biomart/biomart_restful.html
        :param taxid:
        :return:

        """
        # basic stuff for ensembl ids
        cols_to_fetch = [
            "ensembl_gene_id", "external_gene_name", "description",
            "gene_biotype", "entrezgene", "ensembl_peptide_id", "uniprotswissprot"]

        if taxid == '9606':
            cols_to_fetch.append("hgnc_id")

        query_attributes = {
            "virtualSchemaName": "default", "formatter": "TSV", "header": "0",
            "uniqueRows": "1", "count": "0", "datasetConfigVersion": "0.6"}

        qry = etree.Element("Query", query_attributes)

        if taxid in self.localtt:
            object_attributes = {"name": self.localtt[taxid], "interface": "default"}
            dataset = etree.SubElement(qry, "Dataset", object_attributes)
            for col in cols_to_fetch:
                etree.SubElement(dataset, "Attribute", {"name": col})
            # is indent right?
            query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>' \
                + etree.tostring(qry, encoding="unicode")
        else:
            query = None

        return query

    def _process_genes(self, taxid, limit=None):
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        model = Model(g)
        geno = Genotype(g)

        raw = '/'.join((self.rawdir, self.files[taxid]['file']))
        line_counter = 0
        logger.info("Processing Ensembl genes for tax %s", taxid)
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t')
            for row in filereader:
                if len(row) < 4:
                    logger.warning("Too few columns in: " + row)
                    raise ValueError("Data error for file %s", raw)
                (ensembl_gene_id, external_gene_name,
                 description, gene_biotype, entrezgene,
                 peptide_id, uniprot_swissprot) = row[0:7]

                # in the case of human genes, we also get the hgnc id,
                # and is the last col
                if taxid == '9606':
                    hgnc_id = row[7]
                else:
                    hgnc_id = None

                if self.testMode and entrezgene != '' \
                        and int(entrezgene) not in self.gene_ids:
                    continue

                line_counter += 1
                gene_id = 'ENSEMBL:' + ensembl_gene_id
                peptide_curie = 'ENSEMBL:{}'.format(peptide_id)
                uniprot_curie = 'UniProtKB:{}'.format(uniprot_swissprot)
                entrez_curie = 'NCBIGene:{}'.format(entrezgene)

                if description == '':
                    description = None

                # gene_type_id = self._get_gene_type(gene_biotype)
                # this had been punted on. ... see what happens
                gene_type_id = self.resolve(gene_biotype.strip(), False)
                if gene_type_id == gene_biotype.strip():
                    gene_type_id = self.globaltt['polypeptide']
                model.addClassToGraph(
                    gene_id, external_gene_name, gene_type_id, description)
                model.addIndividualToGraph(peptide_curie, None, gene_type_id)
                model.addIndividualToGraph(uniprot_curie, None, gene_type_id)

                if entrezgene != '':
                    if taxid == '9606':
                        # Use HGNC for eq in human data
                        model.addXref(gene_id, entrez_curie)
                    else:
                        model.addEquivalentClass(gene_id, entrez_curie)
                if hgnc_id is not None and hgnc_id != '':
                    model.addEquivalentClass(gene_id, hgnc_id)
                geno.addTaxon('NCBITaxon:'+taxid, gene_id)
                if peptide_id != '':
                    geno.addGeneProduct(gene_id, peptide_curie)
                    if uniprot_swissprot != '':
                        geno.addGeneProduct(gene_id, uniprot_curie)
                        model.addXref(peptide_curie, uniprot_curie)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def getTestSuite(self):
        import unittest
        from tests.test_ensembl import EnsemblTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(EnsemblTestCase)

        return test_suite
