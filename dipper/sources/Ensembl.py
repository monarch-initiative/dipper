import logging
import urllib
from urllib import parse
import http
from http import client
import csv

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Genotype import Genotype

from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper import config
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
import xml.etree.ElementTree as etree

logger = logging.getLogger(__name__)


class Ensembl(Source):
    """
    This is the processing module for Ensembl.

    It only includes methods to acquire the equivalences between NCBIGene and ENSG ids using ENSEMBL's Biomart
    services.

    """

    files = {
        '9606': {'file': 'ensembl_9606.txt'},
        '7955': {'file': 'ensembl_7955.txt'},
        '10090': {'file': 'ensembl_10090.txt'}
    }

    def __init__(self, tax_ids=None, gene_ids=None):
        Source.__init__(self, 'ensembl')

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids
        self.load_bindings()

        self.dataset = Dataset('ensembl', 'ENSEMBL', 'http://www.ensembl.org', None)

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]

        self.gene_ids = []
        if 'test_ids' not in config.get_config() or 'gene' not in config.get_config()['test_ids']:
            logger.warn("not configured with gene test ids.")
        else:
            self.gene_ids = config.get_config()['test_ids']['gene']

        self.properties = Feature.properties

        return

    def fetch(self, is_dl_forced=False):

        for t in self.tax_ids:
            logger.info("Fetching genes for %s", str(t))
            loc_file = '/'.join((self.rawdir, 'ensembl_'+str(t)+'.txt'))
            # todo move into util?
            params = urllib.parse.urlencode({'query': self._build_biomart_gene_query(str(t))})
            conn = http.client.HTTPConnection('www.ensembl.org')
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

        self.load_core_bindings()
        self.load_bindings()

        logger.info("Done parsing files.")

        logger.info("Found %d nodes in graph", len(self.graph))
        logger.info("Found %d nodes in testgraph", len(self.testgraph))

        return

    def _build_biomart_gene_query(self, taxid):
        """
        Building url to fetch equivalent identifiers via Biomart Restful API.
        Documentation at http://www.ensembl.org/info/data/biomart/biomart_restful.html
        :param taxid:
        :return:
        """

        cols_to_fetch = ["ensembl_gene_id", "external_gene_name", "description", "gene_biotype",  #basic stuff for ensembl ids
                        "entrezgene", "hgnc_id"]

        query_attributes = {"virtualSchemaName": "default", "formatter": "TSV", "header": "0",
                 "uniqueRows": "1", "count": "0", "datasetConfigVersion": "0.6"}

        ensembl_taxon_to_db_map = {
            '9606' : 'hsapiens_gene_ensembl',
            '10090' : 'mmusculus_gene_ensembl',
            '7955' : 'drerio_gene_ensembl'
        }
        if taxid not in ensembl_taxon_to_db_map:
            logger.error('Bad taxid specified: %s', str(taxid))
            return None
        q = etree.Element("Query", query_attributes)

        object_attributes = {
            "name" : ensembl_taxon_to_db_map[taxid],
            "interface": "default"
        }
        d = etree.SubElement(q, "Dataset", object_attributes)

        for i in cols_to_fetch:
            a = etree.SubElement(d, "Attribute", {"name": i})

        # prepend the query
        prepend = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'

        x = etree.tostring(q, encoding="unicode")
        query = prepend+x

        return query

    def _process_genes(self, taxid, limit=None):
        gu = GraphUtils(curie_map.get())

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)

        raw = '/'.join((self.rawdir, self.files[taxid]['file']))
        line_counter = 0
        logger.info("Processing Ensembl genes for tax %s", taxid)
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t')
            for row in filereader:
                (ensembl_gene_id, external_gene_name, description, gene_biotype,
                 entrezgene, hgnc_id) = row

                if self.testMode and entrezgene != '' and int(entrezgene) not in self.gene_ids:
                    continue

                line_counter += 1
                gene_id = 'ENSEMBL:'+ensembl_gene_id
                if description == '':
                    description = None
                # gene_type_id = self._get_gene_type(gene_biotype)
                gene_type_id = None
                gu.addClassToGraph(g, gene_id, external_gene_name, gene_type_id, description)

                if entrezgene != '':
                    gu.addEquivalentClass(g, gene_id, 'NCBIGene:'+entrezgene)
                if hgnc_id != '':
                    gu.addEquivalentClass(g, gene_id, hgnc_id)
                geno.addTaxon('NCBITaxon:'+taxid, gene_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        gu.loadProperties(g, Feature.object_properties, gu.OBJPROP)
        gu.loadProperties(g, Feature.data_properties, gu.DATAPROP)
        gu.loadProperties(g, Genotype.object_properties, gu.OBJPROP)
        gu.loadAllProperties(g)

        return

    def getTestSuite(self):
        import unittest
        from tests.test_ensembl import EnsemblTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(EnsemblTestCase)

        return test_suite