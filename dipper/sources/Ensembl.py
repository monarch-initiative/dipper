import logging
import urllib
import csv
import http
import xml.etree.ElementTree as etree

from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.Genotype import Genotype


LOG = logging.getLogger(__name__)
ENS_URL = 'uswest.ensembl.org'    # 'www.ensembl.org'


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
    columns = {
        "bmq_attributes": [  # to be sent
            "ensembl_gene_id",
            "external_gene_name",
            "description",
            "gene_biotype",
            "entrezgene_id",
            "ensembl_peptide_id",
            "uniprotswissprot",
            "hgnc_id",   # human only
        ],
        'bmq_headers': [  # to be recived
            'Gene stable ID',
            'Gene name',
            'Gene description',
            'Gene type',
            'NCBI gene ID',
            'Protein stable ID',
            'UniProtKB/Swiss-Prot ID',
            'HGNC ID',  # human only
        ]
    }

    def __init__(self, graph_type, are_bnodes_skolemized, skip_stats=False, tax_ids=None, gene_ids=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            skip_stats=skip_stats,
            name='ensembl',
            ingest_title='ENSEMBL',
            ingest_url='http://uswest.ensembl.org',
            ingest_logo='https://github.com/monarch-initiative/monarch-ui/blob/master/public/img/sources/source-ensembl.png',
            # license_url=None,
            # data_rights=None,
            # file_handle=None
        )

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = ['9606', '10090', '7955']

        self.tax_ids = [str(x) for x in self.tax_ids]  # for bare ints from commandline

        LOG.info('Have Taxon ID(s)  %s', self.tax_ids)

        self.gene_ids = []
        if 'gene' not in self.all_test_ids:
            LOG.warning("not configured with gene test ids.")
        else:
            self.gene_ids = self.all_test_ids['gene']

        LOG.setLevel(logging.INFO)

    def fetch(self, is_dl_forced=True):  # it is a database query... so no timestamps

        for txid in self.tax_ids:
            LOG.info("Fetching genes for %s", txid)
            loc_file = '/'.join((self.rawdir, 'ensembl_' + txid + '.txt'))

            # todo move into util?
            params = urllib.parse.urlencode(
                {'query': self._build_biomart_gene_query(txid)})
            conn = http.client.HTTPConnection(ENS_URL)
            conn.request("GET", '/biomart/martservice?' + params)
            resp = conn.getresponse()
            with open(loc_file, 'wb') as bin_writer:
                bin_writer.write(resp.read())

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        if self.test_only:
            self.test_mode = True

        LOG.info("Parsing files...")

        for txid in self.tax_ids:
            self._process_genes(txid, limit)

        LOG.info("Done parsing files.")

    # this seems so un specific as to be a bit pointless
    # and will be deleted if there are no objections
    # def fetch_protein_list(self, taxon_id):
    #    """
    #    Fetch a list of proteins for a species in biomart
    #    :param taxid:
    #    :return: list
    #    """
    #    protein_list = list()
    #    att = ['ensembl_peptide_id', ]
    #
    #    params = urllib.parse.urlencode(
    #        {'query': self._build_biomart_gene_query(taxon_id, att)})
    #    conn = http.client.HTTPConnection(ENS_URL)
    #    conn.request("GET", '/biomart/martservice?' + params)
    #    response = conn.getresponse()
    #    header = next(response).rstrip().split('\t')
    #    for line in response:
    #        line = line.decode('utf-8').rstrip()
    #        row = line.split('\t')
    #        if len(row) != len(att):
    #            LOG.warning("Data error for p-list query on %d", taxon_id)
    #            continue
    #        protein_list.append(row[att.index('ensembl_peptide_id')])
    #    conn.close()
    #    return protein_list

    def fetch_protein_gene_map(self, taxon_id):
        """
        Fetch a mapping from proteins to ensembl_gene(S)? for a species in biomart
        :param taxid:
        :return: dict
        """
        # called in StringDB.py
        protein_dict = dict()
        col = ['ensembl_gene_id', 'ensembl_peptide_id']
        col_exp = [
            self.columns['bmq_headers'][
                self.columns['bmq_attributes'].index('ensembl_gene_id')],
            self.columns[
                'bmq_headers'][
                    self.columns['bmq_attributes'].index('ensembl_peptide_id')],
        ]
        raw_query = self._build_biomart_gene_query(taxon_id, col)
        params = urllib.parse.urlencode({'query': raw_query})

        conn = http.client.HTTPConnection(ENS_URL)
        conn.request("GET", '/biomart/martservice?' + params)
        response = conn.getresponse()
        buf = ""
        line_num = 0
        for line in response:
            line = line.decode('utf-8')
            buf = buf + line
            line = line.rstrip()
            row = line.split('\t')
            line_num += 1
            if len(row) != len(col_exp) or row[col.index('ensembl_peptide_id')] == '':
                # LOG.warning('line %i is: %s', line_num, row)
                # ... many rows have no protein id
                pass
            else:
                protein_dict[row[col.index('ensembl_peptide_id')]] = row[
                    col.index('ensembl_gene_id')]

        # hard to know what is there if we never get to check ...
        with open(self.rawdir + '/' + 'gene_prot' + taxon_id + '.resp', 'w') as writer:
            writer.write(str(buf))
        # observed (for human):
        #  - no protien appears more than once
        #  - so no protein could be associated with more than one gene
        #  - so there is no list of genes associated with a protein
        # may have to revisit if other species behave differently
        # (now writing out per species)

        conn.close()
        LOG.info(
            "length gene-protien dict for taxon: %s is %i", taxon_id, len(protein_dict))
        return protein_dict

    def fetch_uniprot_gene_map(self, taxon_id):
        """
        Fetch a dict of uniprot-gene  for a species in biomart
        :param taxid:
        :return: dict
        """
        # unused
        protein_dict = dict()
        col = ['uniprotswissprot', 'ensembl_gene_id']
        params = urllib.parse.urlencode(
            {'query': self._build_biomart_gene_query(taxon_id, col)})
        conn = http.client.HTTPConnection(ENS_URL)
        conn.request("GET", '/biomart/martservice?' + params)
        response = conn.getresponse()
        for line in response:
            line = line.decode('utf-8').rstrip()
            row = line.split('\t')
            if len(row) != len(col):
                continue
            protein_dict[
                row[col.index('uniprotswissprot')]] = row[col.index('ensembl_gene_id')]
        conn.close()
        return protein_dict

    def _build_biomart_gene_query(
            self, taxid, cols_to_fetch=list(columns['bmq_attributes'])):
        """
        Building url to fetch equivalent identifiers via Biomart Restful API.
        Documentation at
        http://uswest.ensembl.org/info/data/biomart/biomart_restful.html
        :param taxid:
        :param array of ensembl biomart attributes to include
        :return:

        """

        if taxid != '9606' and 'hgnc_id' in cols_to_fetch:
            cols_to_fetch.remove('hgnc_id')

        LOG.info('Fetching columns %s', cols_to_fetch)

        query_attributes = {
            "virtualSchemaName": "default", "formatter": "TSV", "header": "1",
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
            LOG.warning("not finding taxon  %s in the local translation table", taxid)
            query = None

        return query

    def _process_genes(self, taxid, limit=None):
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        model = Model(graph)
        geno = Genotype(graph)

        raw = '/'.join((self.rawdir, self.files[taxid]['file']))
        col = list(self.columns['bmq_attributes'])
        if taxid != '9606' and 'hgnc_id' in col:
            col.remove('hgnc_id')
        col_exp = [
            self.columns['bmq_headers'][self.columns['bmq_attributes'].index(x)]
            for x in col]

        LOG.info("Processing Ensembl genes for NCBITaxon:%s", taxid)
        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            row = next(reader)
            if not self.check_fileheader(col_exp, row):
                pass
            for row in reader:
                ensembl_gene_id = row[col.index('ensembl_gene_id')]
                external_gene_name = row[col.index('external_gene_name')]
                description = row[col.index('description')].strip()
                gene_biotype = row[col.index('gene_biotype')].strip()
                entrezgene = row[col.index('entrezgene_id')].strip()
                ensembl_peptide_id = row[col.index('ensembl_peptide_id')].strip()
                uniprotswissprot = row[col.index('uniprotswissprot')].strip()
                hgnc_curie = None
                # in the case of human genes, we also get the hgnc id,
                if taxid == '9606' and 'hgnc_id' in col:
                    hgnc_curie = row[col.index('hgnc_id')].strip()

                if self.test_mode and entrezgene != '' and \
                        entrezgene not in self.gene_ids:
                    continue

                gene_id = 'ENSEMBL:' + ensembl_gene_id
                entrez_curie = 'NCBIGene:{}'.format(entrezgene)

                if description == '':
                    description = None

                gene_type_id = self.resolve(
                    gene_biotype, mandatory=False,
                    default=self.globaltt['polypeptide'])

                model.addClassToGraph(
                    gene_id, external_gene_name, gene_type_id, description)

                if entrezgene != '':
                    if taxid == '9606':
                        # Use HGNC for eq in human data
                        model.addXref(gene_id, entrez_curie)
                    else:
                        model.addEquivalentClass(gene_id, entrez_curie)

                if hgnc_curie is not None and hgnc_curie != '':
                    model.addEquivalentClass(gene_id, hgnc_curie)
                geno.addTaxon('NCBITaxon:' + taxid, gene_id)
                if ensembl_peptide_id is not None and ensembl_peptide_id != '':
                    peptide_curie = 'ENSEMBL:{}'.format(ensembl_peptide_id)
                    model.addIndividualToGraph(peptide_curie, None, gene_type_id)
                    geno.addGeneProduct(gene_id, peptide_curie)
                    if uniprotswissprot != '':
                        uniprot_curie = 'UniProtKB:{}'.format(uniprotswissprot)
                        model.addIndividualToGraph(uniprot_curie, None, gene_type_id)
                        geno.addGeneProduct(gene_id, uniprot_curie)
                        model.addXref(peptide_curie, uniprot_curie)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

    def getTestSuite(self):
        import unittest
        from tests.test_ensembl import EnsemblTestCase
        return unittest.TestLoader().loadTestsFromTestCase(EnsemblTestCase)
