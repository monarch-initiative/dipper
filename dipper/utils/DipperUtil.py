__author__ = 'nlw'


import logging
import re
import unicodedata
import urllib
from urllib import request, parse
from Bio import Entrez
import json


logger = logging.getLogger(__name__)


class DipperUtil:
    """
    Various utilities and quick methods used in this application
    """

    def flatten(self, l):
        """
        Remove None from an array or list
        :param l: An array
        :return:  An array with the None elements removed
        """
        l = list(filter(None.__ne__, l))

        return l

    def remove_control_characters(self, s):
        return "".join(ch for ch in s if unicodedata.category(ch)[0] != "C")


    @staticmethod
    def get_ncbi_taxon_num_by_label(label):
        """
        Here we want to look up the NCBI Taxon id using some kind of label.

        It will only return a result if there is a unique hit.
        :return:
        """
        domain = 'http://eutils.ncbi.nlm.nih.gov'
        path = 'entrez/eutils/esearch.fcgi'

        params = {
            'db':'taxonomy',
            'retmode': 'json',
            'term': label,
        }

        p = urllib.parse.urlencode(params)
        url = '/'.join((domain, path))+'?%s' % p
        logger.info('fetching: %s', url)

        d = urllib.request.urlopen(url)
        resp = d.read().decode()

        myjson = json.loads(resp)

        result = myjson['esearchresult']

        tax_num = None
        if str(result['count']) == '1':
            tax_num = result['idlist'][0]
        else:
            # TODO throw errors
            pass

        return tax_num

    @staticmethod
    def get_homologene_by_gene_num(gene_num):

        Entrez.email = "info@monarchinitiative.org"
        Entrez.tool = "Dipper"
        # first, get the homologene id from the gene id
        # gene_id = '1264'  for testing
        gid = str(gene_num)
        handle = Entrez.esearch(db="homologene",term=gid+"[Gene ID]", retmode="json")
        record = handle.read()
        j = json.loads(record)
        homologene_ids = j["esearchresult"]["idlist"]

        if len(homologene_ids) != 1:
            return

        hid = homologene_ids[0]
        # now, fetch the homologene record
        handle = Entrez.esummary(db="homologene", id=hid, retmode="json")
        record = handle.read()
        j = json.loads(record)

        if 'result' in j and hid in j['result']:
            homologs = j['result'][hid]
        else:
            homologs = None

        return homologs
