import logging
import unicodedata
import requests

__author__ = ('nlw', 'tec')
logger = logging.getLogger(__name__)

session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)
session.mount('http://', adapter)


class DipperUtil:
    """
    Various utilities and quick methods used in this application

    """

    def remove_control_characters(self, s):
        '''
        Filters out charcters in any of these unicode catagories
        [Cc] 	Other, Control      ( 65 characters) \n,\t ...
        [Cf] 	Other, Format       (151 characters)
        [Cn] 	Other, Not Assigned (  0 characters -- none have this property)
        [Co] 	Other, Private Use  (  6 characters)
        [Cs] 	Other, Surrogate    (  6 characters)
        '''
        return "".join(ch for ch in s if unicodedata.category(ch)[0] != "C")

    @staticmethod
    def get_ncbi_taxon_num_by_label(label):
        """
        Here we want to look up the NCBI Taxon id using some kind of label.
        It will only return a result if there is a unique hit.
        :return:

        """
        domain = 'http://eutils.ncbi.nlm.nih.gov'
        path = '/entrez/eutils/esearch.fcgi'
        req = {
            'db': 'taxonomy',
            'retmode': 'json',
            'term': label}

        request = session.get(domain + path, params=req)
        logger.info('fetching: %s', request.url)
        request.raise_for_status()
        result = request.json()['esearchresult']

        tax_num = None
        if str(result['count']) == '1':
            tax_num = result['idlist'][0]
        else:
            # TODO throw errors
            pass

        return tax_num

    @staticmethod
    def get_homologene_by_gene_num(gene_num):
        # first, get the homologene id from the gene id
        # gene_id = '1264'  for testing
        host = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
        esearch = host + '/esearch.fcgi'
        esummary = host + '/esummary.fcgi'
        req = {
            'email': 'info@monarchinitiative.org',
            'tool': 'Dipper',
            'db': 'homologene',
            'retmode': 'json',
            'term': str(gene_num) + "[Gene ID]"
            # 'usehistory': None  # y/n
            # 'rettype':  [uilist|count]
        }
        r = requests.get(esearch, params=req)
        homologene_ids = r.json()['esearchresult']['idlist']
        if len(homologene_ids) != 1:
            return
        hid = homologene_ids[0]
        # now, fetch the homologene record
        req = {'db': 'homologene', 'id': hid, 'retmode': 'json'}
        request = session.get(esummary, params=req)
        data = request.json()

        if 'result' in data and hid in data['result']:
            homologs = data['result'][hid]
        else:
            homologs = None

        return homologs
