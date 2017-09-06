import logging
import unicodedata
import requests

__author__ = ('nlw', 'tec')
logger = logging.getLogger(__name__)

session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=3)
session.mount('https://', adapter)
session.mount('http://', adapter)

EUTIL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
ESEARCH = EUTIL + '/esearch.fcgi'
ESUMMARY = EUTIL + '/esummary.fcgi'
EREQ = {'email': 'info@monarchinitiative.org', 'tool': 'Dipper'}


class DipperUtil:
    """
    Various utilities and quick methods used in this application

    (A little too quick)
    Per: https://www.ncbi.nlm.nih.gov/books/NBK25497/
    NCBI recommends that users post
     no more than three URL requests per second and
     limit large jobs to either weekends or
     between 9:00 PM and 5:00 AM Eastern time during weekdays

    restructuring to make bulk queries
    is less likely to result in another ban for peppering them with one offs

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
        req = {'db': 'taxonomy', 'retmode': 'json', 'term': label}
        req.update(EREQ)

        request = session.get(ESEARCH, params=req)
        logger.info('fetching: %s', request.url)
        request.raise_for_status()
        result = request.json()['esearchresult']

        # Occasionally eutils returns the json blob
        # {'ERROR': 'Invalid db name specified: taxonomy'}
        if 'ERROR' in result:
            request = session.get(ESEARCH, params=req)
            logger.info('fetching: %s', request.url)
            request.raise_for_status()
            result = request.json()['esearchresult']

        tax_num = None
        if 'count' in result and str(result['count']) == '1':
            tax_num = result['idlist'][0]
        else:
            # TODO throw errors
            logger.warning(
                'ESEARCH for taxon label "%s"  returns %s', label, str(result))
            pass

        return tax_num

    @staticmethod
    def get_homologene_by_gene_num(gene_num):
        # first, get the homologene id from the gene id
        # gene_id = '1264'  for testing
        req = {
            'db': 'homologene',
            'retmode': 'json',
            'term': str(gene_num) + "[Gene ID]"
            # 'usehistory': None  # y/n
            # 'rettype':  [uilist|count]
        }
        req.update(EREQ)
        r = requests.get(ESEARCH, params=req)
        homologene_ids = r.json()['esearchresult']['idlist']
        if len(homologene_ids) != 1:
            return
        hid = homologene_ids[0]
        # now, fetch the homologene record
        req = {'db': 'homologene', 'id': hid, 'retmode': 'json'}
        req.update(EREQ)
        request = session.get(ESUMMARY, params=req)
        data = request.json()

        if 'result' in data and hid in data['result']:
            homologs = data['result'][hid]
        else:
            homologs = None

        return homologs

    @staticmethod
    def is_omim_disease(gene_id):
        """
        Process omim equivalencies by examining the monarch ontology scigraph
        As an alternative we could examine mondo.owl, since the ontology
        scigraph imports the output of this script which creates an odd circular
        dependency (even though we're querying mondo.owl through scigraph)

        :param graph: rdfLib graph object
        :param gene_id: ncbi gene id as curie
        :param omim_id: omim id as curie
        :return: None
        """
        SCIGRAPH_BASE = 'https://scigraph-ontology-dev.monarchinitiative.org/scigraph/graph/'

        session = requests.Session()
        adapter = requests.adapters.HTTPAdapter(max_retries=10)
        session.mount('https://', adapter)

        requests_log = logging.getLogger("requests.packages.urllib3")
        requests_log.setLevel(logging.ERROR)

        isOmimDisease = False
        url = SCIGRAPH_BASE + gene_id + '.json'
        response = session.get(url)
        try:
            results = response.json()
            if 'nodes' in results and len(results['nodes']) > 0:
                if 'meta' in results['nodes'][0] \
                        and 'category' in results['nodes'][0]['meta'] \
                        and 'disease' in results['nodes'][0]['meta']['category']:
                    logger.info("{0} is a disease, skipping".format(gene_id))
                    isOmimDisease = True
        except ValueError:
            pass

        return isOmimDisease

    @staticmethod
    def get_ncbi_id_from_symbol(gene_symbol):
        """
        Get ncbi gene id from symbol using monarch and mygene services
        :param gene_symbol:
        :return:
        """
        monarch_url = 'https://solr.monarchinitiative.org/solr/search/select'
        params = DipperUtil._get_solr_weight_settings()
        params["q"] = "{0} \"{0}\"".format(gene_symbol)
        params["fq"] = ["taxon:\"NCBITaxon:9606\"", "category:\"gene\""]
        gene_id = None
        try:
            monarch_request = requests.get(monarch_url, params=params)
            response = monarch_request.json()
            count = response['response']['numFound']
            if count > 0:
                gene_id = response['response']['docs'][0]['id']
        except requests.ConnectionError:
            print("error fetching {0}".format(monarch_url))

        return gene_id

    @staticmethod
    def _get_solr_weight_settings():
        return {
            "qt": "standard",
            "json.nl": "arrarr",
            "fl": "*,score",
            "start": "0",
            "rows": "5",
            "defType": "edismax",
            "personality": "monarch_search",
            "qf": [
                "label_searchable^1",
                "definition_searchable^1",
                "synonym_searchable^1",
                "label_std^2",
                "synonym_std^1"
            ],
            "wt": "json"
        }
