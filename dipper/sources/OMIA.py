import logging
import xml.etree.ElementTree as ET

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.utils.GraphUtils import GraphUtils
from dipper import config
from dipper import curie_map


logger = logging.getLogger(__name__)


class OMIA(Source):
    """
    """

    files = {
        'data': {
            'file': 'omia.xml',
            'url': 'http://omia.angis.org.au/dumps/omia.xml.gz'}  # TODO change to gz
    }

    def __init__(self):
        Source.__init__(self, 'omia')

        self.load_bindings()

        self.dataset = Dataset('omia', 'Online Mendelian Inheritance in Animals', 'http://omia.angis.org.au',
                                None,
                                None,
                                'http://sydney.edu.au/disclaimer.shtml')

        # check to see if there's any ids configured in the config; otherwise, warn
        if 'test_ids' not in config.get_config() or 'disease' not in config.get_config()['test_ids']:
            logger.warn("not configured with disease test ids.")

        return

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced:
        :return:
        """
        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        # first, process the species, breeds, genes, articles, and other static stuff

        # next process the association data

        self._process_associations(limit)

        self.load_core_bindings()
        self.load_bindings()

        logger.info("Done parsing.")

        return

    def _process_associations(self, limit):
        """
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())

        # names of tables to iterate - probably don't need all these:
        # Article_Breed, Article_Keyword, Article_Gene, Article_Keyword, Article_People, Article_Phene,
        # Articles, Breed, Breed_Phene, Genes_gb, Group_Categories, Group_MPO, Inherit_Type, Keywords,
        # Landmark, Lida_Links, OMIA_Group, OMIA_author, Omim_Xref, People, Phene, Phene_Gene,
        # Publishers, Resources, Species_gb, Synonyms

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = open(myfile, 'r', errors='replace', encoding='utf-8')
        f.readline()  # remove the xml declaration line

        parser = ET.XMLParser(encoding='utf-8')
        # if i use utf-8 encoding, i get to line 1510906 - but there's some "control-B" chars to clean up

        for event, elem in ET.iterparse(myfile, parser=parser):
            if elem.tag == 'table_data':
                # get the element name and id
                # id = elem.get('id') # some internal identifier
                article_breed_data = elem.find("[@name='Article_Breed']")
                breed_data = elem.find("[@name='Breed']")
                article_phene = elem.find("[@name='Article_Phene']")
                phene_data = elem.find("[@name='Phene']")
                # process each row
                if article_breed_data is not None:
                    row = {
                        'article_id': None,
                        'breed_id': None,
                        'added_by': None
                    }
                    for r in article_breed_data.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text
                            # make a mentions between the breed and the article

                        # TODO look up the keys to get the actual identifiers
                        article_id = 'PMID:'+str(row['article_id'])
                        breed_id = self._make_breed_id(row['breed_id'])
                        gu.addTriple(g, article_id, gu.object_properties['is_about'], breed_id)

                elif breed_data is not None:
                    row = {
                        'breed_id': None,
                        'breed_name': None,
                        'gb_species_id': None,
                        'added_by': None,
                        'date_modified': None
                    }
                    for r in breed_data.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        breed_id = self._make_breed_id(row['breed_id'])
                        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
                        # todo add the common name of species in parentheses of name
                        gu.addIndividualToGraph(g, breed_id, row['breed_name'], tax_id)

                elif phene_data is not None:
                    row = {}
                    self.phene_hash = {}
                    for r in phene_data.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        if 'omia_id' not in row:
                            logger.info("omia_id not present for %s", row['phene_id'])
                        omia_id = 'OMIA:'+str(row['omia_id'])
                        omia_label = row['phene_name']
                        if omia_label == '':
                            omia_label = None
                        self.phene_hash[row['phene_id']] = omia_id  # add to internal hash store for later lookup
                        descr = row['summary']
                        if descr == '':
                            descr = None
                        gu.addClassToGraph(g, omia_id, omia_label, None, descr)  # add as subclass of disease?
                        if row['symbol'] is not None:
                            gu.addSynonym(g, omia_id, row['symbol'])
                        # TODO add taxon info/restriction?

                    elem.clear()  # discard the element

            if self.testMode and limit is not None and line_counter > limit:
                return

        return

    def _make_breed_id(self, num):
        # todo make this a proper class
        breed_id = '_breed-'+str(num)
        if self.nobnodes:
            breed_id = ':'+breed_id

        return breed_id

    # def getTestSuite(self):
    #     import unittest
    #     from tests.test_orphanet import OMIATestCase
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(OMIATestCase)
    #
    #     return test_suite
