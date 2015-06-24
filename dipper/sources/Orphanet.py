import logging
import urllib
from urllib import request
import re
import time
from datetime import datetime
import json
from subprocess import call
import xml.etree.ElementTree as ET

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.utils.GraphUtils import GraphUtils
from dipper import config
from dipper import curie_map
from dipper.utils.romanplus import romanNumeralPattern, fromRoman, toRoman

logger = logging.getLogger(__name__)


class Orphanet(Source):
    """
    Orphanet’s aim is to help improve the diagnosis, care and treatment of patients with rare diseases.
    For Orphanet, we are currently only parsing the disease-gene associations.


     Note that
    """

    files = {
        'disease-gene': {
            'file': 'en_product6.xml',
            'url': 'http://www.orphadata.org/data/xml/en_product6.xml'}
    }


    def __init__(self):
        Source.__init__(self, 'orphanet')

        self.load_bindings()

        self.dataset = Dataset('orphanet', 'Orphanet', 'http://www.orpha.net',
                               None,
                               'http://creativecommons.org/licenses/by-nd/3.0/',
                               'http://omim.org/help/agreement')


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

        self._process_diseasegene(limit)

        self.load_core_bindings()
        self.load_bindings()

        logger.info("Done parsing.")

        return

    def _process_diseasegene(self, limit):
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

        myfile = '/'.join((self.rawdir, self.files['disease-gene']['file']))

        for event, elem in ET.iterparse(myfile):
            if elem.tag == 'Disorder':
                # get the element name and id
                # id = elem.get('id') # some internal identifier
                disorder_id = 'Orphanet:'+elem.find('OrphaNumber').text
                disorder_label = elem.find('Name').text
                # TODO make the disorder class

                assoc_list = elem.find('DisorderGeneAssociationList')
                for a in assoc_list.findall('DisorderGeneAssociation'):
                    gene_name = a.find('.//Gene/Name').text
                    gene_id = 'Orphanet:'+a.find('./Gene/OrphaNumber').text
                    # TODO make the gene class

                    rel_id = self._map_rel_id(a.find('DisorderGeneAssociationType').get('id'))

                    # TODO make the association
                    # TODO use DisorderGeneAssociationStatus? id="17991" == Assessed
                    # print(disorder_id, gene_name, gene_id, rel_id)

                    rlist = a.find('./Gene/ExternalReferenceList')
                    eqid = None

                    # TODO synonyms? symbol vs name for gene?
                    for r in rlist.findall('ExternalReference'):
                        if r.find('Source').text == 'Ensembl':
                            eqid = 'ENSEMBL:'+r.find('Reference').text
                        elif r.find('Source').text == 'HGNC':
                            eqid = 'HGNC:'+r.find('Reference').text
                        else:
                            pass  # skip the others for now
                        if eqid is not None:
                            # TODO make the equivalence axioms
                            # print(gene_id, '==', eqid)
                            pass
                elem.clear() # discard the element

        return

    def _map_rel_id(self, orphanet_rel_id):
        rel_id = None
        gu = GraphUtils(curie_map.get())
        id_map = {
            '17949': gu.object_properties['has_phenotype'],  # Disease-causing germline mutation(s) in
            '17955': gu.object_properties['has_phenotype'],  # Disease-causing somatic mutation(s) in
            '17961': gu.object_properties['is_marker_for'],  # Major susceptibility factor in
            '17967': gu.object_properties['contributes_to'],  # Modifying germline mutation in
            '17973': gu.object_properties['contributes_to'],  # Modifying somatic mutation in
            '17979': gu.object_properties['contributes_to'],  # Part of a fusion gene in
            '17985': gu.object_properties['contributes_to'],  # Role in the phenotype of
            '18273': None,  # Candidate gene tested in  FIXME
            '25972': gu.object_properties['has_phenotype'],  # Disease-causing germline mutation(s) (loss of function) in
            '25979': gu.object_properties['has_phenotype']   # Disease-causing germline mutation(s) (gain of function) in
        }

        if orphanet_rel_id in id_map:
            rel_id = id_map[orphanet_rel_id]
        else:
            logger.error('Disease-gene association type (%s) not mapped.', orphanet_rel_id)

        return rel_id

    # def getTestSuite(self):
    #     import unittest
    #     from tests.test_omim import OMIMTestCase
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(OMIMTestCase)
    #
    #     return test_suite
