import csv
import os
from datetime import datetime
from stat import *
import logging

from dipper.utils import pysed
from dipper.utils.GraphUtils import GraphUtils
from dipper.sources.Source import Source
from dipper.models.D2PAssoc import D2PAssoc
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.DispositionAssoc import DispositionAssoc
from dipper.models.Dataset import Dataset
from dipper.models.Assoc import Assoc
from dipper import curie_map
from dipper import config


logger = logging.getLogger(__name__)


class MMRRC(Source):
    """
    The [Human Phenotype Ontology](http://human-phenotype-ontology.org) group curates and assembles
    over 115,000 annotations to hereditary diseases using the HPO ontology.
    Here we create OBAN-style associations between diseases and phenotypic features, together with their
    evidence, and age of onset and frequency (if known).
    The parser currently only processes the "abnormal" annotations.  Association to "remarkable normality"
    will be added in the near future.
    """

    files = {
        'catalog': {'file' : 'mmrrc_catalog_data.csv',
                   'url' : 'http://www.mmrrc.org/about/mmrrc_catalog_data.csv'},
    }


    def __init__(self):
        Source.__init__(self, 'mmrrc')

        self.load_bindings()

        self.dataset = Dataset('mmrrc', 'Mutant Mouse Regional Resource Centers',
                               'https://www.mmrrc.org', None,
                               None)

        if 'test_ids' not in config.get_config() and 'gene' not in config.get_config()['test_ids']:
            logger.warn("not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['gene']

        return

    def fetch(self, is_dl_forced=False):

        #self.get_files(is_dl_forced)
        fname = '/'.join((self.rawdir,self.files['catalog']['file']))
        st = os.stat(fname)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        #TODO note can set the data version to what is in the header
        # first line like: This MMRRC catalog data file was generated on 2015-04-22

        self.dataset.setVersion(filedate)

        return


    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_phenotype_data( limit)


        logger.info("Finished parsing.")

        return


    def _process_phenotype_data(self, limit):
        '''
STRAIN/STOCK_ID,STRAIN/STOCK_DESIGNATION,STRAIN_TYPE,STATE,MGI_ALLELE_ACCESSION_ID,ALLELE_SYMBOL,ALLELE_NAME,MUTATION_TYPE,CHROMOSOME,MGI_GENE_ACCESSION_ID,GENE_SYMBOL,GENE_NAME,SDS_URL,MPT_IDS,ACCEPTED_DATE
"000001-UNC","C57BL/6-Tg(Fga,Fgb,Fgg)1Unc/Mmnc","MSR","","MGI:95526","Fgg","fibrinogen gamma chain","TG","3","MGI:95526","Fgg","fibrinogen gamma chain","http://www.mmrrc.org/strains/1.php","1764 1790 2419 2723 5376 5387 5416 8469 9642 10211 10213","2001-05-01"

        # NOTE: If a Strain carries more than one mutation, each Mutation description,
        i.e., the set: (Mutation Type - Chromosome - Gene Symbol - Gene Name - Allele Symbol - Allele Name)
        will require a separate line.


        :param raw:
        :param limit:
        :return:
        '''
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        gu = GraphUtils(curie_map.get())
        fname = '/'.join((self.rawdir,self.files['catalog']['file']))
        with open(fname, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for row in filereader:
                line_counter += 1
                # skip the first line, which is the version info
                # and skip the second line, which is the header
                if line_counter < 3:
                    continue

                (strain_num, strain_label, strain_type_symbol, strain_state,
                 mgi_allele_id, mgi_allele_symbol, mgi_allele_name, mutation_type, chrom,
                 mgi_gene_id, mgi_gene_symbol, mgi_gene_name, sds_url, mp_nums, accepted_date) = row

                if self.testMode and (mgi_gene_id not in self.test_ids):
                    continue

                #"000001-UNC","C57BL/6-Tg(Fga,Fgb,Fgg)1Unc/Mmnc","MSR","","MGI:95526","Fgg","fibrinogen gamma chain","TG","3",\
                #"MGI:95526","Fgg","fibrinogen gamma chain","http://www.mmrrc.org/strains/1.php","1764 1790 2419 2723 5376 5387 5416 8469 9642 10211 10213","2001-05-01"

                phenotype_ids = []
                # split apart the mp ids
                if mp_nums != '':
                    for i in mp_nums.split(' '):
                        mp_id = 'MP:'+i.zfill(7)
                        phenotype_ids.append(mp_id)

                # TODO the contents of the strain_num after the dash are the "holding center".  strip this out.
                # should map to https://www.mmrrc.org/catalog/sds.php?mmrrc_id=36643

                strain_id = ':'.join(('MMRRC',strain_num))
                # TODO add strain_type
                gu.addIndividualToGraph(g,strain_id,strain_label)
                # TODO 1.  build a genotype from the alleles
                # TODO 2.  pull the background out of the strain

                # TODO 3.  create associations for each of the phenotypes with the strain

                for pid in phenotype_ids:
                    gu.addClassToGraph(g,pid,None)   #assume the phenotype label is in the ontology
                    assoc_id = self.make_id(('mmrrc'+strain_id+pid))
                    assoc = G2PAssoc(assoc_id,strain_id,pid,None,None)
                    assoc.addAssociationNodeToGraph(g)

                # is the date for the associations, or for the strain?

                if not self.testMode and (limit is not None and line_counter > limit):
                    break

            Assoc().loadAllProperties(g)

        return

    def getTestSuite(self):
        import unittest
        from tests.test_mmrrc import MMRRCTestCase
        # TODO add G2PAssoc tests

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MMRRCTestCase)

        return test_suite
