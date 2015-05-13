import csv
import logging
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.OrthologyAssoc import OrthologyAssoc
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


class AnimalQTLdb(Source):

    files = {
        'cattle_btau_bp': {'file': 'QTL_Btau_4.6.gff.txt.gz',
                 'url': 'http://www.animalgenome.org/QTLdb/tmp/QTL_Btau_4.6.gff.txt.gz'},
        'cattle_umd_bp': {'file': 'QTL_UMD_3.1.gff.txt.gz',
                 'url': 'http://www.animalgenome.org/QTLdb/tmp/QTL_UMD_3.1.gff.txt.gz'},
        'cattle_cm': {'file': 'cattle_QTLdata.txt',
                 'url': 'http://www.animalgenome.org/QTLdb/export/KSUI8GFHOT6/cattle_QTLdata.txt'},
        'chicken_bp': {'file': 'QTL_GG_4.0.gff.txt.gz',
                 'url': 'http://www.animalgenome.org/QTLdb/tmp/QTL_GG_4.0.gff.txt.gz'},
        'chicken_cm': {'file': 'chicken_QTLdata.txt',
                 'url': 'http://www.animalgenome.org/QTLdb/export/KSUI8GFHOT6/chicken_QTLdata.txt'},
        'pig_bp': {'file': 'QTL_SS_10.2.gff.txt.gz',
                 'url': 'http://www.animalgenome.org/QTLdb/tmp/QTL_SS_10.2.gff.txt.gz'},
        'pig_cm': {'file': 'pig_QTLdata.txt',
                 'url': 'http://www.animalgenome.org/QTLdb/export/KSUI8GFHOT6/pig_QTLdata.txt'},
        'sheep_bp': {'file': 'QTL_OAR_3.1.gff.txt.gz',
                 'url': 'http://www.animalgenome.org/QTLdb/tmp/QTL_OAR_3.1.gff.txt.gz'},
        'sheep_cm': {'file': 'sheep_QTLdata.txt',
                 'url': 'http://www.animalgenome.org/QTLdb/export/KSUI8GFHOT6/sheep_QTLdata.txt'},
        'horse_bp': {'file': 'QTL_EquCab2.0.gff.txt.gz',
                 'url': 'http://www.animalgenome.org/QTLdb/tmp/QTL_EquCab2.0.gff.txt.gz'},
        'horse_cm': {'file': 'horse_QTLdata.txt',
                 'url': 'http://www.animalgenome.org/QTLdb/export/KSUI8GFHOT6/horse_QTLdata.txt'},
        'rainbow_trout_cm': {'file': 'rainbow_trout_QTLdata.txt',
                 'url': 'http://www.animalgenome.org/QTLdb/export/KSUI8GFHOT6/rainbow_trout_QTLdata.txt'}
    }

    # I do not love putting these here; but I don't know where else to put them
    test_ids = {
    }

    def __init__(self):
        Source.__init__(self, 'animalqtldb')

        # update the dataset object with details about this resource
        # TODO put this into a conf file?
        self.dataset = Dataset('animalqtldb', 'Animal QTL db', 'http://www.genome.jp/kegg/', None, None)

        # source-specific warnings.  will be cleared when resolved.

        return

    def fetch(self, is_dl_forced):
        self.get_files(is_dl_forced)
        #if self.compare_checksums():
            #logger.debug('Files have same checksum as reference')
        #else:
            #raise Exception('Reference checksums do not match disk')
        return

    def parse(self, limit=None):
        """

        :param limit:
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %s rows fo each file", str(limit))

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_cattle_cm(limit)

        logger.info("Finished parsing")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return


    def _process_cattle_cm(self, limit=None):
        """
        This method processes the cattle QTLs in cm format.

        Triples created:

        :param limit:
        :return:
        """

        logger.info("Processing cattle QTLs in cM")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['cattle_cm']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (qtl_id, qtl_symbol, trait_name, assotype, empty, chromosome, position_cm, range_cm,
                 flankmark_a2, flankmark_a1, peak_mark, flankmark_b1, flankmark_b2, exp_id, model, test_base,
                 sig_level, lod_score, ls_mean, p_values, f_statistics, variance, bayes_value, likelihood_ratio,
                 trait_id, dom_effect, add_effect, pubmed_id, gene_id, gene_id_src, gene_id_type, empty2) = row

                #if self.testMode and disease_id not in self.test_ids['disease']:
                    #continue
                #print(row)


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with diseases")
        return
