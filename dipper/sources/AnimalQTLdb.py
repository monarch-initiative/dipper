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

        logger.info("Processing QTLs in cM")
        self._process_QTLs_genetic_location(('/').join((self.rawdir, self.files['cattle_cm']['file'])), 'AQTLCattle:', 'AQTLTraitCattle:', 'NCBITaxon:9913', limit)
        self._process_QTLs_genetic_location(('/').join((self.rawdir, self.files['chicken_cm']['file'])), 'AQTLChicken:', 'AQTLTraitChicken:', 'NCBITaxon:9031', limit)
        self._process_QTLs_genetic_location(('/').join((self.rawdir, self.files['pig_cm']['file'])), 'AQTLPig:', 'AQTLTraitPig:', 'NCBITaxon:9823', limit)
        self._process_QTLs_genetic_location(('/').join((self.rawdir, self.files['sheep_cm']['file'])), 'AQTLSheep:', 'AQTLTraitSheep:', 'NCBITaxon:9940', limit)
        self._process_QTLs_genetic_location(('/').join((self.rawdir, self.files['horse_cm']['file'])), 'AQTLHorse:', 'AQTLTraitHorse:', 'NCBITaxon:9796', limit)
        self._process_QTLs_genetic_location(('/').join((self.rawdir, self.files['rainbow_trout_cm']['file'])), 'AQTLRainbowTrout:', 'AQTLTraitRainbowTrout:', 'NCBITaxon:8022', limit)

        # TODO: Need to bring in the Animal QTL trait map?
        logger.info("Finished parsing")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return


    #TODO: Abstract this into a general function
    # Need to pass in: file, qtl prefix, trait prefix, taxon,

    def _process_QTLs_genetic_location(self, raw, qtl_prefix, trait_prefix, taxon_id, limit=None):
        """
        This method processes the cattle QTLs in cm format.

        Triples created:

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
        #raw = ('/').join((self.rawdir, self.files['cattle_cm']['file']))
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

                #FIXME: Not sure that I like these prefixes. Is there a better approach?
                qtl_id = qtl_prefix+qtl_id
                trait_id = trait_prefix+trait_id

                #FIXME: For assotype, the QTL is indicated either as a QTL or an Association.
                # Should Associations be handled differently?

                # Add QTL to graph
                gu.addIndividualToGraph(g, qtl_id, qtl_symbol, geno.genoparts['QTL'])

                geno.addTaxon(taxon_id,qtl_id)
                # Add trait to graph as a phenotype - QTL has phenotype?


                if re.match('ISU.*', pubmed_id):
                    pub_id = 'AQTLPub:'+pubmed_id
                else:
                    pub_id = 'PMID:'+pubmed_id

                # Add publication
                gu.addIndividualToGraph(g,pub_id,None)
                eco_id = "ECO:0000059"  # Using experimental phenotypic evidence
                assoc_id = self.make_id((qtl_id+trait_id+pub_id))
                assoc = G2PAssoc(assoc_id, qtl_id, trait_id, pub_id, eco_id)
                assoc.addAssociationNodeToGraph(g)

                # Add gene to graph,

                # Add cm data as location?

                # Add publication



                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with diseases")
        return
