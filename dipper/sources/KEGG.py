import csv
import logging
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


class KEGG(Source):

    files = {
        'disease': {'file': 'disease',
                 'url': 'http://rest.genome.jp/list/disease'},
        'pathway': {'file': 'pathway',
                 'url': 'http://rest.genome.jp/list/pathway'},
        'genes': {'file': 'genes',
                 'url': 'http://rest.genome.jp/list/hsa'},
        'orthology': {'file': 'orthology',
                 'url': 'http://rest.genome.jp/list/orthology'},
        'disease_gene': {'file': 'disease_gene',
                 'url': 'http://rest.kegg.jp/link/disease/hsa'}
    }

    # I do not love putting these here; but I don't know where else to put them
    test_ids = {
        "pathway": ["path:map00010", "path:map00195", "path:map00100", "path:map00340"],
        "disease": ["ds:H00015", "ds:H00026", "ds:H00712", "ds:H00736"],
        "genes": ["hsa:100506275", "hsa:285958", "hsa:286410", "hsa:6387"],
        "orthology": ["ko:K00010", "ko:K00027", "ko:K00042", "ko:K00088"]
    }

    def __init__(self):
        Source.__init__(self, 'kegg')

        # update the dataset object with details about this resource
        # TODO put this into a conf file?
        self.dataset = Dataset('kegg', 'KEGG', 'http://www.genome.jp/kegg/', None, None)

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

        self._process_pathways(limit)
        self._process_diseases(limit)
        self._process_genes(limit)
        self._process_disease2gene(limit)




        logger.info("Finished parsing")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return



    def _process_pathways(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing pathways")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['pathway']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pathway_id, pathway_name) = row

                if self.testMode and pathway_id not in self.test_ids['pathway']:
                    continue

                pathway_id = 'KEGG:'+pathway_id.strip()
                # Add the pathway as a class.
                gu.addClassToGraph(g, pathway_id, pathway_name)


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with pathways")
        return

    def _process_diseases(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing diseases")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (disease_id, disease_name) = row

                if self.testMode and disease_id not in self.test_ids['disease']:
                    continue

                disease_id = 'KEGG:'+disease_id.strip()
                # Add the disease as a class.
                gu.addClassToGraph(g, disease_id, disease_name)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with diseases")
        return

    def _process_genes(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing genes")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['genes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, gene_name) = row

                if self.testMode and gene_id not in self.test_ids['genes']:
                    continue

                gene_id = 'KEGG:'+gene_id.strip()
                # Add the disease as a class.
                gu.addClassToGraph(g, gene_id, gene_name)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with genes")
        return

    def _process_orthology(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing orthology")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['orthology']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (orthology_id, gene_name) = row

                if self.testMode and gene_id not in self.test_ids['orthology']:
                    continue

                gene_id = 'KEGG:'+gene_id.strip()
                # Add the disease as a class.
                #gu.addClassToGraph(g, gene_id, gene_name)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with orthology")
        return



    def _process_disease2gene(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing disease to gene")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['disease_gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, disease_id) = row

                if self.testMode and gene_id not in self.test_ids['']:
                    continue

                gene_id = 'KEGG:'+gene_id.strip()
                disease_id = 'KEGG:'+disease_id.strip()

                # Add the disease to gene relationship.


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with disease to gene")
        return