import csv
import logging


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
                 'url': 'http://rest.genome.jp/list/pathway'}
    }

    # I do not love putting these here; but I don't know where else to put them
    test_ids = {
        "pathway": ["path:map00010", "path:map00195", "path:map00100", "ZDB-GENO-050209-30",
                      "ZDB-GENO-051018-1", "ZDB-GENO-070209-80", "ZDB-GENO-070215-11", "ZDB-GENO-070215-12",
                      "ZDB-GENO-070228-3", "ZDB-GENO-070406-1", "ZDB-GENO-070712-5", "ZDB-GENO-070917-2",
                      "ZDB-GENO-080328-1", "ZDB-GENO-080418-2", "ZDB-GENO-080516-8", "ZDB-GENO-080606-609",
                      "ZDB-GENO-080701-2", "ZDB-GENO-080713-1", "ZDB-GENO-080729-2", "ZDB-GENO-080804-4",
                      "ZDB-GENO-080825-3", "ZDB-GENO-091027-1", "ZDB-GENO-091027-2", "ZDB-GENO-091109-1",
                      "ZDB-GENO-100325-3", "ZDB-GENO-100325-4", "ZDB-GENO-100325-5", "ZDB-GENO-100325-6",
                      "ZDB-GENO-100524-2", "ZDB-GENO-100601-2", "ZDB-GENO-100910-1", "ZDB-GENO-111025-3",
                      "ZDB-GENO-120522-18", "ZDB-GENO-121210-1", "ZDB-GENO-130402-5", "ZDB-GENO-980410-268",
                      "ZDB-GENO-080307-1", "ZDB-GENO-960809-7", "ZDB-GENO-990623-3"]
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

