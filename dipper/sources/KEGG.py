
import logging


from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset


logger = logging.getLogger(__name__)


class KEGG(Source):

    files = {
        'disease': {'file': 'disease',
                 'url': 'http://rest.genome.jp/list/disease'}
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


        logger.info("Finished parsing")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return


