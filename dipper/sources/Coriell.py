import logging
import csv

from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Dataset import Dataset
from dipper import config



logger = logging.getLogger(__name__)

class Coriell(Source):
    """
    The Coriell Catalog provided to Monarch includes metadata and descriptions of NIGMS, NINDS, NHGRI, and
    NIA cell lines.  These lines are made available for research purposes.
    Here, we create annotations for the cell lines as models of the diseases from which
    they originate.

    Notice: The Coriell catalog is delivered to Monarch in a specific format, and requires ssh rsa fingerprint
    identification.  Other groups wishing to get this data in it's raw form will need to contact Coriell
    for credentials.

    """

    files = {
        'ninds': {'file': 'NINDS_2014-02-03_13-32-24.csv'},
        'nigms': {'file': 'NIGMS_2014-02-03_13-31-42.csv'},
        'nia': {'file': 'NIA_2015-02-20_16-08-04.csv'}
    }


    def __init__(self):
        Source.__init__(self, 'coriell')

        self.load_bindings()

        self.dataset = Dataset('coriell', 'Coriell', 'http://ccr.coriell.org/nigms/')

        #data-source specific warnings (will be removed when issues are cleared)
        #print()

        #check if config exists; if it doesn't, error out and let user know
        if (not (('keys' in config.get_config()) and ('coriell' in config.get_config()['keys']))):
            logger.error("ERROR: not configured with FTP user/password.")
        return


    def fetch(self, is_dl_forced):
        """
        Be sure to have pg user/password connection details in your conf.json file, like:
        dbauth : {
        "coriell" : {"user" : "<username>", "password" : "<password>", "host" : <host>}
        }

        :param is_dl_forced:
        :return:
        """
        host1 = 'host=\''+ config.get_config()['keys']['coriell']['host']+'\''
        user1 = 'user=\''+ config.get_config()['keys']['coriell']['user']+'\''
        passwd1 = 'passwd=\''+ config.get_config()['keys']['coriell']['password']+'\''
        #print(host1,user1,passwd1)
        #ftp = FTP(config.get_config()['keys']['coriell']['host'],config.get_config()['keys']['coriell']['user'],config.get_config()['keys']['coriell']['password'],timeout=None)
        #ftp = FTP(host1,user1,passwd1,timeout=None)

        #ftp.login()

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:

        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        #pysed.replace("\r", '', ('/').join((self.rawdir,'dv.nlx_157874_1')))

        return

    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")

        for f in ['ninds','nigms']:
            file = ('/').join((self.rawdir,self.files[f]['file']))
            self._process_data(file, limit)

        logger.info("Finished parsing.")


        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        logger.info("Found %s nodes", len(self.graph))
        return

    def _process_data(self, raw, limit=None):
        """
        This function will process the data files from Coriell



        Triples:



        :param raw:
        :param limit:
        :return:
        """
        logger.info("Processing Data from %s",raw)
        #gu = GraphUtils(curie_map.get())

        line_counter = 0
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                line_counter += 1

                (catalog_id,description,omim_number,sample_type,cell_line_available,
                 dna_in_stock,dna_ref,gender,age,race,ethnicity,affected,karyotype,
                 relprob,mutation,gene,fam,collection,url,cat_remark,pubmed_ids) = row



                if (limit is not None and line_counter > limit):
                    break
        return
