import os
from datetime import datetime
from stat import *
import re
import logging
from ftplib import FTP


from rdflib import Literal
from rdflib.namespace import DC, FOAF
from rdflib import URIRef

from sources.Source import Source
from models.Assoc import Assoc
from models.Dataset import Dataset
from utils.CurieUtil import CurieUtil
from conf import config, curie_map
from utils.GraphUtils import GraphUtils

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
        host1 = 'host=\''+config.get_config()['keys']['coriell']['host']+'\''
        user1 = 'user=\''+config.get_config()['keys']['coriell']['user']+'\''
        passwd1 = 'passwd=\''+config.get_config()['keys']['coriell']['password']+'\''
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


        logger.info("Finished parsing.")


        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        logger.info("Found %s nodes", len(self.graph))
        return


