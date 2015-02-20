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
    '''
    Be sure to have pg user/password connection details in your conf.json file, like:
      dbauth : {
        'disco' : {'user' : '<username>', 'password' : '<password>'}
      }
    '''



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

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")


        logger.info("Finished parsing.")


        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        logger.info("Found %s nodes", len(self.graph))
        return







    #TODO generalize this to a set of utils
    def _getcols(self,cur,table):
        query=(' ').join(("SELECT * FROM",table,"LIMIT 0"))  #for testing
        #print("COMMAND:",query)
        cur.execute(query)
        colnames = [desc[0] for desc in cur.description]
        logger.info("COLS (%s): %s", table, colnames)

        return


    def file_len(self,fname):
        with open(fname) as f:
            return sum(1 for line in f)




    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


