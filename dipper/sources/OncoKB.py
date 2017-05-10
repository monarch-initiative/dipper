import csv
import os
from datetime import datetime
from stat import ST_CTIME
import logging
import re
import shutil
from git import Repo
from git import GitCommandError
from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Model import Model

logger = logging.getLogger(__name__)



class OncoKB(Source):
    """
    OncoKB is  a comprehensive and curated precision oncology knowledge base.
    It offers oncologists detailed, evidence-based information about 
    individual somatic mutations and structural alterations present in 
    patient tumors with the goal of supporting optimal treatment decisions.
    """

    files = {
        'allAnnotatedVars': {
            'file': 'allAnnotatedVariants.txt',
            'url': 'http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'oncokb')
        
        self.dataset = Dataset(
            'oncokb', 'OncoKB Precision Oncology Knowledgebase',
            'http://oncokb.org', None,
            'http://oncokb.org/#/dataAccess')

        self.replaced_id_count = 0

        #if 'test_ids' not in config.get_config()\
        #        or 'disease' not in config.get_config()['test_ids']:
        #    logger.warning("not configured with disease test ids.")
        #    self.test_ids = []
        #else:
        #    self.test_ids = config.get_config()['test_ids']['disease']

        # data-source specific warnings to be removed when issues are cleared
        logger.warning(
            "note that OncoKB provides protein-level mutations with no transcript data.")

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)

        self.scrub()

        # get the latest build from jenkins

        # use the files['version'] file as the version
        fname = '/'.join((self.rawdir, self.files['version']['file']))

        with open(fname, 'r', encoding="utf8") as f:
            # 2015-04-23 13:01
            v = f.readline()  # read the first line (the only line, really)
            d = datetime.strptime(
                v.strip(), '%Y-%m-%d %H:%M').strftime("%Y-%m-%d-%H-%M")
        f.close()

        st = os.stat(fname)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        # this will cause two dates to be attached to the dataset
        # (one from the filedate, and the other from here)
        # TODO when #112 is implemented,
        # this will result in only the whole dataset being versioned
        self.dataset.setVersion(filedate, d)

        self.get_common_files()

        return


    def scrub(self):
        """
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:
        * revise errors in identifiers for some OMIM and PMIDs

        :return: None

        """
        ## TODO -- what needs to be scrubbed?
        return
