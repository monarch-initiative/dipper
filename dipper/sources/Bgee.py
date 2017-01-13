import logging
import os
from stat import ST_CTIME, ST_SIZE
import re
from zipfile import ZipFile
from datetime import datetime
import time
from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
import ftplib


logger = logging.getLogger(__name__)


class Bgee(Source):
    """
    Bgee is a database to retrieve and compare gene expression
    patterns between animal species.

    Bgee first maps heterogeneous expression data (currently RNA-Seq,
    Affymetrix, in situ hybridization, and EST data) to anatomy and
    development of different species.

    Then, in order to perform automated cross species comparisons,
    homology relationships across anatomies, and comparison criteria
    between developmental stages, are designed.
    """
    BGEE_FTP = 'ftp.bgee.org'
    ANAT_ENTITY_DIR = 'current/download/ranks/anat_entity/'
    DEFAULT_TAXA = ['10090', '10116', '13616', '28377', '6239',
                    '7227', '7955', '8364', '9031', '9258',
                    '9544', '9593', '9597', '9598', '9606',
                    '9823', '9913']

    def __init__(self, tax_ids=None):
        """
        :param tax_ids: List of taxon
        :return:
        """
        super().__init__('bgee')
        if tax_ids is None:
            self.tax_ids = Bgee.DEFAULT_TAXA
        else:
            self.tax_ids = tax_ids
        self.load_bindings()

        self.dataset = Dataset(
            'Bgee', 'Bgee Gene expression data in animals',
            'http://bgee.org/')

    def fetch(self, is_dl_forced=False):
        """
        :return: None
        """
        files_to_download, ftp = self._get_file_list()
        for name, info in files_to_download:
            localfile = '/'.join((self.rawdir, name))
            if not os.path.exists(localfile)\
                    or self.checkIfRemoteIsNewer(localfile, info['size'],
                                                 info['modify']):
                logger.info("Fetching {}".format(name))
                logger.info("Writing to {}".format(localfile))
                ftp.retrbinary('RETR {}'.format(name), open(localfile, 'wb').write)
                remote_dt = Bgee._convert_ftp_time_to_utc(info['modify'])
                os.utime(localfile, (time.mktime(remote_dt.timetuple()),
                                     time.mktime(remote_dt.timetuple())))

        ftp.quit()

        return

    def parse(self, limit=None):
        """
        Given the input taxa, expects files in the raw directory
        with the name {tax_id}_anat_entity_all_data_Pan_troglodytes.tsv.zip
        :param tax_ids: Limit to top ranked anatomy associations
        :return: None
        """
        return

    # Override
    @staticmethod
    def checkIfRemoteIsNewer(localfile, remote_size, remote_modify):
        """
        Overrides checkIfRemoteIsNewer in Source class
        :param localfile:
        :param remote_size:
        :param remote_modify:
        :return:
        """
        is_remote_newer = False
        status = os.stat(localfile)
        logger.info("Local file date: {0}, size: {1}".format(
                    datetime.fromtimestamp(status.st_atime),
                    status[ST_SIZE]))
        remote_dt = Bgee._convert_ftp_time_to_utc(remote_modify)

        if remote_dt != datetime.fromtimestamp(status.st_atime) \
                or status[ST_SIZE] != int(remote_size):
            is_remote_newer = True
            logger.info(
                "Object on server is has different size {0} and/or "
                "date {1}".format(remote_size, remote_dt))

        return is_remote_newer

    @staticmethod
    def _convert_ftp_time_to_utc(dt):
        """
        Convert datetime in the format 20160705042714 to a datetime object
        :return: Tuple of (Generator object with Tuple(file name, info object), ftp object)
        """
        date_time = datetime(int(dt[:4]), int(dt[4:6]),
                             int(dt[6:8]), int(dt[8:10]),
                             int(dt[10:12]), int(dt[12:14]))
        return date_time

    def _get_file_list(self):
        """
        Get file list from ftp server filtered by taxon
        :return: Tuple of (Generator object with Tuple(file name, info object), ftp object)
        """
        ftp = ftplib.FTP(Bgee.BGEE_FTP)
        ftp.login("anonymous", "info@monarchinitiative.edu")
        ftp.cwd(Bgee.ANAT_ENTITY_DIR)

        directory = ftp.mlsd()
        files = (value for value in directory if value[1]['type'] == 'file')
        files_to_download = (value for value in files
                             if re.match(r'.*_all_data_.*', value[0], )
                             and re.findall(r'^\d+', value[0])[0] in self.tax_ids)

        return files_to_download, ftp