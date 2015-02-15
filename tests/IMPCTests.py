#!/usr/bin/env python3

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from sources.IMPC import IMPC
from utils.TestUtils import TestUtils
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    compare_checksums()

def compare_checksums():
    """
    test to see if fetched file matches checksum from ebi
    :return: True or False
    """
    logger.info('Comparing IMPC files to EBI Checksums')
    is_foo = True
    impc = IMPC()
    impc.fetch(False)
    reference_checksums = impc.parse_checksum_file(
        impc.files['checksum']['file'])
    for md5, file in reference_checksums.items():
        test = TestUtils()
        if test.get_file_md5(impc.rawdir, file) != md5:
            is_match = False
            logger.info('FAILED: ' + file +
                        ' was not downloaded completely')
            return is_match

    logger.info('PASSED: Files have same checksum as reference')

if __name__ == "__main__":
    main()


