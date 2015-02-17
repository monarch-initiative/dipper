from sources.IMPC import IMPC
from utils.TestUtils import TestUtils
import logging

logger = logging.getLogger(__name__)

def compare_checksums():
    """
    test to see if fetched file matches checksum from ebi
    :return: True or False
    """
    is_match = True
    impc = IMPC()
    impc.rawdir = "../raw/impc"
    impc.fetch(False)
    reference_checksums = impc.parse_checksum_file(
        impc.files['checksum']['file'])
    for md5, file in reference_checksums.items():
        test = TestUtils()
        if test.get_file_md5(impc.rawdir, file) != md5:
            is_match = False
            logger.warn('%s was not downloaded completely', file)
            return is_match

    return is_match
