#!/usr/bin/env python3

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from tests import IMPCTests
from tests import SourceTests
from sources.CTD import CTD
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def main():
    logger.info('Comparing IMPC files to EBI Checksums')
    if IMPCTests.compare_checksums():
        logger.info('PASSED: Files have same checksum as reference')
    else:
        logger.warn('FAILED: Reference checksums do not match disk')

    ctd = CTD()
    if SourceTests.compare_local_remote_bytes(ctd):
        logger.info('PASSED: Remote file and local file match')
    else:
        logger.warn('FAILED: Remote file and local file do not match')


if __name__ == "__main__":
    main()