#!/usr/bin/env python3

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from tests import IMPCTests
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    logger.info('Comparing IMPC files to EBI Checksums')
    if IMPCTests.compare_checksums():
        logger.info('PASSED: Files have same checksum as reference')


if __name__ == "__main__":
    main()