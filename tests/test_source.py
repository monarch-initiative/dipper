#!/usr/bin/env python3

from dipper.sources.Source import Source
from tests import test_general

import unittest
import logging
import os
import json

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class SourceTestCase(unittest.TestCase):
    """
    A testing class for generic source processing functions.
    You would never call these tests directly; rather this should be called for any specific source subclasses
    """

    def setUp(self):
        self.source = None
        return

    def tearDown(self):
        self.source = None
        return

    def test_parse(self):
        # TODO figure out how to skip this if we are running this from the source itself
        if self.source != Source():  # don't test the abstract class
            try:
                self.source.parse()
                self.assertTrue(True)
                self.source.write(format='turtle')
                self.assertTrue(True)
            except:
                self.assertFalse(False, "Parsing failed")

        return

    def test_readGraph(self):
        if self.source is not None:  # don't test the abstract class
            f = self.source.testfile
            p = os.path.abspath(f)
            self.assertTrue(os.path.exists(f), "path does not exist for {0}".format(f))
            test_general.GeneralGraphTestCase().readGraphFromTurtleFile(f)

        return

    @unittest.skip
    def test_readGraphIntoOWL(self):
        if self.source is not None:  # don't test the abstract class
            f = self.source.testfile
            p = os.path.abspath(f)
            self.assertTrue(os.path.exists(f), "path does not exist for "+f)
            test_general.GeneralGraphTestCase().readGraphIntoOWL(f)

        return

    def _setDirToSource(self):
        if len(os.listdir(self.source.rawdir)) < 1:
            # reset the raw dir to be the source data if it doesn't exist in the test dir
            self.source.rawdir = '../'+self.source.rawdir
            p = os.path.abspath(self.source.rawdir)
            logging.info("Resetting the rawdir to %s", p)
        return

    def _get_conf(self):
        if os.path.exists(os.path.join(os.path.dirname(__file__), 'test_ids.json')):
            with open(os.path.join(os.path.dirname(__file__),
                           'test_ids.json')) as json_file:
                conf = json.load(json_file)
        return conf


if __name__ == '__main__':
    unittest.main()