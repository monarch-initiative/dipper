#!/usr/bin/env python3

import unittest
import logging
import os
import yaml
from tests import test_general
# from tests import test_dataset
from dipper.utils.GraphUtils import GraphUtils

logging.basicConfig(level=logging.WARNING)
LOG = logging.getLogger(__name__)


class SourceTestCase(unittest.TestCase):
    """
    A testing class for generic source processing functions.
    You would never call these tests directly;
    rather this should be called for any specific source subclasses
    """
    all_test_ids = yaml.safe_load('../resources/test_ids.yaml')

    def setUp(self):
        self.source = None
        # pull in the common test identifiers
        self.all_test_ids = SourceTestCase.all_test_ids
        return

    def tearDown(self):
        self.source = None
        return

    def test_parse(self):
        if self.source is not None:  # don't test the abstract class
            self.source.parse()
            """
            seems we get a better stack trace by not catching the exception
            am I missing something?
            try:
                self.source.parse()
            except Exception as ParseException:  # tec too broad?
                LOG.error(ParseException)
                self.assertFalse(True, "Parsing failed")
            """
            try:
                properties = GraphUtils.get_properties_from_graph(
                    self.source.graph)
                GraphUtils.add_property_axioms(self.source.graph, properties)
                self.source.write()  # default to  fmt='turtle'
                # self.source.write(fmt='nt')
                # self.source.write(fmt='nquads')
            except Exception as WriteException:
                LOG.error(WriteException)
                self.assertFalse(True, "Write failed")

        return

    def test_readGraph(self):
        if self.source is not None:  # don't test the abstract class
            f = self.source.testfile
            self.assertTrue(
                os.path.exists(f), "path does not exist for {0}".format(f))
            test_general.GeneralGraphTestCase().readGraphFromTurtleFile(f)

        return

    @unittest.skip
    def test_readGraphIntoOWL(self):
        if self.source is not None:  # don't test the abstract class
            f = self.source.testfile
            self.assertTrue(os.path.exists(f), "path does not exist for " + f)
            test_general.GeneralGraphTestCase().readGraphIntoOWL(f)

        return

    def _setDirToSource(self):
        if len(os.listdir(self.source.rawdir)) < 1:
            # reset the raw dir to be the source data if it doesn't exist
            # in the test dir
            self.source.rawdir = '../'+self.source.rawdir
            p = os.path.abspath(self.source.rawdir)
            logging.info("Resetting the rawdir to %s", p)
        return

    # no such 'test_ids.json' file exists
    # def _get_conf(self):
    #    if os.path.exists(
    #            os.path.join(os.path.dirname(__file__), 'test_ids.json')):
    #        with open(
    #            os.path.join(
    #                os.path.dirname(__file__), 'test_ids.json')) as json_file:
    #             conf = json.load(json_file)
    #    return conf

    """
    Commenting out as most of our sources do not have licenses
    def test_source_has_license(self):
        if self.source is not None:
            d = self.source.dataset
            test_dataset.DatasetTestCase(d).test_has_license()
        return
    """

if __name__ == '__main__':
    unittest.main()
