#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.GeneReviews import GeneReviews
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class GeneReviewsTestCase(SourceTestCase):

    def setUp(self):
        self.source = GeneReviews('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    # TODO add some specific tests
    # to make sure we are hitting all parts of the code
    # @unittest.skip('test not yet defined')
    # def test_genereviewstest(self):
    #    logger.info("A GeneReviews-specific test")
    #
    #    return


if __name__ == '__main__':
    unittest.main()
