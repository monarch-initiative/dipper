#!/usr/bin/env python3

import unittest
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class DatasetTestCase(unittest.TestCase):
    """
    For testing requirements of the Dataset description
    """

    def __init__(self, dataset):
        super().__init__()
        self.dataset = dataset

    def test_has_license(self):

        self.assertTrue(
            self.dataset.get_license() is not None,
            "Source not configured with a license")

        return

    def test_has_version(self):

        self.assertTrue(
            self.dataset.version is not None,
            "Source not configured with a version")

        return

    def test_has_date(self):

        self.assertTrue(
            self.dataset.date_issued is not None,
            "Source not configured with a date")

        return

    # TODO make a testing suite that has all musts
    # TODO make a testing suite that has all desired

if __name__ == '__main__':
    unittest.main()
