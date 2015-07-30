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

        self.assertTrue(self.dataset.get_license() is not None, "Source not configured with a license")

        return

    # TODO make a testing suite that has all musts
    # TODO make a testing suite tha has all desired

if __name__ == '__main__':
    unittest.main()
