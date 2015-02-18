#!/usr/bin/env python3

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from models.Genotype import Genotype
import unittest
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class GenotypeTestCase(unittest.TestCase):

    def setUp(self):
        genotype_id = 'MGI:5515892'
        genotype_label = 'Pmp22<Tr-2J>/Pmp22<+> [C57BL/6J-Pmp22<Tr-2J>/GrsrJ]'
        self.genotype = Genotype('The widget')

    def tearDown(self):
        self.genotype.dispose()
        self.genotype = None


if __name__ == '__main__':
    unittest.main()
