#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.FlyBase import FlyBase
from tests.test_source import SourceTestCase

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class FlyBaseTestCase(SourceTestCase):
    def setUp(self):
        self.source = FlyBase()
        self.source.settestonly(True)
        self.source.setnobnodes(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    def testcvterms(self):
        """
        There's some cvterms we expect to stay constant that are coded against.
        These are included and intended to be compared against their label,
        and will fail if they change.
        :return:
        """
        cvterm_list = {
            131502: 'derived_sf_assoc_alleles', 
            131495: 'derived_tp_assoc_alleles', 
            59983: 'associated_with',
            27: 'producedby',
            60095: 'gets_expression_data_from',
            60384: 'alleleof',
            131520: 'linked_image'
        }

        # TODO check if cvterm file is there first
        import csv
        line_counter = 0
        raw = '/'.join((self.source.rawdir, 'cvterm'))
        logger.info("processing cvterms")
        cvterms_from_file = {}
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                (cvterm_id, cv_id, definition, dbxref_id, is_obsolete,
                 is_relationshiptype, name) = line

                cvterms_from_file[int(cvterm_id)] = name

        for cvterm in cvterm_list:
            self.assertTrue(
                cvterm in cvterms_from_file,
                "cvterm {0} not in file".format(cvterm))
            self.assertTrue(
                cvterm_list[cvterm] == cvterms_from_file[cvterm],
                "cvterm name mismatch (expected:{0}, has:{1})".format(cvterm_list[cvterm],
                                                                                  cvterms_from_file[cvterm]))

        logger.info("All %d cvterms are as expected", len(cvterm_list))

        return

if __name__ == '__main__':
    unittest.main()
