#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.ZFIN import ZFIN

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class ZFINTestCase(SourceTestCase):

    def setUp(self):
        self.source = ZFIN('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return

    @unittest.skip(
        'Will eventually write test to check if phenotype sextuples' +
        'are mapped to ZP ids')
    def test_allZPAvailable(self):
        """
        This test will identify if there are
        any missing ZP terms in the mapping file
        :return:

        """
        # TODO add this test to check if all phenotype sextuples
        # are mapped to ZP ids

        return

    def test_load_zp_mappings(self):
        if self.source is not None:
            try:
                zp_map = self.source._load_zp_mappings("./tests/resources/zfin/zp-mapping-test.txt")
                self.assertIsInstance(zp_map, dict, "_load_zp_mappings() didn't return dict!")
                self.assertTrue(len(zp_map) == 1, "_load_zp_mappings() didn't return exactly one thing!")
                self.assertDictEqual(zp_map, {'MONARCH:b311b709bf1cf852d770': {'zp_id': 'ZP:0002959', 'subterm1_id': 'ZFA:0009114', 'post_composed_relationship_id_1': 'BFO:0000050', 'superterm1_id': 'ZFA:0001056', 'quality_id': 'PATO:0001453', 'subterm2_id': 'GO:0005927', 'post_composed_relationship_id_2': 'BFO:0000050', 'superterm2_id': 'ZFA:0001056'}}, "_load_zp_mappings() didn't return what I expected!")
            except Exception as e:
                logger.error(e)
                
        return

    
if __name__ == '__main__':
    unittest.main()
