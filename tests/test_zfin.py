#!/usr/bin/env python3

import unittest
import logging
from rdflib.namespace import URIRef
from tests.test_source import SourceTestCase
from dipper.sources.ZFIN import ZFIN

logging.basicConfig(level=logging.WARNING)
LOGGER = logging.getLogger(__name__)


class ZFINTestCase(SourceTestCase):

    def setUp(self):
        self.source = ZFIN('rdf_graph', True)
        self.source.settestonly(True)
        self._setDirToSource()

    def tearDown(self):
        self.source = None

    def test_mapping_of_phenotypes_to_zp_ids(self):
        """
        test that code correctly uses zp_map to map phenotypes to zp ids
        :return:

        """
        mapping_file = "./tests/resources/zfin/zp-mapping-test-map.txt"
        pheno_file = "./tests/resources/zfin/zp-mapping-test-phenotype.txt"

        self.source.zp_map = self.source._load_zp_mappings(mapping_file)
        pheno_dat = open(pheno_file).read().split('\t')

        (fish_num, fish_name, start_stage_id, start_stage_name, end_stage_id,
         end_stage_name, subterm1_id, subterm1_name, postcomp1_rel_id,
         postcomp1_rel_name, superterm1_id, superterm1_name, quality_id,
         quality_name, modifier, subterm2_id, subterm2_name, postcomp2_rel_id,
         postcomp2_rel_name, superterm2_id, superterm2_name, pub_id, env_id) = \
            pheno_dat

        self.assertEqual(
            self.source._map_octuple_to_phenotype(subterm1_id, postcomp1_rel_id,
                                                  superterm1_id, quality_id,
                                                  subterm2_id, postcomp2_rel_id,
                                                  superterm2_id, "abnormal"),
            'ZP:0022140')

    def test_load_zp_mappings(self):
        """
        test correct loading of zp mappings file and construction of zp_map

        """
        if self.source is not None:
            try:
                zp_map = self.source._load_zp_mappings(
                    "./tests/resources/zfin/zp-mapping-test.txt")
                self.assertIsInstance(zp_map, dict,
                                      "_load_zp_mappings() didn't return dict!")
                self.assertTrue(len(zp_map) == 1,
                                "_load_zp_mappings() didn't return exactly one thing!")
                self.assertDictEqual(zp_map,
                                     {'MONARCH:b308a8f1c67793a56d16':
                                          {'post_composed_relationship_id_1':
                                               'BFO:0000050',
                                           'post_composed_relationship_id_2':
                                               'BFO:0000050',
                                           'quality_id': 'PATO:0001453',
                                           'subterm1_id': 'ZFA:0009114',
                                           'subterm2_id': 'GO:0005927',
                                           'superterm1_id': 'ZFA:0001056',
                                           'superterm2_id': 'ZFA:0001056',
                                           'zp_id': 'ZP:0002959',
                                           'modifier': 'PATO:0000460'}},
                                     "_load_zp_mappings() " +
                                     "didn't return what I expected!")
            except Exception as t_except:
                LOGGER.error(t_except)

    def test_make_zpkey(self):
        """
        test that _make_zpkey returns correct id

        """
        if self.source is not None:
            try:
                dummy_args = list(map(str, list(range(1, 9)))) # 1 - 8 as strings
                expected_key = self.source.make_id("_".join(dummy_args))
                self.assertEqual(self.source._make_zpkey(*dummy_args), expected_key)
                self.assertEqual(self.source._make_zpkey(['0'] * 8),
                                 self.source._make_zpkey([''] * 8),
                                 "_make_zpkey() doesn't seem to be replacing empty " +
                                 "strings with zeros before making key," +
                                 "this might cause zp_map lookup issues")

            except Exception as t_except:
                LOGGER.error(t_except)

    def test_genotype_labels(self):
        """
        test that genotype label is set correctly after parse()

        """
        if self.source is not None:
            test_resource_dir = "../../tests/resources/zfin/"
            self.source.files['fish_components']['file'] = test_resource_dir + \
                "genotype-label-test-fish_components_fish.txt"
            self.source.files['backgrounds']['file'] = test_resource_dir + \
                "genotype-label-test-genotype_backgrounds.txt"
            self.source.files['geno']['file'] = test_resource_dir + \
                "genotype-label-test-genotype_features.txt"
            self.source.parse()
            try:
                this_iri = URIRef("http://zfin.org/ZDB-GENO-070228-3")
                expect_label: str = "shha<sup>tbx392/tbx392</sup> [AB]"
                assert(str(self.source.testgraph.label(this_iri, None)) == expect_label)
            except Exception as t_except:
                LOGGER.error(t_except)


if __name__ == '__main__':
    unittest.main()
