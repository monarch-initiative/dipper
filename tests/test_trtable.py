#!/usr/bin/env python3

import unittest
import os
import yaml


class TranslationTestCase(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return
    
    def testIfTableIsAMap(self):

        from yaml.constructor import ConstructorError

        try:
            from yaml import CLoader as Loader
        except ImportError:
            from yaml import Loader
        
        # Credit https://gist.github.com/pypt/94d747fe5180851196eb
        def no_duplicates_constructor(loader, node, deep=False):
            """Check for duplicate keys."""

            mapping = {}
            for key_node, value_node in node.value:
                key = loader.construct_object(key_node, deep=deep)
                value = loader.construct_object(value_node, deep=deep)
                if key in mapping:
                    raise ConstructorError("while constructing a mapping", node.start_mark,
                                           "found duplicate key (%s)" % key, key_node.start_mark)
                mapping[key] = value

            return loader.construct_mapping(node, deep)

        yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, no_duplicates_constructor)

        file_path = '../translationtable/global_terms.yaml'
        if os.path.exists(os.path.join(os.path.dirname(__file__), file_path)):
            tt_file = open(os.path.join(os.path.dirname(__file__), file_path), 'r')
            try:
                translation_table = yaml.load(tt_file)
            except yaml.constructor.ConstructorError as e:
                tt_file.close()
                self.assertTrue(False)
                print(e)
            tt_file.close()

        
        
 
    def testIfTableIsBiMap(self):
        file_path = '../translationtable/global_terms.yaml'
        if os.path.exists(os.path.join(os.path.dirname(__file__), file_path)):
            tt_file = open(os.path.join(os.path.dirname(__file__), file_path), 'r')
            translation_table = yaml.load(tt_file)
            tt_file.close()

        temp_dict = {}
        passed = True
        failed_list = []
        for k,v in translation_table.items():
            if v not in temp_dict:
                temp_dict[v] = k
            else:
                passed = False
                failed_list.append(v)

        if not passed:
            print("Duplicate values in yaml: {}".format(failed_list))        
                        
        self.assertTrue(passed)



