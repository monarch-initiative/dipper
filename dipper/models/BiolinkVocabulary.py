import logging
import os
from pathlib import Path

import yaml

"""
This is a class to support categorizing everything we ingest using biolink 
categories. It reads resources/biolink_vocabulary.yaml, populates an Enum with the
contents of this yaml file, such that each term is a name, and the value for each
term is a curie for the term, formed by prepended the term with `curie_prefix`. 
The CURIEs are used to assign biolink categories to entities we ingest.

If you need to add a biolink category, do so in dipper/biolink_vocabulary.yaml
"""

LOG = logging.getLogger(__name__)


class BioLinkVocabulary:
    bl_file_with_path = os.path.join(Path(__file__).parents[2],
                                     'resources/biolink_vocabulary.yaml')
    if os.path.exists(bl_file_with_path):
        with open(os.path.join(bl_file_with_path)) as yaml_file:
            bl_vocab = yaml.safe_load(yaml_file)
            LOG.debug("Loaded biolink vocabulary: %s", bl_vocab)
            terms = {} # keys are terms, values are bl CURIES
            for key in bl_vocab["terms"]:
                terms[key] = bl_vocab["curie_prefix"] + ":" + key
    else:
        LOG.debug("Cannot find biolink vocab yaml file: %s", bl_file_with_path)
