import logging
import os
from enum import Enum
from pathlib import Path

import yaml

"""
This is a class to support categorizing everything we ingest using biolink 
categories. It reads dipper/biolink_vocabulary.yaml and populates an Enum that
is used to assign biolink categories to entities we ingest.

If you need to add a biolink category, do so in dipper/biolink_vocabulary.yaml
"""

LOG = logging.getLogger(__name__)


class BioLinkVocabulary:
    bl_file_with_path = os.path.join(Path(__file__).parents[1],
                                     'biolink_vocabulary.yaml')
    if os.path.exists(bl_file_with_path):
        with open(os.path.join(bl_file_with_path)) as yaml_file:
            bl_vocab = yaml.safe_load(yaml_file)
            LOG.debug("Loaded biolink vocabulary: %s", bl_vocab)
    else:
        LOG.debug("Cannot find biolink vocab yaml file: %s", bl_file_with_path)
    terms = Enum('DynamicEnum', bl_vocab)

