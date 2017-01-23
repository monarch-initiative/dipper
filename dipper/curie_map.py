'''
    Acroname central

    Load the curie mapping file 'curie_map.yaml', if it exists.
    it isn't required, but is necessary for most resources

'''
import os.path
import logging
import yaml

__author__ = 'nicole'

logger = logging.getLogger(__name__)

# read configuration file
curie_map = None

if os.path.exists(os.path.join(os.path.dirname(__file__), 'curie_map.yaml')):
    with open(os.path.join(os.path.dirname(__file__), 'curie_map.yaml')) \
            as yaml_file:
        curie_map = yaml.load(yaml_file)
        logger.debug("Finished loading curie maps: %s", curie_map)
else:
    logger.debug(
        "Cannot find 'curie_map.yaml' in  %s",
        os.path.dirname(__file__))


def get():
    return curie_map


def get_base():
    return curie_map['']