import json
import os.path
import logging

__author__ = 'nicole'

logger = logging.getLogger(__name__)

# read configuration file
conf = {}

'''
    Load the configuration file 'conf.json', if it exists.
    it isn't always required, but may be for some sources.
    conf.json may contain sensitive info and should not live in a public repo
'''

if os.path.exists(os.path.join(os.path.dirname(__file__), 'conf.json')):
    with open(
        os.path.join(
            os.path.dirname(__file__), 'conf.json')) as json_file:
        conf = json.load(json_file)
        logger.debug("Finished loading dipper/config.json")
else:
    logger.warning("'dipper/conf.json' not found in '%s'", os.path.dirname(__file__))
    logger.warning("Sources that depend on 'conf.json' will fail")


def get_config():
    return conf
