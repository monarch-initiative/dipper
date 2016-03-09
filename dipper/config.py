import json
import os.path
import logging

__author__ = 'nicole'

logger = logging.getLogger(__name__)

#read configuration file
conf = {}

'''
    Load the configuration file 'conf.json', if it exists.
    it isn't always required, but may be for some sources.
'''
 
if os.path.exists(os.path.join(os.path.dirname(__file__), 'conf.json')):
    with open(
        os.path.join(os.path.dirname(__file__),
                     'conf.json')) as json_file:
        conf = json.load(json_file)
        logger.debug("Finished loading config")
else:
    logger.warning("'conf.json' not found in '%s'", os.path.dirname(__file__))
    logger.warning("Sources that depend on 'conf.json' will fail")

def get_config():
    return conf
