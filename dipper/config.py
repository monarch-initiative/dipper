import yaml
import os.path
import logging

__author__ = 'nicole'

logger = logging.getLogger(__name__)

# default configuration
conf = {
  "keys": {
    "omim": ''
  }
}

'''
    Load the configuration file 'conf.json', if it exists.
    it isn't always required, but may be for some sources.
    conf.json may contain sensitive info and should not live in a public repo
'''

if os.path.exists(os.path.join(os.path.dirname(__file__), 'conf.yaml')):
    with open(
        os.path.join(
            os.path.dirname(__file__), 'conf.yaml')) as yaml_file:
        conf = yaml.safe_load(yaml_file)
        logger.debug("Finished loading dipper/conf.yaml")
else:
    logger.warning("'dipper/conf.yaml' not found in '%s'", os.path.dirname(__file__))
    logger.warning("Sources that depend on 'conf.yaml' will fail")


def get_config():
    return conf
