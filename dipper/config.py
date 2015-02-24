__author__ = 'nicole'

import json
import os.path

#read configuration file
conf = None

#load the configuration file, if it exists.
#it isn't required, but may be for some sources
if os.path.exists(os.path.join(os.path.dirname(__file__), 'conf.json')):
    with open(os.path.join(os.path.dirname(__file__),
                           'conf.json')) as json_file:
        conf = json.load(json_file)
        print("INFO: Finished loading config:")

def get_config():
    return conf