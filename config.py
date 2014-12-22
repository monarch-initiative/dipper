__author__ = 'nicole'

import json
import os.path

#read configuration file
conf = None

#load the configuration file, if it exists.
#it isn't required, but may be for some sources
if (os.path.exists("conf.json")):
    with open("conf.json") as json_file:
        conf = json.load(json_file)
        print("Finished loading config:",conf)

def get_config():
    return conf