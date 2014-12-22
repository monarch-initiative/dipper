#first-pass with dipper
#this will eventually control the processing of data sources
__author__ = 'nlw'

import config
from sources.HPOAnnotations import HPOAnnotations
from sources.ZFIN import ZFIN
from sources.OMIM import OMIM


source_to_class_map={
#    'hpoa' : HPOAnnotations,
#    'zfin' : ZFIN,
    'omim' : OMIM
}

#load configuration parameters
#for example, keys


#TODO subset of sources will eventually be configurable on the commandline
#iterate through all the sources
for source in source_to_class_map.keys():
    mysource = source_to_class_map[source]()
    mysource.fetch()
    mysource.parse(100)
    status = mysource.verify()
    if status is not True:
        print('ERROR: Source',source,'did not pass verification tests.')
    print('***** Finished with',source,'*****')
print("All done.")

#TODO command-line args:
# *force re-download
# *specify the source
# *parse only without writing
# *parse only X lines of original file

###########################



