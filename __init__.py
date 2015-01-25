#first-pass with dipper
#this will eventually control the processing of data sources
__author__ = 'nlw'

import config
from sources.HPOAnnotations import HPOAnnotations
from sources.ZFIN import ZFIN
from sources.OMIM import OMIM
from sources.BioGrid import BioGrid
from sources.MGI import MGI
from sources.IMPC import IMPC
from sources.Panther import Panther
from sources.NCBIGene import NCBIGene
from sources.UCSCBands import UCSCBands

source_to_class_map={
#    'hpoa' : HPOAnnotations, # ~3 min
#   'zfin' : ZFIN,
#    'omim' : OMIM,  #full file takes ~15 min, due to required throttling
#    'biogrid' : BioGrid,  #interactions file takes <10 minutes
#    'mgi' : MGI,
#    'impc' : IMPC,
#    'panther' : Panther,  #this takes a very long time, ~1hr to map 7 species-worth of associations
#    'ncbigene' : NCBIGene,  #takes about 4 minutes to process 2 species
    'bands' : UCSCBands
}

#load configuration parameters
#for example, keys


#TODO subset of sources will eventually be configurable on the commandline
#iterate through all the sources
for source in source_to_class_map.keys():
    print()
    print("*******",source,"*******")
    mysource = source_to_class_map[source]()
    mysource.fetch()
    mysource.parse()
    mysource.write(format='turtle')
    #status = mysource.verify()
#    if status is not True:
#        print('ERROR: Source',source,'did not pass verification tests.')
#    print('***** Finished with',source,'*****')

print("All done.")

#TODO command-line args:
# *force re-download
# *specify the source
# *parse only without writing
# *parse only X lines of original file
# *quiet mode (with proper logging methods)

###########################



