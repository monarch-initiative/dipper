#first-pass with dipper
#this will eventually control the processing of data sources
__author__ = 'nlw'

from sources.HPOAnnotations import HPOAnnotations
from sources.ZFIN import ZFIN


#the following is the registry of data to process
all_sources = ('hpoa','other')

#TODO the subset of sources will either be all, or eventually configurable on the commandline
sources=all_sources
#iterate through all the sources

mysource = HPOAnnotations()
#mysource.whoami()
mysource.fetch()
mysource.parse(1000)

mysource = ZFIN()
mysource.fetch()
mysource.parse(10)

print("All done.")

#TODO command-line args:
# *force re-download
# *specify the source
# *parse only without writing
# *parse only X lines of original file

###########################



