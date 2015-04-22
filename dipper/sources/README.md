#Best Practices for writing Source parsers as part of the Dipper system

##Overview
Each Source should typically have the following basic things defined:

  fetch()

  parse()


##Writing the fetcher()
* this method is intended to go and fetch the data from the remote locations.  
* if you can add a dictionary of files with a structure like 

    ```python
    files = {
        key : {'file' : local.file.name,
               'url' : 'http://remote.file.location'
        }, ... }
    ```
then you can simply call the built-in file fetching method to go get all of them.


##Writing the parser()

Typically these are written by looping through the series of files that were obtained in the fetching function.
The goal is to process each file minimally, adding classes and individuals as necessary, and adding triples
to the sources' graph.

There are certain conventions that we follow here:
  1. Genes are a special case of genomic feature that are added as Classes.  But all other genomic features
  are added as individuals
  2. If a source references an external identifier, assume that it has been processed in a another Source script,
  and only add the identifier (but not the label) to it within this Source's file.  This will help prevent label
  collisions related to slightly different versions of the source data when integrating downstream.
  3. You can instantiate a class or individual as many times as you want; they will get merged in the graph and
  will only show up once in the resulting output.


## Special configurations

Add private configuration parameters into your private config.json file.  Examples of items to put into the config
include:
* database connection parameters (in the "dbauth" object)
* ftp login credentials
* api keys (in the "keys" object)
  
These are organized such that within any object (dbauth, keys, etc), they are keyed again by the source's name.
  
##Testing

###Add Test IDs to the config.
Determine what IDs are suitable to make a representative testing subset.  If possible, add them to the config.
The config is meant to have a set of cross-cutting identifiers that occur in multiple sources.  This will permit
cross-source integration testing in downstream platforms like SciGraph.  


###WRITE TEST(s)
Any Source should have a test suite written around it.  You should write a standalone testing suite in the "tests"
directory, following the convention of "test_<source_name>.py", and subclass it from ```GeneralSourceTestCase```.  

The main ```dipper.py``` wrapper will register all tests within the source and run them automatically.


