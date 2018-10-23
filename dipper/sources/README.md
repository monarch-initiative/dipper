# Best Practices for writing source ingests as part of the Dipper system

## Overview
Each Source should typically have the following basic things defined:

  fetch()

  parse()


## Writing the fetcher()
* This method is intended to fetch data from the remote locations
(if it is newer than the local copy).  
* You can add a dictionary of files with a structure like:

    ```python
    files = {
        'somekey': {  
            'file': 'filename.omg',  
            'url': 'http://example.org/filename.omg'  
        }, 
        ...  
    }  
    ```
then call the built-in file fetching method in Source.py ```get_files()```
to get them all.


## Writing the parser()

Typically these are written by looping through the series of files that were
obtained by the  method fetch().
The goal is to process each file minimally,
adding classes and individuals as necessary,
and adding triples to the sources' graph.  

There are certain conventions that we follow here:
  1. Genes are a special case of genomic feature that are added as (OWL) Classes.
    But all other genomic features are added as individuals of an owl class
  2. If a source references an external identifier then,
     assume that it has been processed in a another Source script,
     and only add the identifier (but not the label) to it within this Source's file.
    This will help prevent label collisions related
    to slightly different versions of the source data when integrating downstream.
  3. You can instantiate a class or individual as many times as you want;
  they will get merged in the graph and will only show up once in the resulting output.


## Special configurations

Add private configuration parameters into your private conf.yaml file.
Examples of items to put into the config include:  
* database connection parameters (in the "dbauth" object)  
* account user names
* paswords
  * ftp login credentials  
  * api keys (in the "keys" object)
* ssh key pairs    
  
These are organized such that within any object (dbauth, users, keys, etc),
they are keyed again by the source's name.
  
## Testing

We employ two types of testing, but are always open to new and improved
methods of testing.

### Integration tests
Integration tests can be executed by generating a file that contains a
subset of a sources data in the same format, and running it through the
source.parse() function, serializing the graph, and then testing this file
in some other piece of code or database.

You may see testing code within source classes, but these tests will be
deleted or refactored and moved to the test directory.


### Unit Tests
Unit style tests can be achieved by mocking source classes
(or specific functions) and testing single functions. The [TestUtils](../utils/TestUtils.py)
class contains a function to test graph equality by supplying a string formatted as headless (no prefixes) turtle
and a graph object.  Most dipper methods are not pure functions, and rely on side effects to a graph object.
Therefore it is best to clear the graph object before any testing logic, eg:

   ```python
    source.graph = RDFGraph(True)  # Reset graph
    source.run_some_function()
    expected_triples = """
        foaf:person1 foaf:knows foaf:person2 .
    """
    self.assertTrue(self.test_util.test_graph_equality(
            expected_triples, source.graph))
   ```
