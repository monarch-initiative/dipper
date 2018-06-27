.. _best_practices:

Best Practices for writing source ingests as part of the Dipper system
======================================================================

Overview
--------

Each Source should extend the Source class have and override the fetch and parse functions.


Writing the fetcher
---------------------

*  This method is intended to fetch data from the remote locations (if
   it is newer than the local copy).
*  You can add a dictionary of files with a structure like:

.. code-block:: python

       files = {
           'somekey': {
               'file': 'filename.tsv',
               'url': 'http://example.org/filename.tsv'
           },
           ...
       }

Then call the `fetch <dipper.sources.Source.html#dipper.sources.Source.Source.fetch>`_ function to download each file.
If a the remote file has already been downloaded.  The fetch method checks the remote headers to see if it has been updated.
For sources not served over HTTP, this method may need to be overriden.
For example in `Bgee <dipper.sources.Bgee.html#dipper.sources.Bgee.Bgee.checkIfRemoteIsNewer>`_

Writing the parser
--------------------

Typically these are written by looping through the series of files that
were obtained by the fetch method. The goal is to process each file
minimally, adding classes and individuals as necessary, and adding
triples to the sources' graph.

There are certain conventions that we follow here:

1. Genes are a special case of genomic feature that are added as (OWL) Classes. But all
other genomic features are added as individuals of an owl class.

2. If a source references an external identifier then, assume that it has been
processed in a another Source script, and only add the identifier (but
not the label) to it within this Source's file. This will help prevent
label collisions related to slightly different versions of the source
data when integrating downstream.

3. You can instantiate a class or individual as many times as you want; they will get merged in the graph
and will only show up once in the resulting output.

Special configurations
----------------------

Add private configuration parameters into your private conf.json file.
Examples of items to put into the config include:

* database connection parameters (in the "dbauth" object)
* ftp login credentials
* api keys (in the "keys" object)

These are organized such that within any object (dbauth, keys, etc),
they are keyed again by the source's name.

Testing
-------

Integration tests
~~~~~~~~~~~~~~~~~

Integration tests can be executed by generating a file that contains a
subset of a source's data in the same format, and running it through the
source.parse() method, serializing the graph, and then testing this
file in some other piece of code or database.

You may see testing code within source classes, but these tests will be
deleted or refactored and moved to the test directory.

Unit tests
~~~~~~~~~~

Unit style tests can be achieved by mocking source classes (or specific
functions) and testing single functions. The
`test_graph_equality <dipper.utils.TestUtils.html#dipper.utils.TestUtils.TestUtils.test_graph_equality>`_
function can be used to test graph equality by supplying a string formatted as headless (no prefixes)
turtle and a graph object. Most dipper methods are not pure functions,
and rely on side effects to a graph object. Therefore it is best to
clear the graph object before any testing logic, eg:

.. code-block:: python

   from dipper.utils.TestUtils import TestUtils

   source.graph = RDFGraph(True)  # Reset graph
   test_util = TestUtils()
   source.run_some_function()
   expected_triples = """
       foaf:person1 foaf:knows foaf:person2 .
   """
   self.assertTrue(self.test_util.test_graph_equality(
       expected_triples, source.graph))
