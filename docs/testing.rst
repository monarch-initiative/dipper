.. _testing:

Testing ingests
===============

Unit tests
----------

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

Integration tests
-----------------

Integration tests can be executed by generating a file that contains a
subset of a source's data in the same format, and running it through the
source.parse() method, serializing the graph, and then testing this
file in some other piece of code or database.

You may see testing code within source classes, but these tests will be
deleted or refactored and moved to the test directory.