.. _graphs:

Working with graphs
===================

The Dipper graph package provides two graph implementations, a RDFGraph which is
an extension of the RDFLib [1]_ Graph [2]_, and a StreamedGraph which prints triples
to standard out in the ntriple format.

RDFGraphs
---------
The RDFGraph class reads the curie_map.yaml file and converts strings formatted as curies
to RDFLib URIRefs.  Triples are added via the addTriple method, for example:

.. code-block:: python

   from dipper.graph.RDFGraph import RDFGraph

   graph = RDFGraph()
   graph.addTriple('foaf:John', 'foaf:knows', 'foaf:Joseph')


The graph can then be serialized in a variety of formats using RDFLib [3]_:

.. code-block:: python

   from dipper.graph.RDFGraph import RDFGraph

   graph = RDFGraph()
   graph.addTriple('foaf:John', 'foaf:knows', 'foaf:Joseph')
   print(graph.serialize(format='turtle').decode("utf-8"))

   # Or write to file
   graph.serialize(destination="/path/to/output.ttl", format='turtle')


Prints:

.. code-block:: text

   @prefix OBO: <http://purl.obolibrary.org/obo/> .
   @prefix foaf: <http://xmlns.com/foaf/0.1/> .
   @prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
   @prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
   @prefix xml: <http://www.w3.org/XML/1998/namespace> .
   @prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

   foaf:John foaf:knows foaf:Joseph .

When an object is a literal, set the object_is_literal param to True

.. code-block:: python

   from dipper.graph.RDFGraph import RDFGraph

   graph = RDFGraph()
   graph.addTriple('foaf:John', 'rdfs:label', 'John', object_is_literal=True)

Literal types can also be passed into the method:

.. code-block:: python

   from dipper.graph.RDFGraph import RDFGraph

   graph = RDFGraph()
   graph.addTriple(
       'foaf:John', 'foaf:age', 12,
        object_is_literal=True, literal_type="xsd:integer"
   )

StreamedGraphs
-------------

StreamedGraphs print triples as they are processed by the addTriple method.  This is useful for
large sources where.  The output should be sorted and uniquified as there is no checking for
duplicate triples.  For example:

.. code-block:: python

   from dipper.graph.StreamedGraph import StreamedGraph

   graph = StreamedGraph()
   graph.addTriple('foaf:John', 'foaf:knows', 'foaf:Joseph')

Prints:

.. code-block:: text

   <http://xmlns.com/foaf/0.1/John> <http://xmlns.com/foaf/0.1/knows> <http://xmlns.com/foaf/0.1/Joseph> .


References
----------

.. [1] RDFLib: `<http://rdflib.readthedocs.io/en/stable/>`_
.. [2] RDFLib Graphs: `<https://rdflib.readthedocs.io/en/stable/apidocs/rdflib.html#graph-module>`_
.. [3] RDFLib Serializing: `<http://rdflib.readthedocs.io/en/stable/apidocs/rdflib.html#rdflib.graph.Graph.serialize>`_


