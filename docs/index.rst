Welcome to Dipper's documentation!
==================================

Dipper is a Python package to generate RDF triples from common scientific resources.
Dipper includes subpackages and modules to create graphical models of this data, including:

* Models package for generating common sets of triples, including common OWL axioms, complex genotypes, associations, evidence and provenance models.

* Graph package for building graphs with RDFLib or streaming n-triples

* Source package containing fetchers and parsers that interface with remote databases and web services

Getting started
---------------

Installing, running, and the basics

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   notebooks
   downloads
   status
   applications

Deeper into Dipper
--------------

A look into the structure of the codebase and how to write ingests

.. toctree::
   :maxdepth: 1

   graphs
   models
   writing_ingests
   testing
   configuration
   schemas

For developers
--------------

.. toctree::
   :maxdepth: 1

   modules
   sources


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
