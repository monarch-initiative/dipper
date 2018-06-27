.. Dipper documentation master file, created by
   sphinx-quickstart on Fri Feb 27 15:32:06 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Dipper's documentation!
==================================

Dipper is a Python package to generate RDF triples from common scientific resources.
Dipper includes subpackages and modules to create graphical models of this data, including:

* Models package for generating common sets of triples, including common OWL axioms, complex genotypes, associations, evidence and provenance models.

* Graph package for building graphs with RDFLib or streaming n-triples

* Source package containing fetchers and parsers that interface with remote databases and web services

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   notebooks
   best_practices
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
