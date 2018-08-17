.. _models:

Working with the model API
===========================
The model package provides classes for building common sets of triples based on our modeling patterns.

For an example see the notebook on this topic: `Building graphs with the model API <http://nbviewer.jupyter.org/github/monarch-initiative/dipper/blob/master/docs/notebooks/model-api-tutorial.ipynb>`_


Basics
------
The model class provides methods for building common RDF and OWL statements

For a list of methods, `see the API docs <dipper.models.Model.html>`_.


Building associations
---------------------

We use the RDF Reification [1]_ pattern to create ternary statements, for example, adding frequency
data to phenotype to disease associations.
We utilize the Open Biomedical Association ontology [2]_ to reify statements, and the SEPIO ontology to
add evidence and provenance.

For a list of classes and methods, `see the API docs <dipper.models.assoc.html>`_.

Building genotypes
------------------

We use the GENO ontology [4]_ to build complex genotypes and their parts.

For a list of methods, `see the API docs <dipper.models.Genotype.html>`_.

GENO docs: `The Genotype Ontology (GENO) <https://github.com/monarch-initiative/GENO-ontology/tree/develop/docs>`_

Building complex evidence and provenance graphs
-----------------------------------------------
We use the SEPIO ontology to build complex evidence and provenance.  For an example see the IMPC source ingest.

For a list of methods, see the API docs for `evidence <dipper.models.Evidence.html>`_ and `provenance <dipper.models.Provenance.html>`_.

SEPIO docs: `The Scientific Evidence and Provenance Information Ontology (SEPIO) <https://github.com/monarch-initiative/SEPIO-ontology/tree/master/docs>`_


References
----------

.. [1] RDF Reification: `<https://www.w3.org/TR/rdf-primer/#reification>`_
.. [2] OBAN: `<https://github.com/EBISPOT/OBAN>`_
.. [3] SEPIO: `<https://github.com/monarch-initiative/SEPIO-ontology>`_
.. [4] GENO: `<https://github.com/monarch-initiative/GENO-ontology>`_
