.. _schemas:

Schemas
=======

Although RDF is inherently schemaless, we aim to construct consistent models across sources.  This allows us to
build source agnostic queries and bridge data across sources.

The dipper schemas are documented as directed graphs.  Examples can be found in the
`ingest artifacts repo <https://github.com/monarch-initiative/ingest-artifacts/tree/master/sources>`_.

Some ontologies contain documentation on how to describe data using the classes and properties defined in the ontology:

* `The Scientific Evidence and Provenance Information Ontology (SEPIO) <https://github.com/monarch-initiative/SEPIO-ontology/tree/master/docs>`_
* `The Genotype Ontology (GENO) <https://github.com/monarch-initiative/GENO-ontology/tree/develop/docs>`_

While not yet implemented, in the future we plan on defining our schemas
and constraints using the `BioLink model specification <https://biolink.github.io/biolink-model/>`_.

The cypher queries that we use to cache inferred and direct relationships between entities
is `stored in GitHub <https://github.com/monarch-initiative/monarch-cypher-queries/tree/master/src/main/cypher/golr-loader>`_.

