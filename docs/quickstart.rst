.. _quickstart:

Getting started with Dipper
===========

This guide assumes you have already installed dipper.  If not, then follow the steps in the
:ref:`installation` section.

Command line
------------

You can run the code by supplying a list of one or more sources on the command line. some examples:

::

   dipper-etl.py --sources impc,hpoa

Furthermore, you can check things out by supplying a limit. this will only process the first N number of rows or data elements:

::

   dipper-etl.py --sources hpoa --limit 100

Other command line parameters are explained if you request help:

::

   dipper-etl.py --help

Notebooks
---------

We provide `Jupyter Notebooks <http://nbviewer.jupyter.org/github/monarch-initiative/dipper/tree/master/docs/notebooks/>`_
to illustrate the functionality of the python library. These can also
be used interactively.

See the :ref:`notebooks` section for more details.

Building models
------

This code example shows some of the basics of building RDF graphs using the models API:

.. code-block:: python

   import pandas as pd
   from dipper.graph.RDFGraph import RDFGraph
   from dipper.models.Model import Model

   columns = ['variant', 'variant_label', 'variant_type',
              'phenotype','relation', 'source', 'evidence', 'dbxref']

   data =  [
        ['ClinVarVariant:254143', 'C326F', 'SO:0000694',
         'HP:0000504','RO:0002200', 'PMID:12503095', 'ECO:0000220',
         'dbSNP:886037891']
   ]

   # Initialize graph and model
   graph = RDFGraph()
   model = Model(graph)

   # Read file
   dataframe = pd.DataFrame(data=data, columns=columns)

   for index, row in dataframe.iterrows():
      # Add the triple ClinVarVariant:254143 RO:0002200 HP:0000504
      # RO:0002200 is the has_phenotype relation
      # HP:0000748 is the phenotype 'Inappropriate laughter'
      model.addTriple(row['variant'], row['relation'], row['phenotype'])

      # The addLabel method adds a label using the rdfs:label relation
      model.addLabel(row['variant'], row['variant_label'])

      # addType makes the variant an individual of a class,
      # in this case SO:0000694 'SNP'
      model.addType(row['variant'], row['variant_type'])

      # addXref uses the relation OIO:hasDbXref
      model.addXref(row['variant'], row['dbxref'])

      # Serialize the graph as turtle
      print(graph.serialize(format='turtle').decode("utf-8"))

For more information see the :ref:`models` section.