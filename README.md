* This requires [Python 3](https://www.python.org/downloads/) or higher
* Be sure to create the following directories to hold incoming and created files:
    * ./out
    *  ./raw
* The basic verification step requires (this can be turned off)
[owltools](https://code.google.com/p/owltools/wiki/InstallOWLTools) be available on your path

* Required python packages:
    * [psycopg2](http://initd.org/psycopg/)
    * [rdflib](https://code.google.com/p/rdflib/)
    * isodate
    * roman

* right now this runs completely through the __init__() method.  you can select the sources to process
by commenting them on/off.

* the general strategy here is that there is one class per data source.  we define the files to be fetched,
any file scrubbing, and then the parsing methods.  as the files are parsed, they are turned into triples and
loaded into a graph.  this graph is then typically dumped into triples (currently coupled to the parsers,
but will be moved into a separate method in the future to be controlled from main.

* The following sources have been partially mapped:
    * Human Phenotype Ontology Annotations (Disease-Phenotype associations)
    * BioGrid (Gene-gene interactions)
    * ZFIN (Genotypes (with their parts), and Genotype-Phenotype associations)
    * OMIM (Diseases, labels, and descriptions)
