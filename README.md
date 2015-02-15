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
    * [python-docx](https://github.com/python-openxml/python-docx)
    * pyyaml

* The OMIM source requires the 'compress' and 'uncompress' system commands to be available, for LZW decompression.  
(This may be a problem for windows users.) 

* **To run the software**:
    * you can run the code by supplying a list of one or more sources on the command line.  some examples:
      dipper --sources omim,ncbigene
    * furthermore, you can run this in a testing mode, by supplying a limit for testing.  this will only test the first N number of rows or data elements
      dipper --sources omim --limit 100
    * other commandline parameters are explained in the code
     
* The general strategy here is that there is one class per data source.  We define the files to be fetched,
any file scrubbing, and then the parsing methods.  As the files are parsed, triples are loaded into a graph.  
This graph is then typically dumped into triples in turtle format.

* The following sources have been partially mapped:
    * Human Phenotype Ontology Annotations (Disease-Phenotype abnormal associations)
    * BioGrid (Gene-gene interactions)
    * ZFIN (Genotypes (with their parts), and Genotype-Phenotype associations)
    * OMIM (Diseases, labels, and descriptions; disease-gene associations)
    * Panther (orthologs for ~20 reference genomes)
    * UCSCBands (taking the chromosomal banding patterns and converting them to a nested and genomically located structure)
    * GeneReviews (mapping their ids to OMIM)
    * NCBIGene (Genes, equivalent ids, taxa, and their chr locations)
