* The dipper main wraps all of the source parsers, enabling users to specify one or more sources to process. 
The general strategy is that there is one class per data source.  We define the files to be fetched,
any file scrubbing, and then the parsing methods.  As the files are parsed, triples are loaded into a graph.  
This graph is then typically dumped into triples in turtle format.  A subsetted file *_test.ttl will also be dumped.

* Dipper requires [Python 3](https://www.python.org/downloads/) or higher
* One of the unit tests requires
[owltools](https://code.google.com/p/owltools/wiki/InstallOWLTools) be available on your path.  You could modify
the code to skip this, if necessary

* Required external python packages:
    * [psycopg2](http://initd.org/psycopg/)
    * [rdflib](https://code.google.com/p/rdflib/)
    * isodate
    * roman
    * [python-docx](https://github.com/python-openxml/python-docx)
    * pyyaml
    * pysftp
    * [biopython](https://github.com/biopython/biopython)

* The OMIM source requires the 'compress' and 'uncompress' system commands to be available, for LZW decompression.  
(This may be a problem for windows users.) 

* **To run the software**:
    * you can run the code by supplying a list of one or more sources on the command line.  some examples:
      ```dipper --sources omim,ncbigene```
    * furthermore, you can run this in a testing mode, by supplying a limit for testing.  this will only test the 
    first N number of rows or data elements
      ```dipper --sources omim --limit 100```
    * you can run the stand-alone tests in ```tests/test_*``` to generate subsets of the data
    * other commandline parameters are explained in the code
    
* You can also write your own dipper packages outside of this project, using the framework we've set up here.  Simply
import Dipper as a python package and write your own wrapper.
    Dipper can be installed as an external python package with pip3

        ```pip3 install git+git://github.com/monarch-initiative/dipper.git```
        
    or clone the repository and run:
    
        ```pip3 install -e /path/to/git/dipper```

* The following sources have been mapped:
    * Human Phenotype Ontology Annotations (Disease-Phenotype abnormal associations)
    * BioGrid (Gene-gene interactions)
    * ZFIN (Genotypes (with their parts), and Genotype-Phenotype associations)
    * OMIM (Diseases, labels, and descriptions; disease-gene associations)
    * Panther (orthologs for ~20 reference genomes)
    * UCSCBands (taking the chromosomal banding patterns and converting them to a nested and genomically located structure)
    * GeneReviews (mapping their ids to OMIM)
    * NCBIGene (Genes, equivalent ids, taxa, and their chr locations)
    * IMPC (Genotypes (with their parts) and Genotype-Phenotype associations)
    * MGI (Genotypes (with their parts) and Genotype-Phenotype associations)
    * CTD (Chemicals and Genes)
    * ClinVar (Human Disease and Variant associations)
    * Coriell (Cell lines, many as models of human disease)
    * EOM - Elements of Morphology (links to images of human phenotypes)
    
* You could contribute a Source parser to our suite!  Please see our [best-practies documentation](sources/README.md) for details on
writing new Source parsers using Dipper code.  