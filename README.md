[![Build Status](https://travis-ci.org/monarch-initiative/dipper.svg?branch=master)](https://travis-ci.org/monarch-initiative/dipper)

# DIPPER
Dipper is a pure Python package to generate RDF triples from common scientific resources.
Dipper includes subpackages and modules to create graphical models of this data, including:

* Models for ternary associations and complex partonomies
* Graph utilities to create common relationships using a variety of ontologies
* Fetchers and Parsers that interface with remote databases and web services

## Beware: This code and resulting TTL is still in active development, and is subject to change.
The code in this repository is still under active development while we are refining the models for variant, 
genotype, and genotype-to-phenotype association.  This means that the graph patterns we are employing 
between the different entities is in flux.  **We will indicate here when the modeling has generally stabilized.**

* The dipper main wraps all of the source parsers, enabling users to specify one or more sources to process. 
The general strategy is that there is one class per data source.  We define the files to be fetched,
any file scrubbing, and then the parsing methods.  As the files are parsed, triples are loaded into an in-memory graph.  
This graph is then typically dumped into triples in turtle format.  For testing purposes,
 a subset of the graph is also dumped to *_test.ttl.
* Data generated from this pipeline can be used in a variety of ways downstream.  We recommend
loading the data into a graph database that is optimized for use with ontologies, such as 
[SciGraph](https://github.com/SciGraph).  Smaller .ttl files can be loaded into an ontology editor 
like [Protege](http://protege.stanford.edu/).

## Requirements
* [Python 3](https://www.python.org/downloads/) or higher
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

* Some of the parsers require login and/or connection credentials with the remote system.  In those cases
 you will need to add the credentials to a conf.json file.  Please see individual parsers for details.   

## Running Dipper:
* you can run the code by supplying a list of one or more sources on the command line.  some examples:
      ```dipper --sources omim,ncbigene```
* furthermore, you can check things out by supplying a limit.  this will only process the
first N number of rows or data elements
    ```dipper --sources hpoa --limit 100```
* you can also run the stand-alone tests in ```tests/test_*``` to generate subsets of the data and run unittests
* other commandline parameters are explained if you request help:
    ```./dipper.py --help```

## Installing Dipper as an external python package:
You can also write your own dipper packages outside of this project, using the framework we've set up here.  Simply
import Dipper as a python package, write your own wrapper, and add your own source parsers.
* as an external python package with pip3
    ```pip3 install git+git://github.com/monarch-initiative/dipper.git```
* or clone the repository and run:
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
    
* Don't see a parser you want?  Feel free to request a new one, or you could contribute a Source parser to our suite!  
Please see our [best-practies documentation](sources/README.md) for details on writing new Source parsers 
using Dipper code, and make a Pull request.  

## License (Under consideration)

## About this project
The DIPper data pipeline was born out of the need for a uniform representation of human and model organism
genotype-to-phenotype data, and an easy Extract-Transform-Load (ETL) pipeline to process it all.  
It became too cumbersome to first get all of these data into a single-schema traditional SQL database, 
then transform it into a graph representation.  So, we decided to go straight from each source into triples that 
are semantically captured, using standard modeling patterns.  
Furthermore, we wanted to provide the bioinformatics community with a set of scripts to help anyone 
get started transforming these standard data sources. 

A manuscript is in preparation.  In the mean time, if you use any of our code or derived data, please cite 
this repository and the [Monarch Initiative](http://www.monarchinitiative.org).  

