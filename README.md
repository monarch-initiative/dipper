[![Build Status](https://travis-ci.org/monarch-initiative/dipper.svg?branch=master)](https://travis-ci.org/monarch-initiative/dipper)
[![Coverage Status](https://coveralls.io/repos/monarch-initiative/dipper/badge.svg)](https://coveralls.io/r/monarch-initiative/dipper)

# DIPPER

Dipper is a pure Python package to generate RDF triples from common scientific resources.
Dipper includes subpackages and modules to create graphical models of this data, including:

* Models package for generating common sets of triples, including common OWL axioms, complex genotypes, associations, evidence and provenance models.
* Graph package for building graphs with RDFLib or streaming n-triples
* Source package containing fetchers and parsers that interface with remote databases and web services


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
* [Python 3](https://www.python.org/downloads/) or higher (and therefore pip3 if using pip)
* One of the unit tests requires
[owltools](https://github.com/owlcollab/owltools) be available on your path.  You could modify
the code to skip this, if necessary
* Running make test requires nosetests (if on OS X you may need to `sudo pip3 install nose`)

* Required external python packages:
    * [rdflib](https://pypi.python.org/pypi/rdflib)
    * isodate
    * roman
    * pyyaml

    
* Optional source specific python packages:
    * [psycopg2](http://initd.org/psycopg/)
    * [python-docx](https://github.com/python-openxml/python-docx)
    * beautifulsoup4
    * GitPython
    * intermine
    * pysftp
    * [Requests](http://requests.readthedocs.org/en/master/)
    
Note, Dipper imports source modules dynamically at runtime.  As a result it is possible to build a core set
of requirements and add source specific dependencies as needed.  Presently this only implemented with pip requirements
files. For example to build dependencies for MGI:

        pip3 install -r requirements.txt
        pip3 install -r requirements/mgi.txt

To install dependencies for all sources:

        pip3 install -r requirements.txt
        pip3 install -r requirements/all-sources.txt
    
If you encounter any errors installing these packages using Homebrew, it could be due to [a curent known issue in upgrading to  pip3](https://github.com/Homebrew/homebrew/issues/25752). In this case, first force reinstall pip2 (````pip2 install --upgrade --force-reinstall pip````) and then install the package using pip3 (eg. ````pip3 install psycopg2````.)


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

## Sources:
* The following sources have been mapped:
    * Human Phenotype Ontology Annotations (Disease-Phenotype abnormal associations)
    * BioGrid (Gene-gene interactions)
    * Bgee (Gene expression in anatomical entities)
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
    * MMRRC (Strains and curated phenotypes)
    * ENSEMBL (genes)
    * AnimalQTLdb (QTLs, genomic locations, and phenotypic traits)
    * GWASCatalog (SNPs and associated phenotypes/factors tested)
    * HGNC (genes)
    * KEGG (genes, diseases, drugs)
    * MPD (strains and phenotypes derived from outlier measurements (>2s.d.))
    * OMIA (non-laboratory animal phenotypes)
    * Wormbase (genes, alleles, phenotypes)
    * FlyBase (genotype, phenotype)
    * Reactome (gene-pathway associations)
    * MyDrug (drug-outcome associations)
    * GeneOntology (gene-process/function/subcellular location associations)
    * Monarch (Monarch curated genotype-phenotyp associations)
    * Monochrom (Ontology of chromosomes)
    * Orphanet (gene-disease associations)
    * UCSCBands (RDF representation of chromosomal bands using FALDO an Monochrom)
    * WormBase (genotype-phenotype associations)
    
    Each source has a corresponding script at https://github.com/monarch-initiative/dipper/tree/master/dipper/sources

   ```
   hpoa,zfin,omim,biogrid,mgi,impc,panther,ncbigene,ucscbands,
   ctd,genereviews,eom,coriell,clinvar,monochrom,kegg,animalqtldb,
   ensembl,hgnc,orphanet,omia,flybase,mmrrc,wormbase,mpd,gwascatalog,go
   ``` 
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
this repository and the [Monarch Initiative](https://monarchinitiative.org).

## Identifiers
Throughout the Monarch web application, we display external entities using their human-friendly labels
(eg. ontology term label 'polydactyly' or gene symbol 'KNG1') as issued by the original data sources;
however, while such labels aid human understanding, they often overlap between sources.
Therefore, in Monarch, we never rely on the labels to integrate data and never display labels alone without a
corresponding prefixed identifier (wherein the local part is exactly as issued by the original data sources and the
prefix is as established by convention or as registered. eg. NCBIGene:3827).

For each prefix we display in Monarch, we have [documented a 1-to-1 relationship with a resolving namespace](https://github.com/monarch-initiative/dipper/blob/master/dipper/curie_map.yaml),
and the prefixed notation (aka CURIE) is usually hyperlinked to its HTTP URI.
For more information regarding identifiers terminology and notation, see McMurry et al. https://zenodo.org/record/31765.

More detailed identifier documentation for Monarch is a work in progress, available [here:](https://docs.google.com/document/d/1jJHM0c358T5h2W2qLbpm9fGNcOsTSfhMPmmXQhI8n9Q/edit)
Please feel free to pose any questions or concerns to info@monarchinitiative.org.


![in action](https://github.com/monarch-initiative/dipper/blob/master/docs/curies-and-uris-in-action.png)

