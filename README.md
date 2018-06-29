[![PyPI](https://img.shields.io/pypi/v/dipper.svg)](https://pypi.python.org/pypi/dipper)
[![Build Status](https://travis-ci.org/monarch-initiative/dipper.svg?branch=master)](https://travis-ci.org/monarch-initiative/dipper)

# Dipper

Dipper is a Python package to generate RDF triples from common scientific resources.
Dipper includes subpackages and modules to create graphical models of this data, including:

* Models package for generating common sets of triples, including common OWL axioms, complex genotypes, associations, evidence and provenance models.
* Graph package for building graphs with RDFLib or streaming n-triples
* Source package containing fetchers and parsers that interface with remote databases and web services


* The dipper-etl.py script wraps all of the source parsers, enabling users to specify one or more sources to process.
The general strategy is that there is one class per data source.  We define the files to be fetched,
any file scrubbing, and then the parsing methods.  As the files are parsed, triples are loaded into an in-memory graph.
This graph is then typically dumped into triples in turtle format.  For testing purposes,
 a subset of the graph is also dumped to *_test.ttl.
* Data generated from this pipeline can be used in a variety of ways downstream.  We recommend
loading the data into a triple store or graph database that is optimized for use with ontologies, such as
[BlazeGraph](https://github.com/blazegraph/database).  We also maintain [SciGraph](https://github.com/SciGraph), an application that loads RDF and OWL into Neo4J.
Smaller files can be loaded into an ontology editor like [Protege](http://protege.stanford.edu/).

## Installing Dipper:
Dipper requires [Python 3.5](https://www.python.org/downloads/) or higher.


* To run the dipper pipeline, or use it as a python module, install the latest stable version with pip:

    ```pip3 install dipper```

* To install the development branch, clone the repository and run:

    ```pip3 install -e /path/to/git/dipper```

* Or alternatively without cloning,

    ```pip3 install git+git://github.com/monarch-initiative/dipper.git```

## Getting started:
* you can run the code by supplying a list of one or more sources on the command line.  some examples:

    ```dipper-etl.py --sources impc,hpoa```

* furthermore, you can check things out by supplying a limit.  this will only process the
first N number of rows or data elements

    ```dipper-etl.py --sources hpoa --limit 100```

* you can also run the stand-alone tests in ```tests/test_*``` to generate subsets of the data and run unittests
* other commandline parameters are explained if you request help:

    ```dipper-etl.py --help```


## Building locally
To build locally, clone this repo and install the requirements using pip.

* Required external python packages can be found in the [requirements.txt](requirements.txt)

* Optional source specific python packages can be found in [./requirements/](requirements)
    
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

## Documentation:

The full documentation, including API docs, can be found on [read the docs](https://dipper.readthedocs.io).


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
    * MyDrug (FAERS drug-outcome associations)
    * GeneOntology (gene-process/function/subcellular location associations)
    * Monarch (Monarch curated genotype-phenotyp associations)
    * Monochrom (Ontology of chromosomes)
    * Orphanet (gene-disease associations)
    * UCSCBands (RDF representation of chromosomal bands using FALDO an Monochrom)
    * String (protein-protein interactions)
    * OMA (orthologs from QfO reference proteomes 2017 (79 species))
    * RGD (gene to phenotype)
    * SGD (gene to phenotype)
    * MyChem (targets, interactions, indications)
    
    Each source has a corresponding script at https://github.com/monarch-initiative/dipper/tree/master/dipper/sources

* Many sources have a corresponding **concept map** diagram that documents modeling patterns implemented in SciGraph, via Dipper-mediated transformation into Monarch's common target model. These are stored in the [ingest-artifacts repo](https://github.com/monarch-initiative/ingest-artifacts/tree/master/sources).

* Don't see a parser you want?  Feel free to request a new one, or you could contribute a Source parser to our suite!
Please see our [best-practices documentation](dipper/sources/README.md) for details on writing new Source parsers
using Dipper code, and make a pull request.

## Identifiers
Our identifier documentation as referenced in our recent paper on identifiers(doi:10.1371/journal.pbio.2001414)[https://doi.org/10.1371/journal.pbio.2001414]


## About this project
The Dipper data pipeline was born out of the need for a uniform representation of human and model organism
genotype-to-phenotype data, and an Extract-Transform-Load (ETL) pipeline to process it all.
It became too cumbersome to first get all of these data into a relational schema; so, we decided to go straight from each source into triples that
are semantically captured, using standard modeling patterns.  We are currently working on tooling around
defining, documenting, and constraining our schema as [biolink models](https://github.com/biolink/biolink-model).

A manuscript is in preparation.  In the meantime, if you use any of our code or derived data, please cite
this repository and the [Monarch Initiative](https://monarchinitiative.org).
