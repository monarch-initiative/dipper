## scrape-impc.py
USAGE: Gets ontology parameter code to IRI mappings

        sudo apt-get install libffi-dev
        pip install scrapy==1.1.0rc1 cryptography
        scrapy runspider scrape-impc.py -o parameter_mappings.json
        python3 collapse-parameters.py -i parameter_mappings.json -o impc_parameters.json
        mv impc_parameters.json ../resources/
        

## fetch-gene-ids.py
Convert HGNC gene symbols to entrez gene with curie prefix NCBIGene,
using the mygene API and monarch vocab services and collapsing results

USAGE ./scripts/fetch-gene-ids.py --input ./resources/mappings/gene.tsv --output ./gene.tsv


* Required python packages:
    * [Mygene](http://mygene-py.readthedocs.org/en/latest/)

