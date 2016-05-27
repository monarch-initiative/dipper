# scrape-impc.py
USAGE: Gets ontology parameter code to IRI mappings

        sudo apt-get install libffi-dev
        pip install scrapy==1.1.0rc1 cryptography
        scrapy runspider scrape-impc.py -o parameter_mappings.json
        python3 collapse-parameters.py -i parameter_mappings.json -o impc_parameters.json
        mv impc_parameters.json ../resources/
