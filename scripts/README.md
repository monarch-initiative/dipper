# scrape-impc.py
USAGE: Gets ontology parameter code to IRI mappings

        sudo apt-get install libffi-dev
        pip install scrapy==1.1.0rc1 cryptography
        scrapy runspider scrape-impc.py -o procedure_mappings.json
        python3 collapse-procedures.py -i procedure_mappings.json -o impc_procedures.json
        mv impc_procedures.json ../resources/
