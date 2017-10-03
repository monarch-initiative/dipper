###
### Environment variables.
###

DIPPER_BIN = ./dipper-etl.py
NOSE = nosetests

###
### Tests
###

test: UDP-test IMPC-fetch IMPC-test GWAS-test   \
CTD-test interactions-test reactome-test RGD-test \
string-test

string-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_string.py

UDP-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_udp.py

IMPC-fetch:
	$(DIPPER_BIN) --sources impc --no_verify --fetch_only

IMPC-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_impc.py

GWAS-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_gwascatalog.py

CTD-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_ctd.py

interactions-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_interactions.py

reactome-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_reactome.py

RGD-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_rgd.py


