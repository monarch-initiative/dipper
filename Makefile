###
### Environment variables.
###

DIPPER_BIN = ./dipper-etl.py
TEST = python3 -m unittest

###
### Tests
###

test: UDP-test IMPC-fetch IMPC-test GWAS-test reactome-test RGD-test \
      string-test trans-test

string-test:
	$(TEST) tests/test_string.py

UDP-test:
	$(TEST) tests/test_udp.py

IMPC-fetch:
	$(DIPPER_BIN) --sources impc --no_verify --fetch_only

IMPC-test:
	$(TEST) tests/test_impc.py

GWAS-test:
	$(TEST) tests.test_gwascatalog.TestGwasSNPModel
	$(TEST) tests.test_gwascatalog.TestGwasHaplotypeModel

reactome-test:
	$(TEST) tests/test_reactome.py

RGD-test:
	$(TEST) tests/test_rgd.py

SGD-test:
	$(TEST) tests/test_sgd.py

trans-test:
	$(TEST) tests/test_trtable.py
