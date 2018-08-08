###
### Environment variables.
###

DIPPER_BIN = ./dipper-etl.py --debug
TEST = python3 -m unittest

###
### Tests
###

test: MGI-test UDP-test IMPC-fetch IMPC-test GWAS-test reactome-test \
      string-test trans-test CTD-test mychem-test
# RGD-test
MGI-test:
	$(TEST) tests.test_mgi.EvidenceTestCase

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

# temp remove till TP's pr is in
#RGD-test:
#	$(TEST) tests/test_rgd.py

SGD-test:
	$(TEST) tests/test_sgd.py

CTD-test:
	$(TEST) tests/test_ctd.py

mychem-test:
	$(TEST) tests/test_mychem.py

trans-test:
	$(TEST) tests/test_trtable.py
