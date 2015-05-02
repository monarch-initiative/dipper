###
### Environment variables.
###

DIPPER_BIN = ./dipper.py
NOSE = nosetests

###
### Tests
###

test: BioGrid-fetch BioGrid-test ClinVar-fetch ClinVar-test GeneReviews-fetch GeneReviews-test \
hpoa-fetch hpoa-test ncbi-fetch ncbi-test Panther-fetch Panther-test ucscBands-fetch ucscBands-test \
zfin-fetch zfin-test

BioGrid-fetch:
	$(DIPPER_BIN) --sources biogrid --no_verify --fetch_only

BioGrid-test:
	$(NOSE) tests/test_biogrid.py

ClinVar-fetch:
	$(DIPPER_BIN) --sources clinvar --no_verify --fetch_only

ClinVar-test:
	$(NOSE) tests/test_clinvar.py

GeneReviews-fetch:
	$(DIPPER_BIN) --sources genereviews --no_verify --fetch_only

GeneReviews-test:
	$(NOSE) tests/test_genereviews.py

hpoa-fetch:
	$(DIPPER_BIN) --sources hpoa --no_verify --fetch_only

hpoa-test:
	$(NOSE) tests/test_hpoa.py

ncbi-fetch:
	$(DIPPER_BIN) --sources ncbigene --no_verify --fetch_only

ncbi-test:
	$(NOSE) tests/test_ncbi.py

Panther-fetch:
	$(DIPPER_BIN) --sources Panther --no_verify --fetch_only

Panther-test:
	$(NOSE) tests/test_Panther.py

ucscBands-fetch:
	$(DIPPER_BIN) --sources ucscBands --no_verify --fetch_only

ucscBands-test:
	$(NOSE) tests/test_ucscBands.py

zfin-fetch:
	$(DIPPER_BIN) --sources zfin --no_verify --fetch_only

zfin-test:
	$(NOSE) tests/test_zfin.py
