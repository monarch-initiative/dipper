###
### Environment variables.
###

DIPPER_BIN = ./dipper.py
NOSE = nosetests

###
### Tests
###

test: BioGrid-fetch ClinVar-fetch GeneReviews-fetch hpoa-fetch ncbi-fetch \
Panther-fetch ucscBands-fetch zfin-fetch dipper-test

BioGrid-fetch:
	$(DIPPER_BIN) --sources biogrid --no_verify --fetch_only

ClinVar-fetch:
	$(DIPPER_BIN) --sources clinvar --no_verify --fetch_only

GeneReviews-fetch:
	$(DIPPER_BIN) --sources genereviews --no_verify --fetch_only

hpoa-fetch:
	$(DIPPER_BIN) --sources hpoa --no_verify --fetch_only

ncbi-fetch:
	$(DIPPER_BIN) --sources ncbigene --no_verify --fetch_only

Panther-fetch:
	$(DIPPER_BIN) --sources panther --no_verify --fetch_only

ucscBands-fetch:
	$(DIPPER_BIN) --sources ucscbands --no_verify --fetch_only

zfin-fetch:
	$(DIPPER_BIN) --sources zfin --no_verify --fetch_only

dipper-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_biogrid.py tests/test_clinvar.py \
tests/test_genereviews.py tests/test_hpoa.py tests/test_ncbi.py tests/test_panther.py \
tests/test_ucscbands.py tests/test_zfin.py
