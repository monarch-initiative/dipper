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
	$(NOSE) --with-coverage --cover-package=dipper tests/test_biogrid.py

ClinVar-fetch:
	$(DIPPER_BIN) --sources clinvar --no_verify --fetch_only

ClinVar-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_clinvar.py

GeneReviews-fetch:
	$(DIPPER_BIN) --sources genereviews --no_verify --fetch_only

GeneReviews-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_genereviews.py

hpoa-fetch:
	$(DIPPER_BIN) --sources hpoa --no_verify --fetch_only

hpoa-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_hpoa.py

ncbi-fetch:
	$(DIPPER_BIN) --sources ncbigene --no_verify --fetch_only

ncbi-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_ncbi.py

Panther-fetch:
	$(DIPPER_BIN) --sources panther --no_verify --fetch_only

Panther-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_panther.py

ucscBands-fetch:
	$(DIPPER_BIN) --sources ucscbands --no_verify --fetch_only

ucscBands-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_ucscbands.py

zfin-fetch:
	$(DIPPER_BIN) --sources zfin --no_verify --fetch_only

zfin-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_zfin.py

kegg-fetch:
	$(DIPPER_BIN) --sources kegg --no_verify --fetch_only

kegg-test:
	$(NOSE) --with-coverage --cover-package=dipper tests/test_kegg.py