###
### Environment variables.
###

DIPPER_BIN = ./dipper-etl.py --debug
TEST = python3 -m unittest

all: test prefix_equivalents

###
### Tests
###

test: trans-test IMPC-fetch IMPC-test  reactome-test RGD-test CTD-test mychem-test \
      string-test  UDP-test  # Orphanet-test MGI-test GWAS-test

MGI-test:
	$(TEST) tests.test_mgi.EvidenceTestCase

Orphanet-test:
	$(TEST) tests.test_orphanet.GeneVariantDiseaseTest

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

CTD-test:
	$(TEST) tests/test_ctd.py

mychem-test:
	$(TEST) tests/test_mychem.py

ncbi-test:
	$(TEST) tests/test_ncbi.py

OMIM-test:
	$(TEST) tests/test_omim.py

GTT = "translationtable/GLOBAL_TERMS.yaml"
trans-test:
	# python unit test for duplicate keys and invertablility
	$(TEST) tests/test_trtable.py
	# is $(GTT) table ordered by ontology curie
	@ if sort -k 2,3n -t ':' -s -c $(GTT) ; then echo "Order is okay" ; fi
	# Are there terms found in the source, but not found in $(GTT)
	@ scripts/check_labels_v_gtt.sh

omia-int-test:
	python tests/omia-integration.py --input ./out/omia.ttl

# Generate specalized files from our various mapping files

prefix_equivalents:  translationtable/generated/prefix_equivalents.yaml

clean_prefix_equivalents:
	rm translationtable/generated/prefix_equivalents.yaml
	rm translationtable/generated/curiemap_prefix.txt
	rm /tmp/local_inverse.tab

translationtable/generated/curiemap_prefix.txt: dipper/curie_map.yaml
	@ cut -f1 -d ':' dipper/curie_map.yaml  | tr -d "'" | egrep -v "^$|^ *#" |\
		grep .|sed 's|\(.*\)|"\1"|g'| LANG=en_EN sort -u > \
			translationtable/generated/curiemap_prefix.txt

/tmp/local_inverse.tab: translationtable/[a-z_-]*.yaml
	@ awk -F '"' '/^"[^"]+": "[^":]+".*/\
		{if($$2 != $$4 && ! match($$2, /[0-9]+/))\
			print "\"" $$4 "\"\t\"" $$2 "\""}' \
				translationtable/[a-z_-]*.yaml | LANG=en_EN sort -u > \
					/tmp/local_inverse.tab

translationtable/generated/prefix_equivalents.yaml: translationtable/generated/curiemap_prefix.txt /tmp/local_inverse.tab
	@ LANG=en_EN join translationtable/generated/curiemap_prefix.txt  /tmp/local_inverse.tab  |\
		awk '{v=$$1;$$1="";print $$0 ": " v}'| sort -u  > \
			translationtable/generated/prefix_equivalents.yaml
