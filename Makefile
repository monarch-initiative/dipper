###
### Environment variables.
###

DIPPER_BIN = ./dipper-etl.py --debug
TEST = python3 -m unittest

all: test prefix_equivalents dot_to_svg

###
### Tests
###

test: ClinVar-test FlyBase-test WormBase-test trans-test IMPC-test reactome-test \
      RGD-test CTD-test mychem-test string-test UDP-test Orphanet-test MGI-test \
      GWAS-test # IMPC-fetch

MGI-test:
	$(TEST) tests.test_mgi.EvidenceTestCase

ClinVar-test:
	$(TEST) tests/test_clinvar.py

FlyBase-test:
	$(TEST) tests/test_flybase.py

WormBase-test:
	$(TEST) tests/test_wormbase.py

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

biogrid-test:
	$(TEST) tests/test_biogrid.py

GTT = "translationtable/GLOBAL_TERMS.yaml"
trans-test:
	@ echo  "Is yamllint installed?"
	@ yamllint --version
	@ echo  "Are TT files valid?"
	@ # note:  sudo apt-get install yamllint
	@ yamllint -c translationtable/.yamllint translationtable/ && echo "TT files Lint OKAY."
	@ echo "python unit test for duplicate keys and invertablility"
	@ $(TEST) tests/test_trtable.py
	@ echo "Is $(GTT) table ordered by ontology curie?"
	@ sort -k2,2 -k3,3n -t ':' --stable --check $(GTT) && echo "GLOBALTT order is OKAY"
	@ echo "Are terms found in dipper/source/Ingest.py, \nbut not in $(GTT)?"
	@ scripts/check_labels_v_gtt.sh

omia-int-test:
	python tests/omia-integration.py --input ./out/omia.ttl

# Generate specalized files from our various mapping files@

prefix_equivalents:  translationtable/generated/prefix_equivalents.yaml

clean_prefix_equivalents:
	rm translationtable/generated/prefix_equivalents.yaml
	rm translationtable/generated/curiemap_prefix.txt
	rm /tmp/local_inverse.tab

translationtable/generated/curiemap_prefix.txt: dipper/curie_map.yaml
	@ cut -f1 -d ':' dipper/curie_map.yaml  | tr -d "'" | egrep -v "^$|^ *#" |\
		grep .|sed 's|\(.*\)|"\1"|g' | LC_ALL=C sort -u > $@

/tmp/local_inverse.tab: translationtable/[a-z_-]*.yaml
	@ awk -F '"' '/^"[^"]+": "[^":]+".*/\
		{if($$2 != $$4 && ! match($$2, /[0-9]+/))\
			print "\"" $$4 "\"\t\"" $$2 "\""}' \
				translationtable/[a-z_-]*.yaml | LC_ALL=C sort -u > $@

translationtable/generated/prefix_equivalents.yaml: \
		translationtable/generated/curiemap_prefix.txt /tmp/local_inverse.tab
	@ echo "---\n# prefix_equivalents.yaml" > $@;
	@ LC_ALL=C join translationtable/generated/curiemap_prefix.txt /tmp/local_inverse.tab|\
	awk '{v=$$1;$$1="";print substr($$0,2) ": " v}' | sort -u >> $@
