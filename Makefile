###
### Environment variables.
###

DIPPER_BIN = ./dipper-etl.py --debug
TEST = python -m unittest

VENV = $(shell python ./scripts/within_pip_env.py)
# $(warning python virtual env detection report $(VENV))

NP=4

PYUTEST = @([ $(VENV) -eq 0 ] && $(TEST) tests/$@.py)||echo "Not in python virtual enviroment"
GTT = "translationtable/GLOBAL_TERMS.yaml"

all: lint_error test tt_generated

###
### Tests
###

test: test_translationtable test_impc test_reactome test_clinvar test_wormbase \
	test_rgd test_ctd test_string test_udp test_orphanet test_mgi test_gwascatalog \
	test_dataset test_source_metadata # test_impc_fetch

test_dataset:
	$(PYUTEST)

test_source_metadata:
	$(PYUTEST)

test_mgi:
	$(PYUTEST)

test_clinvar:
	$(PYUTEST)

test_flybase:
	$(PYUTEST)

test_wormbase:
	$(PYUTEST)

test_string:
	$(PYUTEST)

test_udp:
	$(PYUTEST)

test_impc:
	$(PYUTEST)

test_gwascatalog:
	$(PYUTEST)

test_reactome:
	$(PYUTEST)

test_rgd:
	$(PYUTEST)

test_sgd:
	$(PYUTEST)

test_ctd:
	$(PYUTEST)

test_mychem:
	$(PYUTEST)

test_ncbi:
	$(PYUTEST)

test_omim:
	$(PYUTEST)

test_biogrid:
	$(PYUTEST)

test_orphanet:
	$(TEST) tests.test_orphanet.GeneVariantDiseaseTest

test_impc_fetch:
	$(DIPPER_BIN) --sources impc --no_verify --fetch_only

test_omia-integration:
	python tests/omia-integration.py --input ./out/omia.ttl

###################################################################################
###  checks on supporting artefacts

test_translationtable:
	@ echo  "lint curie_map yaml file"
	@ yamllint -c .yamllint dipper/curie_map.yaml
	@ echo "----------------------------------------------------------------------"
	@ echo  "lint translation table yaml files"
	@ yamllint -c .yamllint translationtable/
	@ echo "----------------------------------------------------------------------"
	@ echo "python unit test for duplicate keys and invertablility in global tt"
	@ $(TEST) tests/test_trtable.py
	@ echo "----------------------------------------------------------------------"
	@ echo "Is the _entire_ $(GTT) file ordered by ontology curie?"
	@ sort -f -k2,2 -k3,3n -t ':' --stable --check $(GTT)
	@ echo "----------------------------------------------------------------------"
	@ echo "Orphan labels in dipper/source/Ingest.py  w.r.t. $(GTT)?"
	@ scripts/check_labels_v_gtt.sh
	@ echo "----------------------------------------------------------------------"


# lint if within a python virtual env
# make can be called without being in a venv and produce unhelpful results
lint_error:
	@ # runs single threaded just under a minute with 4 cores ~ 20 seconds
	@ echo "Lint for errors"
	@ ([ $(VENV) -eq 0 ] && time pylint -j $(NP) -E ./dipper)||\
		echo "Not in python virtual enviroment"
	@ echo "----------------------------------------------------------------------"

lint_warn:
	@ echo "Lint for warnings"
	@ ([ $(VENV) -eq 0 ] && time pylint  -j $(NP) --disable=C,R ./dipper)||\
		echo "Not in python virtual enviroment"
	@ echo "----------------------------------------------------------------------"

lint:
	@ echo "Lint for everything"
	@ ([ $(VENV) -eq 0 ] && time pylint  -j $(NP) ./dipper)||\
		echo "Not in python virtual enviroment"
	@ echo "----------------------------------------------------------------------"

# Generate specialized files from our various mapping files

tt_generated:  translationtable/generated/prefix_equivalents.yaml

clean_tt_generated:
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
				$^ | LC_ALL=C sort -u > $@

translationtable/generated/prefix_equivalents.yaml: \
		translationtable/generated/curiemap_prefix.txt /tmp/local_inverse.tab
	@ echo "---\n# prefix_equivalents.yaml" > $@;
	@ LC_ALL=C join translationtable/generated/curiemap_prefix.txt /tmp/local_inverse.tab|\
	awk '{v=$$1;$$1="";print substr($$0,2) ": " v}' | sort -u >> $@
