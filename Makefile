###
### Environment variables.
###

DIPPER_BIN = ./dipper.py
NOSE = nosetests

testFetch:
	$(DIPPER_BIN) --sources biogrid --no_verify --fetch_only

###
### Tests
###

test:
	$(NOSE) tests/test_biogrid.py
