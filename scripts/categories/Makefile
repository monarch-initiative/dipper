MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

# Monarch release version
ARCHIVE ?= https://archive.monarchinitiative.org/beta

WGET = /usr/bin/wget --timestamping --no-verbose

JENA = https://apache.osuosl.org/jena/binaries/apache-jena-3.15.0.tar.gz
ROBOT_JAR = https://github.com/ontodev/robot/releases/latest/download/robot.jar

BIOLINK_MODEL = https://raw.githubusercontent.com/biolink/biolink-model/master/biolink-model.ttl
MONARCH_OWL = http://purl.obolibrary.org/obo/upheno/monarch.owl

ROBOT := java -jar robot/robot.jar
ARQ := jena/apache-jena-3.15.0/bin/arq

SOURCES =  animalqtldb \
		biogrid \
		bgee \
		clinvar \
		coriell \
		ctd \
		ensembl \
		eom \
		flybase \
		genereviews \
		go \
		gwascatalog \
		hgnc \
		hpoa \
		impc \
		kegg \
		mgi \
		mgislim \
		mmrrc \
		monarch \
		monochrom \
		mpd \
		ncbigene \
		omia \
		omim \
		orphanet \
		panther \
		reactome \
		rgd \
		sgd \
		string \
		ucscbands \
		udp \
		wormbase \
		zfin \
		zfinslim

MONARCH_SOURCES = $(foreach src, $(SOURCES), rdf/$(src).nt)
BL_SHIM = $(foreach src, $(SOURCES), blshim/$(src).nt)
QUERIES = $(foreach src, $(SOURCES), sparql/$(src).rq)
OUTPUT_RDF = $(foreach src, $(SOURCES), output/$(src).nt)

.PHONY: all clean

all: $(MONARCH_SOURCES) rdf/class-category.ttl $(QUERIES) $(BL_SHIM) $(OUTPUT_RDF)

robot/: ; mkdir $@

jena/: ; mkdir $@

rdf/: ; mkdir $@

blshim/: ; mkdir $@

output/: ; mkdir $@

robot/robot.jar: robot/
	cd $< && $(WGET) $(ROBOT_JAR)

monarch.owl:
	$(WGET) $(MONARCH_OWL)

biolink-model.ttl:
	$(WGET) $(BIOLINK_MODEL)

rdf/class-category.ttl: rdf/ robot/robot.jar monarch.owl biolink-model.ttl
	$(ROBOT) \
	  merge \
        --input monarch.owl \
      remove \
        --axioms "disjoint" \
      remove \
        --term owl:Nothing \
      reason \
        --reasoner ELK \
      remove \
        --select anonymous \
        --signature false \
      filter \
        --axioms logical \
      merge \
        --input biolink-model.ttl \
      query \
        --query ./sparql/bl-categories.sparql \
        rdf/class-category.ttl

apache-jena-3.15.0.tar.gz:
	$(WGET) $(JENA)

jena/bin/arq: apache-jena-3.15.0.tar.gz jena/
	tar -C jena/ -xf $<

$(MONARCH_SOURCES): rdf/
	cd rdf && $(WGET) $(ARCHIVE)/$@

$(QUERIES):
	sed "s/{DATASET}/$(basename $(notdir $@)).nt/" sparql/infer-categories.sparql.tmpl > $@

$(BL_SHIM): $(QUERIES) $(MONARCH_SOURCES) blshim/ rdf/class-category.ttl jena/bin/arq
	JVM_ARGS="-Xmx50G" $(ARQ) \
	  --namedgraph rdf/class-category.ttl \
	  --namedgraph rdf/$(basename $(notdir $@)).nt \
	  --query sparql/$(basename $(notdir $@)).rq \
	  --results nt > $@

$(OUTPUT_RDF): $(BL_SHIM) output/
	cat rdf/$(basename $(notdir $@)).nt blshim/$(basename $(notdir $@)).nt > $@

clean:
	rm --force --recursive --verbose --preserve-root --one-file-system \
	  robot/ \
	  jena/ \
	  rdf/ \
	  blshim/ \
	  sparql/*.rq \
	  output/ \
	  monarch.owl \
	  biolink-model.ttl \
	  apache-jena-3.15.0.tar.gz
