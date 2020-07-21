MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

# Monarch release version
ARCHIVE ?= https://archive.monarchinitiative.org/beta

WGET = /usr/bin/wget --timestamping --no-verbose

JENA = http://archive.apache.org/dist/jena/binaries/apache-jena-3.15.0.tar.gz
ROBOT_JAR = https://github.com/ontodev/robot/releases/latest/download/robot.jar

BIOLINK_MODEL = https://raw.githubusercontent.com/kshefchek/biolink-model/mappings-for-monarch/biolink-model.ttl
MONARCH_OWL = http://purl.obolibrary.org/obo/upheno/monarch.owl
MONARCH_OWL_MERGED = https://data.monarchinitiative.org/owl/monarch-merged.owl

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

.PHONY: all infer_categories clean

all: infer_categories qc

infer_categories: $(MONARCH_SOURCES) monarch-slim.owl rdf/class-category.ttl $(QUERIES) $(BL_SHIM) $(OUTPUT_RDF)

qc: infer_categories

robot/: ; mkdir $@

jena/: ; mkdir $@

rdf/: ; mkdir $@

blshim/: ; mkdir $@

output/: ; mkdir $@

robot/robot.jar: robot/
	cd $< && $(WGET) $(ROBOT_JAR)

monarch-merged.owl:
	$(WGET) $(MONARCH_OWL_MERGED)

biolink-model.ttl:
	$(WGET) $(BIOLINK_MODEL)

monarch-slim.owl: rdf/ robot/robot.jar monarch-merged.owl biolink-model.ttl
	$(ROBOT) \
      reason \
        --input monarch-merged.owl \
        --reasoner ELK \
      remove \
        --select anonymous \
        --signature false \
      filter \
        --axioms logical \
        --output monarch-slim.owl

rdf/class-category.ttl: monarch-slim.owl
	$(ARQ) \
          --namedgraph biolink-model.ttl \
          --namedgraph monarch-slim.owl \
          --query sparql/bl-categories.sparql \
          --quiet \
          --results ttl > $@

apache-jena-3.15.0.tar.gz:
	$(WGET) $(JENA)

jena/bin/arq: apache-jena-3.15.0.tar.gz jena/
	tar -C jena/ -xf $<

$(MONARCH_SOURCES): rdf/
	cd rdf && $(WGET) $(ARCHIVE)/$@

$(QUERIES):
	sed "s/{DATASET}/$(basename $(notdir $@)).nt/" sparql/infer-categories.sparql.tmpl > $@

$(BL_SHIM): $(QUERIES) $(MONARCH_SOURCES) blshim/ rdf/class-category.ttl jena/bin/arq
	JVM_ARGS="-Xmx60G" $(ARQ) \
	  --namedgraph rdf/class-category.ttl \
	  --namedgraph rdf/$(basename $(notdir $@)).nt \
	  --query sparql/$(basename $(notdir $@)).rq \
	  --quiet \
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
	  monarch-merged.owl \
	  monarch-slim.owl \
	  biolink-model.ttl \
	  apache-jena-3.15.0.tar.gz
