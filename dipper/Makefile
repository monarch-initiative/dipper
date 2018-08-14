# DIPper Ingest  predicate-label check makefile

GH = 
ONTDIR = ../raw/ontology
FORCE = "when web servers do not set last modified headers timestamping is sad"
WGET   = wget --timestamping --directory-prefix=$(ONTDIR)  --tries=3

all: predicate_label.yaml

PHONY: all ONTOLOGIES FORCE

ONTOLOGIES: $(ONTDIR)/sepio.owl $(ONTDIR)/geno.owl $(ONTDIR)/oban_core.ttl \
	$(ONTDIR)/ro.owl  $(ONTDIR)/iao.owl $(ONTDIR)/pco.owl $(ONTDIR)/xco.owl \
	$(ONTDIR)/ero.owl

$(ONTDIR) :
	if [ ! -d $(ONTDIR) ]; then mkdir -p $(ONTDIR);fi

$(ONTDIR)/sepio.owl	: $(ONTDIR)
	$(WGET) 'https://raw.githubusercontent.com/monarch-initiative/SEPIO-ontology/master/src/ontology/sepio.owl'
	rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/sepio.nt

$(ONTDIR)/geno.owl	:  $(ONTDIR)
	$(WGET) 'https://raw.githubusercontent.com/monarch-initiative/GENO-ontology/develop/src/ontology/geno.owl'
	rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/geno.nt

$(ONTDIR)/oban_core.ttl	:  $(ONTDIR)
	$(WGET) 'https://raw.githubusercontent.com/jamesmalone/OBAN/master/ontology/oban_core.ttl'
	rapper -i turtle -o ntriples $@ > $(ONTDIR)//oban_core.nt

$(ONTDIR)/ro.owl	:  $(ONTDIR)
	$(WGET) 'http://purl.obolibrary.org/obo/ro.owl'
	rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/ro.nt
# the ro imports sub-ontologies


$(ONTDIR)/iao.owl	:  $(ONTDIR)
	$(WGET) 'http://purl.obolibrary.org/obo/iao.owl'
	rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/iao.nt

# pco has turned into a (broken) turtle file
$(ONTDIR)/pco.owl	:  $(ONTDIR)
	$(WGET) 'http://purl.obolibrary.org/obo/pco.owl'
	# rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/pco.nt
	rapper -i turtle -o ntriples $@ > $(ONTDIR)/pco.nt

$(ONTDIR)/xco.owl	:  $(ONTDIR)
	$(WGET) 'http://purl.obolibrary.org/obo/xco.owl'
	rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/xco.nt

$(ONTDIR)/ero.owl	:  $(ONTDIR)
	@ # $(WGET) 'http://data.monarchinitiative.org/owl/ero.owl'
	@ # $(WGET) https://www.eagle-i.net/ero/latest/ero.owl
	$(WGET) 'http://purl.obolibrary.org/obo/ero.owl' 
	rapper -i rdfxml -o ntriples $@ > $(ONTDIR)/ero.nt

# owl files may recursively import other owl files
# i.e. the second level starting with the direct imports are
# grep "<owl:imports rdf:resource=" *.owl
#
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/annotations.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/pato_import.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/temporal-intervals.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/go_mf_import.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/core.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/rohom.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/go-biotic.owl"/>
# ro.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/ro/el-constraints.owl"/>
# geno.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/iao/dev/ontology-metadata.owl"/>
# geno.owl:	<owl:imports rdf:resource="http://www.geneontology.org/formats/oboInOwl"/>
# sepio.owl:	<owl:imports rdf:resource="http://purl.obolibrary.org/obo/iao/2015-02-23/ontology-metadata.owl"/>
# sepio.owl:	 <owl:imports rdf:resource="http://purl.obolibrary.org/obo/bfo/2.0/bfo.owl"/>
# 
# for completeness it would be necessary to bring them all in recursively 
# but all I care about here is the subset of labels actually used in dipper.
# 

# isolate the set of labels from the ontologies we use
# and deal with the multiple un authoritative labels 
# (this will get worse with more imports)
# by including both the owl file name and the term prefix

$(ONTDIR)/owl_iri_pfx_trm_lab.unl	:  ONTOLOGIES
	grep  '> <http://www.w3.org/2000/01/rdf-schema#label> "'  $(ONTDIR)/*.nt |\
	cut -f1,3- -d' '| cut --output-delimiter='|' -f1,2 -d'"' | sort -u |\
	sed "s#^$(ONTDIR)/\([a-z_]*\)\.nt:<#\1|#g;s/> |/|/g;s/$$/|/g" |\
	awk -F'|' 'BEGIN{OFS=FS}\
		{split($$2,i,"/");split(i[length(i)],p,"_");\
		$$2=$$2"|"tolower(p[1])"|"p[2];print}' > $@

# other cleanup?  label case? decorations?

# preferred labels will have matching owl and prefix tags
# (nearly half are preferred labels)
$(ONTDIR)/iri_pfx_trm_lab_1.unl : $(ONTDIR)/owl_iri_pfx_trm_lab.unl
	awk -F'|' 'BEGIN{OFS=FS}$$1==$$3{print $$2,$$3,$$4,$$5}' $< > $@

# discounting all terms with preferred labels, 
# what remains will be available as secondary labels (if necessary)
# the terms which only occur once must be fine as is ( > 3k)
# the ones which occur with different labels are a drag (if used in dipper) (< 60)
# arbitrarily choosing one now is not viable because 
# a different choice may already be embedded somewhere.
 
$(ONTDIR)/iri_pfx_trm_lab_2.unl	: $(ONTDIR)/iri_pfx_trm_lab_1.unl \
								  $(ONTDIR)/owl_iri_pfx_trm_lab.unl
	awk -F'|' 'BEGIN{OFS=FS}FNR==NR{a[$$1]++}FNR!=NR{if(!($$2 in a))print}' $^ |\
	cut -f2- -d\| | sort -u > $@



# cut -f1,4 -d \| iri_pfx_trm_lab_2.unl | sort -u | cut -f1 -d \| | uniq -c | sort -nr | grep -v "^ *1 "


predicate_label.yaml	:	$(ONTDIR)/iri_pfx_trm_lab_1.unl \
							$(ONTDIR)/iri_pfx_trm_lab_2.unl
	echo validate_predicate_label.foo $@ $^ 
