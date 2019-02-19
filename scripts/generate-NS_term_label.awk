#! /usr/bin/awk -f
# ./generate-NS_term_label.awk \
# ${DIPPER}/dipper/curie_map.yaml iri_label-final.tab \
#  > ns_curie_label_iri.tab

# Recreates a format that neo4J produced which became accidentally adopted.
# Used for creating readable labels in GraphViz summaries of
# Monarch's Dipper RDF output.

BEGIN {OFS="\t"}
# dipper/curie_map.yaml    'pfx': 'buri'
NR == FNR && !/^#/{
	FS="'"
	pfx[$4]=$2
}

# iri_label-final.tab  iriterm label
NR != FNR{
	FS="\t"
	split($1, t,"/");
	term = t[length(t)];
	# some terms are page anchors
	term = substr(term, index(term,"#")+1)
	buri=substr($1,1,index($1,term)-1)
	print pfx[buri] ":", term, $2, $1
}
