#! /usr/bin/gawk -f

# label_dot_edge.awk
# Augment ontological terms with their label
# example usage:
# ${DIPPER}/scripts/label_dot_edge.awk
#	${DIPPER}/translationtable/GLOBAL_TERMS.yaml  in.dot > out.dot
#
# Arguments: two files
# Loads the first file as a map from  to curie to label :
 # (expects: "label": "curie",)
# Reads the second Graphviz digraph dot format file:
# looks for <term> as edge-text
# and appends "! <label> (Iff <term> is found in map)
# note: This "<term> ! <label>" structure is a convention used by Ontologists.


BEGIN{FS="\""; OFS="\t"}
# load labels
FNR == NR && /^"/{label[$4]=$2}

# output lines that do not get labels
FNR != NR && !/.* -> .*label=.*/ {print}

# output lines w/labels if they are available

FNR != NR && /.* -> .*label=.*/ {
	line="";
	split($0, row, / /)

	split(row[4], term, /"/)
	split(term[2], lclid, ":")

	if (term[2] in label && lclid[2] != label[term[2]]){
		term[2] = term[2] " ! " label[term[2]];
		row[4] = term[1];
		for(i=2;i<=length(term);i++)
			row[4] = row[4] "\"" term[i]
	}
	line = row[1]
	for(i=2;i<=length(row);i++)
		line = line " " row[i]
	print  line
}
