#! /usr/bin/gawk -f

# label_dot_edge.awk
# Attempts to augment opaque ontological terms with their current label
# example usage:
# label_dot_edge.awk predicate_curie_term_label.tab  in.dot > out.dot
#
# Arguments: two files
# Loads the first file as a map from terms to labels:
# (expects: [namespace:]<tab><term><tab><label>[<tab><whatever>])
# the first set was generated via cypher queries against our scigraph
# (neo4j) instance loaded with the ontologies we use.
# Reads the second Graphviz digraph dot format file:
# looks for <term> as edge-text 
# and appends "! <label> (Iff <term> is found in map)
# note: This "<term> ! <label>" structure is a convention used by Ontologists.


BEGIN{FS="\t"; OFS="\t"}
# load labels
FNR == NR && $3 != "null"{label[$1 $2]=$3}

# output lines that do not get labels
FNR != NR && !/.* -> .*label=.*/ {print}

# output lines w/labels if they are available
FNR != NR && /.* -> .*label=.*/ {
	line="";
	split($0, row, / /)
	split(row[4], term, /"/)
	if (term[2] in label){
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
