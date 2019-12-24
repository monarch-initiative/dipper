#! /usr/bin/awk -f
# deltadot.awk <encumbant.gv> <candidate.gv>
# Use to illuminate the commonality and differences between dot files created
# by 'ntriples2dot.awk' but might work for other dot files
#
# To ignore counts
# deltadot.awk -v"NOCOUNTS=1" <encumbant.gv> <candidate.gv>
#
# expects two graphviz dot digraph files;
# by convention "old"  then "new"
# for color to mean;  orange lost in "new" and blue gained in "new"
# with simple node names ('C-lang' variable naming conventions)
# edge's unadorned predicate label
# concludes with a space then count in parens
#
# For example:
# NODE1 -> NODE2 [label="curie:term (42)"];
#
# filter for edges in each file associate the count with the edge
# if an edge is common between the files
# then change its "count" to the difference. (second minus first)
# otherwise leave count as is
# UNLESS there is a parameter NOCOUNTS not equal to 0
# in which case set all counts the same
# this is for use when the counts have no meaning
# such as with random samples


# associate the count with an edge
function parse(str, arr){
	split(str, a, "(")
	arr[a[1]] = substr(a[2], 1, index(a[2], ")")-1)
}

# strip metadata from filename
function de_path_ext(pth){
	split(pth, a, "/")
	return substr(a[length(a)], 1, index(a[length(a)], ".")-1)
}
BEGIN { NOCOUNTS=0}

# collect the edges
NR==FNR && / -> / {
	label1 = FILENAME;
	parse($0, edge1)
}
NR!=FNR && / -> / {
	label2 = FILENAME;
	parse($0, edge2)
}

END{
	for(e in edge1){
		if(!(e in edge2)) edges[e "(" edge1[e] ")\", color=\"orange\"];"];
		else {
			# edge is in both graphs
			if(NOCOUNTS)
				edges[e "()\", color=\"black\"];"]
			else {
				diff = edge2[e] - edge1[e];
				if(diff < 0)
					diff = "<font color=\"red\">" diff "</font>";
				edges[e "(" diff ")\", color=\"black\"];"]
			}
			delete edge2[e]
		}
	}
	for(e in edge2) edges[e "(" edge2[e] ")\", color=\"blue\"];"];

	print "digraph {"
	print "rankdir=LR;"
	print "charset=\"utf-8\";"
	print "LITERAL [shape=record];"
	print "labelloc=\"t\";"
	print "label=\"" de_path_ext(label1) " <=> " de_path_ext(label2) "\""

	n = asorti(edges,output)
	for(i=1;i<=n;i++ ){print output[i]}
	print "\n}\n"
}

