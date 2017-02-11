#! /usr/bin/awk -f
#deltadot.awk
# Use to illuminate the commonality and differences between dot files created 
# by 'ntriples2dot.awk' but might work for other dot files 

# expects two graphviz dot digraph files 
# with simple node names ('C' variable naming conventions)
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
			# TODO flag losses?
			diff =  edge2[e] - edge1[e];
			#if(diff < 0) diff = "<b>" diff "</b>";
			edges[e "(" diff ")\", color=\"black\"];"]
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

