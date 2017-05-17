#! /usr/bin/gawk -f

#  Reduce the subject and object of RDF triples (ntriples format)
#  down to their @prefix (or literal object class)
#  Reduce predicates to the specific term identifier
#  and the ontology they are from. (curi form)
#  Express the subject and object as nodes with a directed edge
#  labeled with the predicate in the graphviz dot format
#  Include a tally of each combination of nodes and edge

#  This may over generalize in some cases because I do not have
#  a handy way to differentiate uri for subjects and objects
#  which may belong to a "structural" ontology as opposed to
#  a data uri.

#  Post processing to augment predicate identifiers
#  with their labels seems to improve usefulness.


function usage(){
	print "usage: ntriple2dot.awk curie_map.yaml rdf.nt > rdf.dot"
}



##########################################################
# remove first and last chars from input string <>
function trim(str){
	return substr(str,2,length(str)-2)
}

########################################################
# remove final(ish) element in paths with various delimiters
# this leaves the namespace of a removed local identifier

# this may need to be called more than once since
# some identifiers may have embedded delimiters other uri use differently
function stripid(uri){
	# <_:blanknode> are not allowed in ntriples (but may happen anyway)
	if(0 < match(uri, /^_:|^https:\/\/monarchinitiative.org\/.well-known\/genid\//))
		return "BNODE"
	# perspective endpoints, choose the longest
	delim["_"]=0; delim["/"]=0;
	delim["="]=0; delim[":"]=0; delim["#"]=0;
	char=""; max=-1;
	for(c in delim){
		l = match(uri, char)  # side effect is RLENGTH
		if(l > max){
			char = c;
			delim[char] = RLENGTH
		}
	}
	if(max<=0 || char=="")
		# we don't know what it is, a literal perhaps or uri fragment
		return uri
	else
		return substr(uri,1,delim[char])  # the probably truncated uri
}

# keep underscore, letters & numbers
# change the rest to (a single) underscore sans leading & trailing
# this passes valid node labels from dot's perspective
function simplify(str){
	gsub(/[^[:alpha:][:digit:]_]+/,"_",str)
	gsub(/^_*|_*$/, "",str)
	gsub(/__*/,"_",str)
	return str
}

# if possible, find a shorter form for the input
function contract(uri){
	u = uri
	# shorten till longest uri in curi map is found (or not)
	while(!(u in prefix) &&  length(u))
		u = stripid(substr(u,1,length(u)-1))
	if(u in prefix)
		return prefix[u]
	else{
		# printf("WARN  %s\n", uri) > "/dev/stderr"
		return simplify(uri)
	}
}

# get the final (incl fragment identifier) portion of a slashed path
function final(uri){
	split(uri, b, "/")
	p = b[length(b)]
	anchor = match(p, "#")
	if(anchor > 0)
		p = substr(p, anchor+1)
	return p
}

BEGIN{
	# exceptions
	prefix["BNODE"]="BNODE"  # is a fixed point
	# prefix["https://monarchinitiative.org/_"]="BNODE"
	# just until skolemized bnodes get in the curie map?
	prefix["https://monarchinitiative.org/.well-known/genid"]="BNODE"
	############################################################
	# revisit whether these exceptions are still necessary often
	# they may have been moved into the  prefix mapping file
	prefix["http://www.w3.org/1999/02/22-rdf-syntax-ns#"]="rdf"
	prefix["http://www.w3.org/2000/01/rdf-schema#"]="rdfs"
	prefix["http://www.w3.org/2002/07/owl#"]="owl"
	# in mgi
	prefix["https://www.mousephenotype.org"]="IMPC"
    # http://www.mousephenotype.org/data/genes
	# in panther
	prefix["http://identifiers.org/wormbase"]="WormBase"
	# in IMPC
	prefix["http://www.mousephenotype.org"]="IMPC"
	# in GWAScatalog
	prefix["http://www.ncbi.nlm.nih.gov/SNP/"]="dbSNP"
	# till all non httpS: are purged
	prefix["http://monarchinitiative.org/"]="BASE"
	prefix["http://www.monarchinitiative.org/"]="BASE"
	prefix["http://monarchinitiative.org/MONARCH_"]="MONARCH"
	prefix["http://www.monarchinitiative.org/MONARCH"]=""
	prefix["http://www.monarchinitiative.org/MONARCH_"]="MONARCH"
	prefix["http://data.monarchinitiative.org/ttl"]="MonarchData"
	prefix["http://archive.monarchinitiative.org/ttl"]="MonarchArchive"
}

# main loop
# parse and stash the curie yaml file (first file)
# YAML format is tic delimited word (the curi prefix)
# optional whitespace, a colon, optional whitespace
# then a tic delimited url
# (FNR == NR) && /^'[^']*' *: 'http[^']*'.*/ { # loosing some?
(FNR == NR) && /^'.*/ {
	split($0, arr, "'")
	if(arr[2]=="")
		arr[2]="BASE"
	prefix[arr[4]]=arr[2]
}

# process the ntriple file(s)  which are not the first file

### subject predicate and object are all uri
(FNR != NR) && /^<[^>]*> <[^>]*> <[^>]*> \.$/ {
	s = contract(stripid(trim($1)))
	p =  final(trim($2))
	ns = contract(stripid(trim($2)))
	o = contract(stripid(trim($3)))
	edgelist[s,ns ":" p,o]++
	done[NR]=1
}
### subject & predicate are uri but the object is a literal
(FNR != NR) && /^<[^>]*> <[^>]*> "[^"]*".*\.$/ {
	s = contract(stripid(trim($1)))
	p = final(trim($2))
	ns = contract(stripid(trim($2)))
	edgelist[s, ns ":" p, "LITERAL"]++
	done[NR]=1
}
### subject is a bare blank node, predicate and object are uri
(FNR != NR) && /^_:[^ ]* <[^>]*> <[^>]*> \.$/ {
	# once our bnode syntax is uniform we may want more info out of it
	# s = contract(stripid($1))
	p =  final(trim($2))
	ns = contract(stripid(trim($2)))
	### Object (like subject)
	o = contract(stripid(trim($3)))
	edgelist["BNODE",ns ":" p,o]++
	done[NR]=1
}
### subject & predicate are uri & object is a bare blank node
(FNR != NR) && /^<[^>]*> <[^>]*> _:[^ ]* \.$/ {
	s = contract(stripid(trim($1)))
	p =  final(trim($2))
	ns = contract(stripid(trim($2)))
	#o = contract(stripid($3))
	edgelist[s,ns ":" p,"BNODE"]++
	done[NR]=1
}
### subject and object are a bare blank nodes, predicate is a uri
(FNR != NR) && /^_:[^ ]* <[^>]*> _:[^ ]* \.$/ {
	# s = contract(stripid($1))
	p =  final(trim($2))
	ns = contract(stripid(trim($2)))
	# o = contract(stripid($3))
	edgelist["BNODE",ns ":" p,"BNODE"]++
	done[NR]=1
}
### subject is a bare blank node, predicate is a uri & the object is a literal
(FNR != NR) && /^_:[^ ]* <[^>]*> "[^"]*".*\.$/ {
	# s = contract(stripid($1))
	p = final(trim($2))
	ns = contract(stripid(trim($2)))
	edgelist["BNODE", ns ":" p, "LITERAL"]++
	done[NR]=1
}

# print anything else so it makes itself known
(FNR != NR) && (done[NR]<1) && (NF>1){
	printf("ERROR? in %s at line %i:  %s\n",FILENAME,FNR,$0) > "/dev/stderr"
}

# output dot file, include edge counts
END{
	print "digraph {"
	print "rankdir=LR;"
	print "charset=\"utf-8\";"
	# sort the edges to facilitate future inter dot comparisons
	n=asorti(edgelist, sortlist)
	for (i=1; i<=n; i++) {
		edge=sortlist[i]
		split(edge ,spo, SUBSEP);
		print simplify(spo[1]) " -> " simplify(spo[3]) \
		" [label=\"" spo[2] " (" edgelist[edge] ")\"];"
	}
	print "LITERAL [shape=record];"
	print "labelloc=\"t\";"
	title = final(FILENAME)
	datestamp = strftime("%Y%m%d", systime())
	print "label=\"" substr(title,1,length(title)-3) " (" datestamp ")\";"
	print "}"
}
