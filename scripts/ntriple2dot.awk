#! /usr/bin/awk -f

#  Reduce the subject and object of RDF triples (ntriples format)
# down to their @prefix (or literal object class)
#  Reduce predicates to the specific identifier and
# the ontology they are from.(curi form)
#  Express the subject and object as nodes with a directed edge
# labeled with the predicate in the graphviz dot format
#  Include a tally of each combination of nodes and edge

#  This may over generalize in some cases because I do not have
# a handy way to differentiate uri for subjects and objects
# which may belong to a "structural" ontology as opposed to
# a data uri.

#  Perhaps an improvement would be to also
# express subject and objects as individual curies (ala predicates)
# iff their namespace is also used by a predicate.

##########################################################
# remove first and last chars from input string <>
function trim(str){
	return substr(str,2,length(str)-2)
}

########################################################
# remove final(ish) element in paths with various delimiters
# this leaves the namespace of a removed identifier

# this may need to be called more than once since
# some identifiers may have embedded delimiters
function stripid(uri){
	# <_:blanknode> are not allowed in ntriples (but may happen anyway)
	if(1 == match(uri, /^_:/))
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
	gsub(/^_+|_+$/, "",str)
	gsub(/__*/,"_",str)
	return str
}

# if possible, find a shorter form for the input
function contract(uri){
	u = uri
	# shorten till longest uri in curi map is found (or not)
	while(!(u in prefix))
		u = stripid(substr(u,1,length(u)-1))
	if(u in prefix)
		return prefix[u]
	else
		return simplify(uri)
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
	prefix["https://monarchinitiative.org/_"]="BNODE"
	# just until skolemized bnodes get in the curie map?
	prefix["https://monarchinitiative.org/.well-known/genid"]="BNODE"
	# revisit if exceptions are still necessary
	prefix["http://www.w3.org/1999/02/22-rdf-syntax-ns#"]="rdf"
	prefix["http://www.w3.org/2000/01/rdf-schema#"]="rdfs"
	prefix["http://www.w3.org/2002/07/owl#"]="owl"
	# in mgi
	prefix["https://www.mousephenotype.org"]="IMPC"
	# in panther
	prefix["http://identifiers.org/wormbase"]="WormBase"
	# till all non httpS: are purged
	prefix["http://monarchinitiative.org/"]="BASE"
	prefix["http://www.monarchinitiative.org/"]="BASE"
	prefix["http://monarchinitiative.org/MONARCH_"]="MONARCH"
	prefix["http://www.monarchinitiative.org/MONARCH_"]="MONARCH"
}

# main loop
# parse and stash the curie yaml file (first file)
# YAML format is tic delimited word (the curi prefix)
# optional whitespace, a colon, then a tic delimited url
# (FNR == NR) && /^'[^']*' *: 'http[^']*'.*/ { # loosing some?
(FNR == NR) && /^'.*/ {
	split($0, arr, "'")
	if(arr[2]=="")
		arr[2]="BASE"
	prefix[arr[4]]=arr[2]
}
# process the ntriple file(s)  which are not the first file
### case when subject predicate and object are all uri
(FNR != NR) && /^<[^>]*> <[^>]*> <[^>]*> \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	p =  final(trim($2))
	ns = contract(stripid(trim($2)))
	### Object (like subject)
	o = contract(stripid(trim($3)))
	edgelist[s,ns ":" p,o]++
}
### case when the object is a literal
(FNR != NR) && /^<[^>]*> <[^>]*> "[^"]*" \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	p = final(trim($2))
	ns = contract(stripid(trim($2)))
	### not a uri
	o = "LITERAL"
	edgelist[s, ns ":" p, o]++
	nodelist[o " [shape=record];"]++
}

# output dot file, include edge counts
END{
	print "digraph {"
	print "rankdir=LR;"
	print "charset=\"utf-8\";"
	for(edge in edgelist){
		split(edge ,spo, SUBSEP);
		print simplify(spo[1]) " -> " simplify(spo[3]) \
		" [label=\"" spo[2] " (" edgelist[edge] ")\"];"
	}
	for(node in nodelist) print node
	print "labelloc=\"t\";"
	title = final(FILENAME)
	datestamp = strftime("%Y%m%d", systime())
	print "label=\"" substr(title,1,length(title)-3) " (" datestamp ")\";"
	print "}"
}

