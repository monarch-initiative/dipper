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

#  Perhaps an improvement would be to to also
# express subject and objects as individual curies (ala predicates)
# iff their namespace is also used by a predicate.

##########################################################
# remove first and last chars from input
function trim(str){ 
	return substr(str,2,length(str)-2)
}

# remove final element in slashed path
# this leaves the namespace of the removed identifier
function stripid(uri){
	l = match(uri,/.*\//)  # side effect is RLENGTH
	if(l<=0)
		return uri
	else 
		return substr(uri,1,RLENGTH-1)
}

# keep underscore, letters & numbers 
# change the rest to (a single) underscore sans leading & trailing
function simplify(str){
	gsub(/[^[:alpha:][:digit:]_]+/,"_",str)
	gsub(/^_+|_+$/, "",str)
	gsub(/__*/,"_",str)
	return str
}

# if possible, find a shorter form for the input
function contract(uri){
	if(uri in prefix)
		return prefix[uri]
	else 
		for(ex in exception){
			start = match(uri, ex)
			if(0 < start)
				return exception[substr(uri, start, RLENGTH)]
		}
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
	# Monarchinitiave's dipper ingest files have cases 
	# which escape the general rules
	# in everything
	exception["http://www.w3.org/2002/07"]="OWL_"
	exception["http://www.w3.org/1999/02"]="RDF"
	exception["http://www.w3.org/2000/01"]="RDFS"
	# in go,mgi,wormbase
	exception["http://dx.doi.org"]="DOI"
	# in mgi
	exception["https://www.mousephenotype.org"]="IMPC_"
	# in panther
	exception["http://identifiers.org/wormbase"]="WormBase"
	# just until skolemized bnodes get in the curie map?
	exception["https://monarchinitiave.org/.well-known/genid"]="BNODE"
	# till all non httpS: are purged
	exception["http://monarchinitiave.org"]="MONARCH"
}

# main loop
# parse and stash the curie file (first file)
(FNR == NR) && /^'[^']*' *: 'http[^']*'.*/ {
	split($0,a,"'")
	prefix[stripid(a[4])]=a[2]
}
# process the ntriple file(s)  which are not the first file
### case when subject predicate and object are all uri
(FNR != NR) && /^<[^>]*> <[^>]*> <[^>]*> \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	q = final(trim($2))
	p = contract(stripid(trim($2)))
	### Object (like subject)
	o = contract(stripid(trim($3)))
	edgelist[s " -> " o " [label=\"" p ":" q]++
}
### case when the object is a literal
(FNR != NR) && /^<[^>]*> <[^>]*> "[^"]*" \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	q = final(trim($2))
	p = contract(stripid(trim($2)))
	### not a uri
	o = "LITERAL"
	edgelist[s " -> " o " [label=\"" p":"q ]++
	nodelist[o " [shape=record];"]++
}

# output dot file, include edge counts
END{
	print "digraph {"
	print "rankdir=LR;"
	print "charset=\"utf-8\";"
	for(edge in edgelist) print edge " (" edgelist[edge] ")\"];"
	for(node in nodelist) print node
	print "labelloc=\"t\";"
	title = final(FILENAME)
	now = systime()
	# return date(1)-style output
	datestamp = strftime("%Y%m%d", now)
	print "label=\"" substr(title,1,length(title)-3) " (" datestamp ")\";"
	print "}"
}

