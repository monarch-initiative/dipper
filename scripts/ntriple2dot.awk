#! /usr/bin/awk -f

# remove first and last chars
function trim(str){ 
	return substr(str,2,length(str)-2)
}

# remove final element in slashed path
function stripid(uri){
	l = match(uri,/.*\//)  # side effect is RLENGTH
	if(l<=0)
		return uri
	else 
		return substr(uri,1,RLENGTH -1)
}

# keep underscore, letters & numbers 
# change the rest to (a single) underscore sans leading & trailing
function simplify(str){
	gsub(/[^[:alpha:][:digit:]_]+/,"_",str)
	gsub(/^_+|_+$/, "",str)
	gsub(/__*/,"_",str)
	return str
}

function contract(uri){
	if(uri in prefix)
		return prefix[uri]
	else 
		for(ex in exception){
			start = match(uri, ex)
			if(0 < start)
				return exception[substr(uri, start, RLENGTH)]
		}
	# else
		return simplify(uri)
}

# get the final (identifier) portion of a slashed path
function final(uri){
	split(trim(uri),b,"/")
	p = b[length(b)]
	anchor = match(p, "#")
	if(anchor > 0)
		p = substr(p, anchor+1)
	return
}

BEGIN{
	# Monarchinitiave's dipper ingest files have cases 
	# which escape the general rules
	# in everything
	exception["http://www.w3.org/2002/07"]="OWL_"
	# in go,mgi,wormbase
	exception["http://dx.doi.org"]="DOI"
	# in mgi
	exception["https://www.mousephenotype.org"]="IMPC_"
	# in panther
	exception["http://identifiers.org/wormbase"]="WormBase"
	# just until skolemized bnodes get in the curie map?
	exception["https://monarchinitiave.org/.well-known/genid"]="BNODE"
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
	p = final($2)
	### Object (like subject)
	o = contract(stripid(trim($3)))
	edgelist[s " -> " o " [label=\"" p ]++
}
### case when the object is a literal
(FNR != NR) && /^<[^>]*> <[^>]*> "[^"]*" \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	p = final($2)
	### not a uri
	o = "LITERAL"
	edgelist[s " -> " o " [label=\"" p ]++
	nodelist[o " [shape=record];"]++
}

END{
	print "digraph {"
	print "rankdir=LR;"
	print "charset=\"utf-8\";"
	for(edge in edgelist){
		print edge " (" edgelist[edge] ")\"];"
	}
	for(node in nodelist){
		print node
	}
	print "}"
}

