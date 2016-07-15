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

function contract(uri){
	if(uri in prefix)
		return prefix[uri]
	else {
		gsub(/[:./-]+/,"_",uri)
		return uri
	}
}

# parse and stash the curie file (first file)
(FNR == NR) && /^'[^']*' *: 'http[^']*'.*/ {
	split($0,a,"'")
	prefix[stripid(a[4])]=a[2]
}

# process the ntriple file(s)  which are not the first file
### all are uri
(FNR != NR) && /^<[^>]*> <[^>]*> <[^>]*> \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	split(trim($2),b,"/")
	p = b[length(b)]
	### Object (like subject)
	o = contract(stripid(trim($3)))
	edgelist[s " -> " o " [label=\"" p "\"];"]++
}
### object is a literal
(FNR != NR) && /^<[^>]*> <[^>]*> "[^"]*" \.$/ {
	### Subject (uri)
	s = contract(stripid(trim($1)))
	### Predicate (uri)
	split(trim($2),b,"/")
	p = b[length(b)]
	### not a uri
	o = "LITERAL"
	edgelist[s " -> " o " [label=\"" p "\"];"]++
	nodelist[o " [shape=record];"]++
}

END{
	print "digraph {"
	print "rankdir=LR;"
	print "charset=\"utf-8\";"
	for(edge in edgelist){
		print edge
	}
	for(node in nodelist){
		print node
	}
	print "}"
}

