#! /usr/bin/gawk -f

#  Reduce the subject and object of RDF triples (ntriples format)
#  down to their curie prefix (or literal object)
#  Reduce predicates to the specific term identifier
#  and the ontology namespace they are from. (curi prefix )
#  Express the subject and object as nodes with a directed edge
#  labeled with the predicate in the graphviz dot format
#  Include a tally of each combination of nodes and edge

#  This may over generalize in some cases because I do not have
#  a handy way to differentiate uri for subjects and objects
#  which may belong to a "structural" ontology (t-box) as opposed to
#  a data uri (a-box).

#  Post processing to augment predicate identifiers
#  with their labels seems to improve usefulness.

function usage(){
    print "usage: ntriple2dot.awk prefix_baseurl.yaml rdf.nt > dot.gv"
}

##########################################################
# remove first and last chars from input string
function trim(str){
    return substr(str,2,length(str)-2)
}

# break URI into tokens to use a a path in
# multi-dim array pointing to a curie-prefix
# base URI include "!#&-=?_"
# which are all also present in localID as well
function tokenize (str){
    # b/c found in local IDs leave .:
    gsub(/[!&:#=?/_-]+/, SUBSEP, str)
    # discard the first token it is only low info content  http|https|ftp
    # discard any empty trailing tokens
    tail = match(str, ".*" SUBSEP)  # side effects RSTART & RLENGTH
    # print "tail = " tail " tail length = " RLENGTH
    finish = 0
    start = index(str, SUBSEP) + 1
    if (tail > 0)
        finish = RLENGTH - start
    else
        finish = length(str) - start + 1
    return substr(str, start, finish)
}

# remove final token from token path delimited with SUBSEP
# this ends up matching chars included in localID
# ._-;():%#,?=@*\&+   (yea ... nope)

function detail(key,   start) {
    start = match(key, ".+" SUBSEP)
    if (start == 1)
        return substr(key, 1, RLENGTH -1)
    else
        return key
}

# for output, keep underscore, letters & numbers
# change the rest to (a single) underscore sans leading & trailing
# this passes valid node labels from dot's perspective (C-lang identifier rules)
function simplify(str){
    gsub(/[^[:alpha:][:digit:]_]+/,"_",str)
    gsub(/^_*|_*$/, "",str)
    gsub(/__*/,"_",str)
    return str
}

# Replace IRI with CuriePrefix when possible
# Otherwise whinge and return a gv printable version of the original
function contract(uri,  u){
    u = tokenize(uri)

    # sometimes uri is already correct as is
    # if (index(u, SUBSEP) > 0)
    #    u = detail(u)
    # if (u == "") {
    #    printf("ERROR contracting %s\n", uri) > "/dev/stderr"
    #    printf("ERROR tokenizeing %s\n", tokenize(uri)) > "/dev/stderr"
    #    printf("ERROR detailing   %s\n", detail(tokenize(uri))) > "/dev/stderr"
    #    # exit(1)  # while testing
    # }

    # shorten till longest uri in curi map is found (or not)
    while(!(u in prefix) && index(u, SUBSEP)>0)
        u = detail(u)

    if(u in prefix)
        return prefix[u]
    else{
        printf("WARN  %s\n", uri) > "/dev/stderr"
        # printf("BASE  %s\n", detail(tokenize(uri))) > "/dev/stderr"
        return simplify(uri)
    }
}

# get the final (incl fragment identifier) portion of a slashed path
function final(uri,  b, p, anchor){
    split(uri, b, "/")
    p = b[length(b)]
    anchor = match(p, "#")
    if(anchor > 0)
        p = substr(p, anchor+1)
    return p
}

# fix curie namespace when redundant due to OBO  i.e RO:RO_123
# is this still needed circa 2020?
function deoboify(curie,  a){
    split(curie, a , ":")
    if(match(a[2], a[1]"_") > 0){
        sub("_", ":", a[2]);
        curie = a[2]
    }
    return curie
}

BEGIN{
    # exceptions
    prefix["BNODE"]="BNODE"  # is a fixed point
    prefix[tokenize("https://monarchinitiative.org/.well-known/genid")]="BNODE"
    prefix[tokenize("https://archive.monarchinitiative.org/")]="MonarchArchive"
    ############################################################
    # not expected to be added to curie map
    # in HPOA
    prefix[tokenize("http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=gene&part=")]="NCBIBook"
    prefix[tokenize("http://www.ncbi.nlm.nih.gov/books/NBK")]="NCBINBK"

    ############################################################
    # Often re-visit whether these exceptions are still necessary
    # they may have been moved into the curie prefix mapping file
    # or addressed in a local translation table
    # https://raw.githubusercontent.com/monarch-initiative/dipper/master/dipper/curie_map.yaml

    # in mgi
    prefix[tokenize("https://www.mousephenotype.org")]="IMPC"
    # in IMPC (not httpS  --- hmm all three IRI exist)
    prefix[tokenize("http://www.mousephenotype.org")]="IMPC"
    # note in curie_map IMPC is:
    # http://www.mousephenotype.org/data/genes/

    # in panther
    # prefix["http://identifiers.org/wormbase/"]="WormBase"
    # note in curie_map WormBase is
    # 'https://www.wormbase.org/get?name='

    # various, transient?
    prefix[tokenize("http://www.ncbi.nlm.nih.gov/gene/")]="NCBIGene"
    prefix[tokenize("http://www.ncbi.nlm.nih.gov/genome/")]="NCBIGenome"
    prefix[tokenize("http://www.ncbi.nlm.nih.gov/assembly?term=")]="NCBIAssembly"
    # kegg
    prefix[tokenize("http://www.genome.jp/kegg/pathway/map/")]="KEGG-img"
    # in EOM
    prefix[tokenize("https://elementsofmorphology.nih.gov/images/terms/")]="EOM_IMG"
    prefix[tokenize("http://elementsofmorphology.nih.gov/index.cgi?tid=")]="EOM"  # w/o httpS

    # playing with the idea of a LOCAL identifier
    # not the same as a bnode in that tools won't rewrite identifier
    # and their IRI are unroutable, i.e in same machine (or local net)
    # this is motivated in part by not wanting to publish
    # urls to our own servers which go nowhere.
    # prefix["https://127.0.0.1/.well-known/genid"]="BNODE"  # exists, phase out
    # prefix["http://0.0.0.0/.well-known/genid"]="BNODE"     # prefer
    # prefix["http://0.0.0.0/"]="LCL"                              # ID halfway house

    # LCL's obvious down side is zero options for global collisions
    # best you can do is go the uuid which is ugly, or
    # belive your obscure but readable string _is_ a special snowflake
}

# main loop
# parse and stash the curie yaml file (first file)
# YAML format is tic delimited word (the curi prefix)
# a colon, whitespace
# then a tic delimited url
# convert the base url to a path in a multi-dimensional array
# (FNR == NR) && /^'[fht]\+tp[^']*'.*/ { # loosing some?
(FNR == NR) && /^'.*/ {
    split($0, arr, "'")
    if(arr[2]=="")
        arr[2]="BASE"
    prefix[tokenize(arr[4])]=arr[2]
}

# process the ntriple file(s)  which are not the first file

### subject predicate and object are all uri
(FNR != NR) && /^<[^>]*> <[^>]+> <[^>]+> \.$/ {
    s = contract(trim($1))
    p =  final(trim($2))
    ns = contract(trim($2))
    o = contract(trim($3))
    edgelist[s,ns ":" p,o]++
    done[NR]=1
}
### subject & predicate are uri but the object is a literal
(FNR != NR) && /^<[^>]*> <[^>]*> "[^"]*".*\.$/ {
    s = contract(trim($1))
    p = final(trim($2))
    ns = contract(trim($2))
    edgelist[s, ns ":" p, "LITERAL"]++
    done[NR]=1
}
### subject is a bare blank node, predicate and object are uri
(FNR != NR) && /^_:[^ ]* <[^>]*> <[^>]*> \.$/ {
    # once our bnode syntax is uniform we may want more info out of it
    # s = contract(($1)
    p =  final(trim($2))
    ns = contract(trim($2))
    ### Object (like subject)
    o = contract(trim($3))
    edgelist["BNODE",ns ":" p,o]++
    done[NR]=1
}
### subject & predicate are uri & object is a bare blank node
(FNR != NR) && /^<[^>]*> <[^>]*> _:[^ ]* \.$/ {
    s = contract(trim($1))
    p =  final(trim($2))
    ns = contract(trim($2))
    #o = contract($3))
    edgelist[s,ns ":" p,"BNODE"]++
    done[NR]=1
}
### subject and object are a bare blank nodes, predicate is a uri
(FNR != NR) && /^_:[^ ]* <[^>]*> _:[^ ]* \.$/ {
    # s = contract($1)
    p =  final(trim($2))
    ns = contract(trim($2))
    # o = contract($3)
    edgelist["BNODE",ns ":" p,"BNODE"]++
    done[NR]=1
}
### subject is a bare blank node, predicate is a uri & the object is a literal
(FNR != NR) && /^_:[^ ]* <[^>]*> "[^"]*".*\.$/ {
    # s = contract($1)
    p = final(trim($2))
    ns = contract(trim($2))
    edgelist["BNODE", ns ":" p, "LITERAL"]++
    done[NR]=1
}

# print anything else so it makes itself known
(FNR != NR) && (done[NR]<1) && (NF>1){
    printf("ERROR? in %s at line %i:  %s\n",FILENAME,FNR,$0) > "/dev/stderr"
}

# output graphviz dot file, include edge counts
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
        " [label=\"" deoboify(spo[2]) " (" edgelist[edge] ")\"];"
    }
    print "LITERAL [shape=record];"
    print "labelloc=\"t\";"
    title = final(FILENAME)
    datestamp = strftime("%Y%m%d", systime())
    print "label=\"" substr(title,1,length(title)-3) " (" datestamp ")\";"
    print "}"
}
