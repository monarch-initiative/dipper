#! /usr/bin/gawk -f
# spote_count.awk
# RDF ntriple counts for:
# subject predicate object triples entity
#
# spote_count.awk -v"RELEASE=YYYYMM" path/file.nt


function sum(a, t){
    for(i in a)
        t += a[i]
    return t
}

BEGIN{
    if(! RELEASE)
        "date +%Y%m" | getline RELEASE;
    FS = " "
    archive = "https://archive.monarchinitiative.org/"
}

# iri objects
/^.*> .$/ {
    s[$1]++;p[$2]++;o[$3]++
}

# entities  (typed subjects)
$2 == "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>" {
    e[$1]++
}

# literal objects
/^.*" .$/ {
    s[$1]++;p[$2]++;
    # not counting literals for now
    # $1="";$2="";l[$0]++
}

# not counting bnodes for now
# blank node subject
#/^[^ ]*\.well-known\/genid/ {
#    bs[$1]++
#    b[$1]++
#}

# blank node object
#/ <[^ ]*\.well-known\/genid/ {
#    bo[$3]++
#    b[$3]++
#}

END{
    split(FILENAME, path, "/")
    graphname = path[length(path)]
    # still has extension '.nt' might need to truncate

    subject_iri = "<" archive RELEASE "/rdf/" graphname ">"
    subject_curie = "MonarchArchive:" RELEASE "/rdf/" graphname

    # turtle like output for now b/c _datasets output is always turtle
    # not sure if or how this will get intregrated.

    print "@prefix MonarchArchive: <https://archive.monarchinitiative.org/> ."
    print "@prefix void: <http://rdfs.org/ns/void#> ."

    print subject_curie " a dctypes:Dataset ;"  # redundant typing
    print "\tvoid:distinctSubjects " length(s)  " / " sum(s) " ;"
    print "\tvoid:distinctObjects "  length(o)  " / " sum(o) " ;"
    print "\tvoid:entities "         length(e)  " / " sum(e) " ;"
    print "\tvoid:properties "       length(p)  " / " sum(p) " ;"

    # not aware of what predicates to use. (have not researched)
    # print "\texample:literals "      length(l)  " / " sum(l) " ;"
    # print "\texample:blanknode "     length(b)  " / " sum(b) " ;"
    # print "\texample:bnodesub "      length(bs) " / " sum(bs) " ;"
    # print "\texample:bnodeobj "      length(bo) " / " sum(bo) " ;"

    print "\tvoid:triples " FNR " ."

}
