#! /usr/bin/awk -f
# turtle_merge.awk  t_1.ttl t_2.ttl > t_3.ttl
#
# Notes:
# You could merge more than two files
#
# Written to merge some small turtle formatted _METADATA_ files
# that is "small", as in maybe a page long files (i.e. a few kb)

# This script does not know anything about IRI or curies
# just plain ol' string pattern matching.

# If multiple prefixes have different base-uri then the last one wins
#
# Order is neither preserved nor reliably repeated.
#
# There may be room to reduce redundancy in object lists in some cases

BEGIN{RS= " \\.\n*"}

# rdf prefix declarations
/^@prefix /{prefix[$2]=$3}

# rdf record
/^<[^>].*> / {
    # isolate subject
    s = $1; $1 = "";
    # strip newlines
    gsub("\n", "")
    # collapse space
    gsub("  ", " ")
    # partition predicates list
    split(substr($0,2), plist, " ; ")
    for(pos in plist){
        prdobjs = plist[pos]
        p = substr(prdobjs, 1, index(prdobjs, " "))
        prdobjs = substr(prdobjs, length(p)+1)
        # what remains are objects
        record[s,p] = prdobjs
    }
}

END {
    for(x in prefix){print "@prefix " x " " prefix[x] " ."}
    print ""
    # regroup by subjects
    for(r in record){
        split(r, sp, SUBSEP)
        s = sp[1]; p = sp[2]
        subject[s] = subject[s] "\t" p " " record[r] " ;\n"
    }

    # resorting could happen here

    for(s in subject){
        print s
        print substr(subject[s],1, length(subject[s])-3) " .\n"
    }
}



