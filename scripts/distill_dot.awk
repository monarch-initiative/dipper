#! /usr/bin/gawk -f
# distill_dot.awk <numbered.gv>  > <filtered_weighted.gv>
# 

# I have dot files with edge lables that terminate with a number in parens `(nnn)`.
# the numbers may be simple counts (n > 0) or may be differences [nil neg zero pos ]
#
# I'm initally intrested in the case where the numbers are a difference 
# and will filter out small differences and colorize/weight larger changes.
#
# To normalize, will need to read through first then output

BEGIN {hi = -2^PREC; lo = 2^PREC; THRESHOLD = 10}

# there is an unfortunate case where the label contains parens e.g.
# ZFIN -> LITERAL [label=<dc:title ! title (dce) (30783)>];

# AspGD -> GO [label=<RO:0002331 ! involved in (-12)>, color="pink"];


/.*label=<.*\(-?[0-9]+\).*;$/{

    n = split($0, a, /[()]/)
    # for(i in  a) print a[i]
    val = a[length(a)-1]
    sfx = a[length(a)]
    # print "Val=" val
    if(val > hi) hi = val
    if(val < lo) lo = val
    aval = (val < 0.0) ? -val : val
    if(aval > THRESHOLD){
        a[length(a)] = ""
        a[length(a)-1] = ""
        key=""
        for(x in a) key = key a[x]
        value[key] = val
        sufix[key] = substr(sfx,2)
    }
}

# not processed
!/.*label=<.*\(-?[0-9]+\).*;$/ {
    pass_on[pss++] = $0
}

END {
    SCALE =  (hi > -lo)? hi: -lo
    # Threshold could extended as a fraction of scale ...
    THRESHOLD =  SCALE / 200

#    if (SCALE < 1){
#       print "ERROR not to scale: " SCALE
#        SCALE = 1
#    }


    for(i = 0; i < pss; i++)  # not the last line
        if( (pass_on[i]!="") && (pass_on[i]!="}")) print pass_on[i]

    for(key in value)   {
        val = value[key]
        aval = (val < 0.0) ? -val : val
        if(aval >= THRESHOLD)
            print key ">, penwidth = \""(aval/SCALE*5)"\", weight = \"" aval "\"" sufix[key]
    }

    print "}"
}
