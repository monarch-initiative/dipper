#! /usr/bin/awk -f
# distill_dot.awk <numbered.gv>  > <filtered_weighted.gv>
# 

# I have dot files with edge lables that terminate with a number in parens `(nnn)`.
# the numbers may be simple counts (n > 0) or may be differences [nil neg zero pos ] 
# 
# I'm initally intrested in the case where the numbers are a difference 
# and will filter out small differences and colorize/weight larger changes.    


BEGIN: {hi=-1000000; lo=1000000}

# there is an unfortunate case where the label cantains parens e.g. 
# ZFIN -> LITERAL [label=<dc:title ! title (dce) (30783)>];
# will need to 

/^[^ ]+ -> [^ ]+ \[label=<.*([^)]+>.*];$/ {

    # 
}
