#!/bin/bash

# ../dipper.py --sources hpoa,biogrid,impc,panther,ncbigene,\
# ucscbands,ctd,genereviews,eom,zfin,clinvar,coriell,mgi

# dont want the hardcoded sources

# don't want the random __init__.py file
for source in ../dipper/sources/[^_]*.py; do
    file=${source##*/}
    name=${file%.py}
    src=${name,,}
    # don't want the resident super classes: Sources & PostgresqlSource
    if ! [[ ${src} =~ ".*sources$" ]]; then
        # echo "../dipper.py --sources ${src}";  # --no_verify
        ../dipper.py --sources ${src}"
    fi
done
