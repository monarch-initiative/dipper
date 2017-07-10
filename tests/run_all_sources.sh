#!/bin/bash

# don't want the hardcoded sources
# ../dipper.py --sources hpoa,biogrid,impc,panther,ncbigene,\
# ucscbands,ctd,genereviews,eom,zfin,clinvar,coriell,mgi


# about 20 hours sequentially
# might want to run them in parallel
# (i.e. panther first, along with everything else on other cores )

# note this places the 'out' dir in the directory it located in

# don't want the __init__.py file
for source in ../dipper/sources/[^_]*.py; do
    file=${source##*/}
    name=${file%.py}
    src=${name,,}
    # don't want the resident super classes: Sources & PostgresqlSource
    if ! [[ ${src} =~ ".*sources*$" ]]; then
        # echo "../dipper.py --sources ${src}";  # --no_verify
        ../dipper-etl.py --sources "${src}"
    fi
done
