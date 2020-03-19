#! /bin/bash -e

# mdma.sh
# Move Data to Monarch Archive
# to be run (estaticly) by jenkins ci
    DTSTMP=$(date +%Y%m%d)
    YYYYMM="${YYYYMM:-$(date +%Y%m)}"
    MONARCHIVE="monarch@monarch-archive:/var/www/data/$YYYYMM"

    # copy main turtle files to ntriple format
    # https://github.com/drobilla/serd
    mkdir -p ntriples
    if [ -e  out/*.nt ] ; then
        mv -u out/*.nt ntriples
    fi

    # enable extended gobbing for selecting filenames without an underbar
    shopt -s extglob
    stat --printf="%s\t%n\n" ./out/+([^_]).ttl |sort -nr |cut -f2 |
        parallel -j0 "serdi -i turtle -o ntriples {} > ntriples/{/.}.nt"

    # compress the ntriples as machine readable shouldn't care (humans use zcat)
    tar -czf ntriples_$YYYYMM.tgz ntriples/
    # rename; as 'out' makes no sense downstream. todo: change in dipper
    mv -u  out rdf
    mv -u ntriples_$YYYYMM.tgz rdf/
    rm -fr ntriples/
    tar -czf rdf_$YYYYMM.tgz rdf

    # move the yamal files along with the data they were applied to
    scp dipper/curie_map.yaml $MONARCHIVE/translationtable/
    scp translationtable/GLOBAL_TERMS.yaml $MONARCHIVE/translationtable/
    scp rdf_$YYYYMM.tgz $MONARCHIVE/
    scp -r ./rdf/* $MONARCHIVE/rdf/
    # move the dir back for now
    mv -u rdf out

    # need a higher resolution than year-month in this case
    git tag -a $DTSTMP -m "release candidate $DTSTMP"
    git push origin --tag

