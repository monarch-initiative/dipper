#!/bin/bash
# generate source names as a yaml map
# which is a bit dodgy since they can be generated anyway
for source in dipper/sources/[^_]*.py; do
    file=${source##*/}
    name=${file%.py}
    src=${name,,}
    # don't want the resident super classes: Sources & PostgresqlSource
    if ! [[ ${src} =~ 'sources' ||  ${src} =~ 'postgresqlsource' ]] ; then
        echo -e "\t'${src}':\t'${name}',"
    fi
done

