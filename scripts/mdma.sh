#! /bin/bash

set -e
set -u
set -o pipefail

# mdma.sh
# Move Data to Monarch Archive
# to be run (estaticly) by jenkins ci
#    DTSTMP=$(date +%Y%m%d)  # higher resolution date stamp
    YYYYMM="${YYYYMM:-$(date +%Y%m)}"
    MONARCHIVE="monarch@monarch-archive:/var/www/data/$YYYYMM"

    # copy main turtle files to ntriple format
    # https://github.com/drobilla/serd
    mkdir -p ntriples

    # avoid having redundant ntriple files in both ntriples/ and rdf/ dirs
    for t in  ./out/*.nt ; do
        sort -u "$t" > "ntriples/${t##out/}" && rm -f "$t"
    done

    # enable extended gobbing for selecting filenames without an underbar
    shopt -s extglob
    stat --printf="%s\t%n\n" ./out/+([^_]).ttl |sort -nr |cut -f2 |
        parallel -j0 "serdi -i turtle -o ntriples {} | sort -u > ntriples/{/.}.nt"

    # keeping the new .nt files group writable (same as the .ttl are) helps downstream
    chmod g+rw ./ntriples/*.nt

    # Scigraph requires the RDF data files be self-declared as OWL ontologies
    # p="<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>"
    # o="<http://www.w3.org/2002/07/owl#Ontology>"
    # (
    #     cd ntriples
    #     for nt in *.nt; do
    #         echo "<https://archive.monarchinitiative.org/#${nt%%.nt}> $p $o ." >> "$nt"
    #     done
    # )

    # compress the ntriples as machine readable shouldn't care (& humans can use zcat)
    tar -czf ntriples_"$YYYYMM.tgz" ntriples/
    # rename; as 'out' makes no sense downstream. todo: change in dipper
    mv -u  out rdf
    # rdf/ dir cannot be the complete rekease without the .nt files so re adding them
    mv -u ntriples/*.nt rdf/
    rm -fr ntriples/
    chmod -R g+rw ./rdf
    tar -czf "rdf_$YYYYMM.tgz" rdf

    # copy the yaml files along with the data they were applied to
    rsync dipper/curie_map.yaml "$MONARCHIVE"/translationtable/
    rsync translationtable/GLOBAL_TERMS.yaml "$MONARCHIVE"/translationtable/
    # copy the compressed  files to archive
    rsync "rdf_$YYYYMM.tgz" "$MONARCHIVE"/
    rsync "ntriples_$YYYYMM.tgz" "$MONARCHIVE"/

    # move the dir back for now
    mv -u rdf out

    # need a higher resolution than year-month in this case
    # hold off  for now.
    # git tag -a $DTSTMP -m "release candidate $DTSTMP"
    # git push origin --tag
    # Kent is suggesting adding commit hash to metadata which has the same effect
    # without creating tag noise
