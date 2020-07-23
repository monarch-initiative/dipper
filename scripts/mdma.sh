#! /bin/bash -e

# mdma.sh
# Move Data to Monarch Archive
# to be run (estaticly) by jenkins ci
#    DTSTMP=$(date +%Y%m%d)  # higher resolution date stamp
    YYYYMM="${YYYYMM:-$(date +%Y%m)}"
    MONARCHIVE="monarch@monarch-archive:/var/www/data/$YYYYMM"

    # copy main turtle files to ntriple format
    # https://github.com/drobilla/serd
    mkdir -p ntriples
    # avoid redundant ntriple files in botb ntriples/ and rdf/ dirs
    #if [ -e  out/*.nt ] ; then
    mv -u out/*.nt ntriples/
    #fi

    # enable extended gobbing for selecting filenames without an underbar
    shopt -s extglob
    stat --printf="%s\t%n\n" ./out/+([^_]).ttl |sort -nr |cut -f2 |
        parallel -j0 "serdi -i turtle -o ntriples {} > ntriples/{/.}.nt"

    # Scigraph requires the RDF data files be self-delared as OWL ontologies
    (
        cd ntriples
        for nt in *.nt; do
            echo "<https://archive.monarchinitiative.org/#${nt%%.nt}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#Ontology> ." >> "$nt"
        done
    )
    # compress the ntriples as machine readable shouldn't care (humans use zcat)
    tar -czf ntriples_"$YYYYMM.tgz" ntriples/
    # rename; as 'out' makes no sense downstream. todo: change in dipper
    mv -u  out rdf
    mv -u ntriples/*.nt rdf/
    rm -fr ntriples/
    tar -czf "rdf_$YYYYMM.tgz" rdf

    # copy the yaml files along with the data they were applied to
    scp dipper/curie_map.yaml "$MONARCHIVE"/translationtable/
    scp translationtable/GLOBAL_TERMS.yaml "$MONARCHIVE"/translationtable/
    # copy the compressed  files to archive
    scp "rdf_$YYYYMM.tgz" "$MONARCHIVE"/
    scp "ntriples_$YYYYMM.tgz" "$MONARCHIVE"/

    # move the uncompressed files over ...or maybe just decompress a copy there
    # scp -r ./rdf/* $MONARCHIVE/rdf/

    # move the dir back for now
    mv -u rdf out

    # need a higher resolution than year-month in this case
    # hold off  for now.
    # git tag -a $DTSTMP -m "release candidate $DTSTMP"
    # git push origin --tag
    # Kent is suggesting adding commit hash to metadata which has the same effect
    # without creating tag noise
