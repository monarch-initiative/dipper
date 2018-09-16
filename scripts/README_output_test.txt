
for ttl in  /var/lib/jenkins/workspace/build-*ttl/out/[a-z]*.ttl; do 
	o=${ttl##*/};rapper -iturtle  -o ntriples ${ttl} > out/${o%%.ttl}.nt; 
done

for nt in /var/lib/jenkins/workspace/build-*ttl/out/[a-z]*.nt; do 
	cp ${nt} out/${nt##*/}; 
done


ls -lat /var/lib/jenkins/workspace/build-clinvar-ttl/out/clinvar.ttl
-rw-r--r-- 1 jenkins jenkins 300830534 Jun 28 09:08 
/var/lib/jenkins/workspace/build-clinvar-ttl/out/clinvar.ttl

nano  +254506 /var/lib/jenkins/workspace/build-clinvar-ttl/out/clinvar.ttl

:MONARCH_033744041f63e434690733d5cfcc51c5 a OBAN:association ;
    OBAN:association_has_object <http://purl.obolibrary.org/obo/OMIM_^@> ;
    OBAN:association_has_object_property OBO:RO_0002200 ;
    OBAN:association_has_subject <http://www.ncbi.nlm.nih.gov/clinvar/variation/215826> .


################################################################################

for nt in out/*.nt; do scripts/ntriple2dot.awk dipper/curie_map.yaml ${nt} > ${nt%.nt}.dot ; done
# four hours later

for dot in out/*.dot; do dot -Tpng ${dot} > ${dot%.dot}.png ; done
Warning: out/ctd.dot: syntax error in line 4 near '-'
Warning: out/impc.dot: syntax error in line 27 near '-'
Warning: out/kegg.dot: syntax error in line 4 near '-'
Warning: out/mpd.dot: syntax error in line 22 near '-'
Warning: out/omia.dot: syntax error in line 17 near '-'



time for nt in out/*.nt;do if [ "dataset.nt" == ${nt##*_} ] ; then continue ;fi; echo -e ${nt} "\t" `date`;scripts/ntriple2dot.awk dipper/curie_map.yaml ${nt} > ${nt%.nt}.dot ; done



time for nt in out/*.nt;do if [ "dataset.nt" == ${nt##*_} ] ; then continue ;fi; echo ${nt} " " `date` ;nohup scripts/ntriple2dot.awk dipper/curie_map.yaml ${nt} > ${nt%.nt}.dot ; done


for nt in out/*.nt;do if [ "dataset.nt" == ${nt##*_} ] ; then continue ;fi; echo ${nt}; done > nt.list


date
nohup cat nt.list | parallel -j11 'scripts/ntriple2dot.awk dipper/curie_map.yaml {} > {.}.dot'
date

