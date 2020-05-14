
## Orphanet XML

To (re)survey:

PREVREL=201811
RELEASE=202005

    xmlstarlet el -a  ../../raw/orphanet/en_product6.xml | sort -u > en_product6_$RELEASE.xpath
    xpath2dot.awk -v ORIENT="UD"  en_product6_$RELEASE.xpath > en_product6_$RELEASE.gv

    dot -Tpng en_product6_$RELEASE.gv > en_product6_$RELEASE.png
    dot -Tsvg en_product6_$RELEASE.gv > en_product6_$RELEASE.svg

    ![schema](https://github.com/monarch-initiative/dipper/blob/master/resources/orphanet/en_product6_$RELEASE.png)


note: 2018 Fall, there is no `GeneList` element in the XML
although it is refered to in the code and exists in the dipper/test/resources/orphanet
examples.

    grep GeneList en_product6.xpath

with that change the GeneType@id is also gone?

moved? may now be:

    JDBOR/DisorderList/Disorder/DisorderGeneAssociationList
        /DisorderGeneAssociation/DisorderGeneAssociationType/@id

but those are

```
xmlstarlet sel -t -v './JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/DisorderGeneAssociation/DisorderGeneAssociationType/Name'   en_product6.xml | sort | uniq -c | sort -nr
   4542 Disease-causing germline mutation(s) in
   1084 Disease-causing germline mutation(s) (loss of function) in
    495 Major susceptibility factor in
    303 Candidate gene tested in
    228 Part of a fusion gene in
    226 Role in the phenotype of
    206 Disease-causing germline mutation(s) (gain of function) in
    177 Disease-causing somatic mutation(s) in
     91 Biomarker tested in
     38 Modifying germline mutation in

xmlstarlet sel -t -v './JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/DisorderGeneAssociation/DisorderGeneAssociationStatus/Name'   en_product6.xml | sort | uniq -c | sort -nr
   7006 Assessed
    384 Not yet assessed

```

------------------------------------------------------------------------------

there is a "source of validation" which if it exists is,
a list of PMIDs (underscore separated)  i.e.
    `<SourceOfValidation>27018472[PMID]_27018475[PMID]</SourceOfValidation>`

```
fgrep  '<SourceOfValidation' en_product6.xml | fgrep "[" | wc -l
6560
fgrep  '<SourceOfValidation' en_product6.xml | fgrep "[PMID" | wc -l
6560
```

 all source of validations have at least the first as a PMID

-------------------------------------------------------------------


Orphanet's  ORPHA:nnn  like their internal gene identifiers, are too sub optimal to propagate
(not even obvious how to search for non-disease-orpha-numbers on their site)

Idealy human gene's will be HGNC but no source covers all their genes


```
 grep -E '<Source>|Gene id' en_product6.xml | cut -f1 --d '"' |sort | uniq -c | sort -nr
   7390           <Gene id=
   7387                 <Source>HGNC</Source>
   7343                 <Source>Ensembl</Source>
   7337                 <Source>SwissProt</Source>
   7329                 <Source>OMIM</Source>
   7316                 <Source>Genatlas</Source>
   6224                 <Source>Reactome</Source>
   1047                 <Source>IUPHAR</Source>

```

I will take the first stable curie by rank (as HGNC,Ensembl,SwissProt,OMIM...)
as there are only a few handfuls out of 7k, this is mostly moot.

ahh nope. default to the gene's orphanet local identifier
some may not have any external id's

 use something like this to (re) generate files for the tests

```
xmlstarlet sel -t -c ./JDBOR/DisorderList/Disorder[@id="17604"]  en_product6.xml
<Disorder id="17604">
      <OrphaNumber>166035</OrphaNumber>
      <Name lang="en">Brachydactyly-short stature-retinitis pigmentosa syndrome</Name>
      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <SourceOfValidation>28285769[PMID]</SourceOfValidation>
          <Gene id="26792">
            <OrphaNumber>510666</OrphaNumber>
            <Name lang="en">CWC27 spliceosome associated protein homolog</Name>
            <Symbol>CWC27</Symbol>
            <SynonymList count="2">
              <Synonym lang="en">NY-CO-10</Synonym>
              <Synonym lang="en">SDCCAG-10</Synonym>
            </SynonymList>
            <ExternalReferenceList count="0">
            </ExternalReferenceList>
            <LocusList count="0">
            </LocusList>
          </Gene>
          <DisorderGeneAssociationType id="17949">
            <Name lang="en">Disease-causing germline mutation(s) in</Name>
          </DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus id="17991">
            <Name lang="en">Assessed</Name>
          </DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>

```

My plan at the moment is to drop the blank node
and instead make the subject of assocciations be a gene (hgnc when possible)
and make any other dbxrefs equivilent

a WIP model of intent (incomplete)


################################################################
# 2020 May

Their model is changed and expanded

To see just how, I tried deltadot but it is to specific to our dipper output.
So here I will look at it long hand.

comm -23 en_product6_$PREVREL.xpath en_product6_$RELEASE.xpath > xp_dropped
comm -13 en_product6_$PREVREL.xpath en_product6_$RELEASE.xpath > xp_added
comm -12 en_product6_$PREVREL.xpath en_product6_$RELEASE.xpath > xp_continue 


the path

JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/DisorderGeneAssociation/Gene/OrphaNumber

is dropped  which is fine, see previous discussion for why.

An OrphaNumber field is now attached as:  JDBOR/DisorderList/Disorder/OrphaNumber
which may be more appropiate, especially if it resolves at a relevent page there.

16 paths are added, but some may just be moved ... no it does not seem that way
as only  `Name` `lang` pairs are conserved.

```
echo "digraph G {" > orpha_update_$RELEASE.gv
awk -F'/' '{print "\t" $(NF-1) " -> " $(NF) " [color=\"red\"];"}' <(tr -d '@' <xp_dropped) >> orpha_update_$RELEASE.gv
awk -F'/' '{print "\t" $(NF-1) " -> " $(NF) " [color=\"blue\"];"}' <(tr -d '@' <xp_added) >> orpha_update_$RELEASE.gv
awk -F'/' '{print "\t" $(NF-1) " -> " $(NF) " [color=\"black\"];"}' <(tr -d '@' <xp_continue) >> orpha_update_$RELEASE.gv
echo "}" >> orpha_update_$RELEASE.gv

```
xdot orpha_update_202005.gv

not pretty, but it is helpful.

added edges are 


	Disorder	DisorderGroup (Name & @id)
	Disorder	DisorderType (Name & id)
	Disorder	ExpertLink (@lang)
	Gene	    GeneType (Name & @id)

	Locus	    GeneLocus
	Locus	    LocusKey


ExpertLink:
The expert fields are all unique URI and look like:
``` 
http://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&amp;Expert=99880
```
but  resolve to an error page saying:  

```exec cExpIdData2 , 'en' SQL query failed```


.../DisorderGroup/Name:
   3163 Disorder
    670 Subtype of disorder
      6 Group of disorders


.../DisorderType/Name:
   2308 Disease
    745 Malformation syndrome
    506 Clinical subtype
    141 Etiological subtype
     92 Morphological anomaly
     23 Histopathological subtype
     10 Biological anomaly
      6 Clinical syndrome
      5 Clinical group
      2 Particular clinical situation in a disease or syndrome
      1 Category

.../Gene/GeneType/Name:
   7678 gene with protein product
     84 Non-coding RNA
     31 Disorder-associated locus

.../Gene/LocusList/Locus/GeneLocus   
    these are Cytogenic locations (think we should get these)

/Gene/LocusList/Locus/LocusKey:
  all but five are the number one and the rest are the number two.
  no idea what it represents.
