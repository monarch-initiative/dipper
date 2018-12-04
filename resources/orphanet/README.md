
## Orphanet XML

To survey:

    xmlstarlet el -a  ../../raw/orphanet/en_product6.xml | sort -u > [orphanet6.xpath](https://github.com/monarch-initiative/dipper/blob/master/resources/orphanet/orphanet6.xpath)

    xpath2dot.awk -v ORIENT="UD"  en_product6.xpath > en_product6.gv

    dot -Tpng en_product6.gv > en_product6.png

    ![schema](https://github.com/monarch-initiative/dipper/blob/master/resources/orphanet/en_product6.png)


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


Orphanet's  ORPHA:nnn  like their internal gene identifiers,
are too sub optimal to propagate
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
