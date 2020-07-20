##### Infer biolink categories from monarch.owl and rdf:type triples

A robot -> arq pipeline that construct biolink:category triples from 
either monarch.owl and the biolink-model.ttl file alone, or along with
rdf:type triples in dipper

For example, given the input turtle

```
@prefix OMIM: <https://omim.org/entry/> .
@prefix HP: <http://purl.obolibrary.org/obo/HP_> .
@prefix SO: <http://purl.obolibrary.org/obo/SO_> .

<http://www.example.com/GeneA> a SO:0000704 .
<http://www.example.com/GeneA> <http://www.example.com/someRelation> OMIM:154700 .
OMIM:154700 <http://www.example.com/hasPhenotype> HP:0002705 .
```

We expect the output:

```
@prefix owl:   <http://www.w3.org/2002/07/owl#> .
@prefix meta:  <https://w3id.org/biolink/biolinkml/meta/> .
@prefix rdfs:  <http://www.w3.org/2000/01/rdf-schema#> .
@prefix biolink: <https://w3id.org/biolink/vocab/> .

<https://omim.org/entry/154700>
        biolink:category  biolink:Disease .

<http://purl.obolibrary.org/obo/SO_0000704>
        biolink:category  biolink:GenomicEntity , biolink:Gene .

<http://www.example.com/GeneA>
        biolink:category  biolink:GenomicEntity , biolink:Gene .

<http://purl.obolibrary.org/obo/HP_0002705>
        biolink:category  biolink:PhenotypicFeature .
```

make requires 60G of memory (for ClinVar), in the future these queries may
be chunked to lessen the memory requirement, although robot reason is also
memory intensive prohibiting this from running on most laptops.

On monarch4 this runs in 40 minutes running make in parallel mode via
```make -j 4```