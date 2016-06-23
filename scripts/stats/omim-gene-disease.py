from rdflib.graph import ConjunctiveGraph

graph = ConjunctiveGraph()

graph.parse('http://data.monarchinitiative.org/ttl/omim.ttl', format='turtle')
query = """
            SELECT DISTINCT ?gene ?disease
            WHERE {
                ?anonVariant a OBO:GENO_0000002 ;
                    OBO:GENO_0000408 ?gene .

                {?anonVariant OBO:RO_0002200 ?disease . }
                UNION {?anonVariant OBO:RO_0002326 ?disease . }
                UNION {?anonVariant OBO:RO_0002607 ?disease . }
            }
        """
count = 0
query_result = graph.query(query)
count += len(list(query_result))


query = """
            SELECT DISTINCT ?gene ?disease
            WHERE {
                ?gene rdfs:subClassOf OBO:SO_0000704 .
                ?gene owl:equivalentClass ?eqGene .

                {?eqGene OBO:RO_0002200 ?disease . }
                UNION {?eqGene OBO:RO_0002326 ?disease . }
                UNION {?eqGene OBO:RO_0002607 ?disease . }
            }
        """

query_result = graph.query(query)
count += len(list(query_result))
print("Gene-disease count: {0}".format(count))

query = """
            SELECT DISTINCT ?heritable_marker ?disease
            WHERE {
                ?heritable_marker a OBO:SO_0001500 .

                {?heritable_marker OBO:RO_0002200 ?disease . }
                UNION {?heritable_marker OBO:RO_0002326 ?disease . }
                UNION {?heritable_marker OBO:RO_0002607 ?disease . }
            }
        """
pcount = 0
query_result = graph.query(query)
pcount += len(list(query_result))
print("heritable phenotypic marker-disease count: {0}".format(pcount))
