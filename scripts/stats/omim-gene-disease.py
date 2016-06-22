from rdflib.graph import ConjunctiveGraph

graph = ConjunctiveGraph()

graph.parse('http://data.monarchinitiative.org/ttl/omim.ttl', format='turtle')

count = 0

query = """
            SELECT DISTINCT ?gene ?disease
            WHERE {
                ?gene rdfs:subClassOf OBO:SO_0000704 .
                ?gene owl:equivalentClass ?eqGene .

                ?eqGene OBO:RO_0002200 ?disease .
            }
        """

query_result = graph.query(query)
#for i in query_result:
#    print("%s\t%s" % i)
print(list(query_result))
