from dipper.graph.RDFGraph import RDFGraph
from rdflib.namespace import RDFS
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map
from rdflib import URIRef
import copy
import json
from collections import OrderedDict

def main():

    hpo = RDFGraph()
    root = "HP:0000118"
    hpo_terms = OrderedDict()

    hpo.parse("http://purl.obolibrary.org/obo/hp.owl", format='xml')
    hpo.bind_all_namespaces()
    hpo.bind("oboInOwl", "http://www.geneontology.org/formats/oboInOwl#")

    tree = {}
    tree[root] = {}
    path = []

    hpo_to_tree(root, hpo_terms, hpo, tree, path)

    with open('hpo-tree.json', 'w') as outfile:
        json.dump(tree, outfile)

    with open('hpo-terms.tsv', 'w') as outfile:
        for key, value in hpo_terms.items():
            outfile.write("{0}\t{1}\t{2}\t{3}\n".format(
                key, value['label'], "|".join(value['lay_person']), value['parents']
            ))

def hpo_to_tree(cls, hpo_terms, hpo_graph, tree, path):
    tree_path = copy.copy(path)
    tree_path.append(cls)
    curie_util = CurieUtil(curie_map.get())
    if cls not in hpo_terms:
        hpo_terms[cls] = {
            'label': hpo_graph.label(URIRef(curie_util.get_uri(cls)))
        }
        parents = hpo_graph.objects(URIRef(curie_util.get_uri(cls)), RDFS.subClassOf)
        hpo_terms[cls]['parents'] = len(list(parents))

        lay_person = get_lay_person(cls, hpo_graph)
        hpo_terms[cls]["lay_person"] = lay_person

    # Traverse the tree to get to the input class
    position = tree[tree_path[0]]
    for term in tree_path[1:]:
        position = position[term]

    for sub_class in hpo_graph.subjects(RDFS.subClassOf, URIRef(curie_util.get_uri(tree_path[-1]))):
        curie = curie_util.get_curie(sub_class).replace("OBO:HP_", "HP:")
        position[curie] = {}
        hpo_to_tree(curie, hpo_terms, hpo_graph, tree, tree_path)


def get_lay_person(cls, hpo_graph):
    lay_synonyms = []
    lay_person_query = """
        SELECT ?lp
        WHERE {{
            ?bnode owl:annotatedSource {0} .
            ?bnode oboInOwl:hasSynonymType <http://purl.obolibrary.org/obo/hp.owl#layperson> .
            ?bnode owl:annotatedTarget ?lp .
        }}
    """.format(cls)
    query_result = hpo_graph.query(lay_person_query)
    if len(list(query_result)) > 0:
        for res in query_result:
            lay_synonyms.append(str(res[0]))
    return lay_synonyms

if __name__ == "__main__":
    main()
