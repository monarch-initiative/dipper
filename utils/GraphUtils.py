__author__ = 'nicole'

from rdflib import Literal, URIRef
from rdflib.namespace import DC, RDF, RDFS, OWL
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil

class GraphUtils:

    def __init__(self,curie_map):
        self.curie_map = curie_map
        self.cu = CurieUtil(curie_map)
        return

    def addClassToGraph(self, g, id, label, type=None, description=None):

        '''
        Any node added to the graph will get at least 3 triples:
        *(node,type,owl:Class) and
        *(node,label,literal(label))
        *if a type is added, then the node will be an OWL:subclassOf that the type
        *if a description is provided, it will also get added as a dc:description
        :param id:
        :param label:
        :param type:
        :param description:
        :return:
        '''
        n = URIRef(self.cu.get_uri(id))

        g.add((n, RDF['type'], Assoc.OWLCLASS))
        g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(self.cu.get_uri(type))
            g.add((n, Assoc.OWLSUBCLASS, t))
        if (description is not None):
            g.add((n, DC['description'], Literal(description)))
        return g

    def addEquivalentClass(self,g,id1,id2):
        n1 = URIRef(self.cu.get_uri(id1))
        n2 = URIRef(self.cu.get_uri(id2))
        g.add((n1,OWL['equivalentClass'],n2))
        return