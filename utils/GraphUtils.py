__author__ = 'nicole'

from rdflib import Literal, URIRef
from rdflib.namespace import DC, RDF, RDFS, OWL, XSD
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
        if (label is not None):
            g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(self.cu.get_uri(type))
            g.add((n, Assoc.OWLSUBCLASS, t))
        if (description is not None):
            g.add((n, DC['description'], Literal(description)))
        return g

    def addEquivalentClass(self,g,id1,id2):
        if (self.cu.get_uri(id1) is not None and self.cu.get_uri(id2) is not None):
            n1 = URIRef(self.cu.get_uri(id1))
            n2 = URIRef(self.cu.get_uri(id2))
            g.add((n1,OWL['equivalentClass'],n2))
        return

    def addDeprecatedClass(self,g,oldid,newids=None):
        '''
        Will mark the oldid as a deprecated class.
        if one newid is supplied, it will mark it as replaced by.
        if >1 newid is supplied, it will mark it with consider properties
        :param g:
        :param oldid: the class id to deprecate
        :param newids: the class idlist that is the replacement(s) of the old class.  Not required.
        :return:
        '''
        #print("INFO: adding deprecated class:",oldid, "with newids",newids)
        consider = URIRef(self.cu.get_uri(Assoc.relationships['consider']))
        replaced_by = URIRef(self.cu.get_uri(Assoc.relationships['replaced_by']))

        n1 = URIRef(self.cu.get_uri(oldid))
        g.add((n1,RDF['type'],Assoc.OWLCLASS))
        g.add((n1,OWL['deprecated'],Literal(True,datatype=XSD[bool])))

        if (newids is not None):
            if (len(newids) == 1):
                n = URIRef(self.cu.get_uri(newids[0]))
                g.add((n1,replaced_by,n))
            elif (len(newids) > 0):
                for i in newids:
                    n = URIRef(self.cu.get_uri(i.strip()))
                    g.add((n1,consider,n))
        return

    def addSynonym(self,g,cid,synonym,synonym_type=None):
        n = URIRef(self.cu.get_uri(cid))
        if (synonym_type is None):
            synonym_type = URIRef(self.cu.get_uri(Assoc.relationships['hasExactSynonym'])) #default
        g.add((n, synonym_type, Literal(synonym)))

    def addDefinition(self,g,cid,definition):
        n = URIRef(self.cu.get_uri(cid))
        p = URIRef(self.cu.get_uri(Assoc.relationships['definition']))
        g.add((n,p,Literal(definition)))

        return

    def write(self, graph, format=None, file=None):
        '''
         a basic graph writer (to stdout) for any of the sources.  this will write
         raw triples in rdfxml, unless specified.
         to write turtle, specify format='turtle'
         an optional file can be supplied instead of stdout
        :return: None
        '''
        if (format is None):
            format = 'rdfxml'
        if (file is not None):
            filewriter = open(file, 'w')
            print("INFO: Writing triples in",format,"to",file)
            print(graph.serialize(format=format).decode(), file=filewriter)
            filewriter.close()
        else:
            print(graph.serialize(format=format).decode())
        return