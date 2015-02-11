__author__ = 'nicole'

from rdflib import Literal, URIRef, BNode
from rdflib.namespace import DC, RDF, RDFS, OWL, XSD
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil
import re

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

        n = self._getNode(id)

        g.add((n, RDF['type'], Assoc.OWLCLASS))
        if (label is not None):
            g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(self.cu.get_uri(type))
            g.add((n, Assoc.SUBCLASS, t))
        if (description is not None):
            g.add((n, DC['description'], Literal(description)))
        return g

    def addIndividualToGraph(self, g, id, label, type=None, description=None):
        n = self._getNode(id)

        if (label is not None):
            g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(self.cu.get_uri(type))
            g.add((n, RDF['type'], t))
        else:
            g.add((n, RDF['type'], Assoc.OWLIND))
        if (description is not None):
            g.add((n, DC['description'], Literal(description)))
        return g

    def addEquivalentClass(self,g,id1,id2):
        n1 = self._getNode(id1)
        n2 = self._getNode(id2)

        if (n1 is not None and n2 is not None):
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

    def addSubclass(self,g,parentid,childid):
        p = URIRef(self.cu.get_uri(parentid))
        c = URIRef(self.cu.get_uri(childid))
        g.add((c,Assoc.SUBCLASS,p))

        return

    def addSynonym(self,g,cid,synonym,synonym_type=None):
        '''
        Add the synonym as a property of the class cid.  Assume it is an exact synonym, unless
        otherwise specified
        :param g:
        :param cid: class id
        :param synonym: the literal synonym label
        :param synonym_type: the CURIE of the synonym type (not the URI)
        :return:
        '''
        n = self._getNode(cid)
        if (synonym_type is None):
            synonym_type = URIRef(self.cu.get_uri(Assoc.relationships['hasExactSynonym'])) #default
        else:
            synonym_type = URIRef(self.cu.get_uri(synonym_type))

        g.add((n, synonym_type, Literal(synonym)))
        return

    def addDefinition(self,g,cid,definition):
        if (definition is not None):
            n = self._getNode(cid)
            p = URIRef(self.cu.get_uri(Assoc.relationships['definition']))
            g.add((n,p,Literal(definition)))

        return

    def addXref(self,g,cid,xrefid):
        n1 = self._getNode(cid)
        n2 = self._getNode(xrefid)
        p = URIRef(self.cu.get_uri(Assoc.relationships['has_xref']))
        if (n1 is not None and n2 is not None):
            g.add((n1,p,n2))

        return

    def write(self, graph, format=None, file=None):
        """
         a basic graph writer (to stdout) for any of the sources.  this will write
         raw triples in rdfxml, unless specified.
         to write turtle, specify format='turtle'
         an optional file can be supplied instead of stdout
        :return: None
        """
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


    def _getNode(self,id):
        """
        This is a wrapper for creating a node with a given identifier.  If an id starts with an
        underscore, it assigns it to a BNode, otherwise it creates it with a standard URIRef.
        This will return None if it can't map the node properly.
        :param id:
        :return:
        """
        n=None
        if (re.match('^_',id)):
            n = BNode(id)
        else:
            u = self.cu.get_uri(id)
            if (u is not None):
                n = URIRef(self.cu.get_uri(id))
            else:
                print("ERROR: couldn't make URI for ",id)
        return n