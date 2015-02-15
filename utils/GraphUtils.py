__author__ = 'nicole'

from rdflib import Literal, URIRef, BNode, Namespace
from rdflib.namespace import DC, RDF, RDFS, OWL, XSD
from utils.CurieUtil import CurieUtil
import re

class GraphUtils:

    #FIXME - i've duplicated relationships in Assoc and here - pick one or the other and refactor
    #TODO refactor using the getNode() method to clear out the URIRef(cu.get_uri(<id>)) nonsense

    OWLCLASS=OWL['Class']
    OWLIND=OWL['NamedIndividual']
    OWLPROP=OWL['ObjectProperty']
    SUBCLASS=RDFS['subClassOf']

    relationships = {
        'has_disposition':'GENO:0000208',
        'has_phenotype':'RO:0002200',
        'replaced_by' : 'IAO:0100001',
        'consider' : 'OIO:consider',
        'hasExactSynonym' : 'OIO:hasExactSynonym',
        'hasRelatedSynonym' : 'OIO:hasRelatedSynonym',
        'definition' : 'IAO:0000115',
        'in_taxon' : 'RO:0002162',
        'has_quality' : 'RO:0000086',
        'towards' : 'RO:0002503',
        'has_xref' : 'OIO:hasDbXref'
    }


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

        g.add((n, RDF['type'], self.OWLCLASS))
        if (label is not None):
            g.add((n, RDFS['label'], Literal(label)))
        if (type is not None):
            t = URIRef(self.cu.get_uri(type))
            g.add((n, self.SUBCLASS, t))
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
            g.add((n, RDF['type'], self.OWLIND))
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
        consider = URIRef(self.cu.get_uri(self.relationships['consider']))
        replaced_by = URIRef(self.cu.get_uri(self.relationships['replaced_by']))

        n1 = URIRef(self.cu.get_uri(oldid))
        g.add((n1,RDF['type'],self.OWLCLASS))
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
        g.add((c,self.SUBCLASS,p))

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
            synonym_type = URIRef(self.cu.get_uri(self.relationships['hasExactSynonym'])) #default
        else:
            synonym_type = URIRef(self.cu.get_uri(synonym_type))

        g.add((n, synonym_type, Literal(synonym)))
        return

    def addDefinition(self,g,cid,definition):
        if (definition is not None):
            n = self._getNode(cid)
            p = URIRef(self.cu.get_uri(self.relationships['definition']))
            g.add((n,p,Literal(definition)))

        return

    def addXref(self,g,cid,xrefid):
        n1 = self._getNode(cid)
        n2 = self._getNode(xrefid)
        p = URIRef(self.cu.get_uri(self.relationships['has_xref']))
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

    def getNode(self,id):

        return self._getNode(id)


    def loadObjectProperties(self,graph,op):
        """
        Given a graph, it will load the supplied object properties
        as owl['ObjectProperty'] types
        :param g: a graph
        :param op: a dictionary of object properties
        :return: None
        """
        cu = self.cu
        for k in op:
            graph.add((URIRef(cu.get_uri(op[k])),RDF['type'],self.OWLPROP))
        return