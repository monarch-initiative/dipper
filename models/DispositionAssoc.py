__author__ = 'nicole'

from rdflib import BNode, ConjunctiveGraph, Literal, RDF, OWL, extras, Namespace, URIRef
from rdflib.namespace import FOAF, RDFS, DC

from models.Assoc import Assoc
import re
import urllib

#The first of many associations
#This one is specific for making a disease-to-phenotype
class DispositionAssoc(Assoc):
        #FIXME init might not require all these elements
        has_disposition="GENO:0000208"
        CLASS="owl:Class"
        IND="owl:NamedIndividual"
        monarch_namespace = Namespace('http://www.monarchinitiative.org/')

        prefixes={'OMIM','http://www.omim.org/',
                  'DECIPHER','http://www.decipher.org/',
                  'ORPHANET','http://www.orpha.net/',
                  'HP','http://purl.obolibrary.org/obo/'}

        def __init__(self,annot_id,entity_id,heritability_id,pub,evidence_code):
            self.annot_id = annot_id
            self.entity_id = entity_id
            self.heritability_id = heritability_id
            self.pub_id = pub
            self.evidence = evidence_code
            self.rel = self.has_disposition  #default to has_disposition
            return

        def set_relationship(self,rel):
            self.rel = rel
            return

        def printTypes(self,filewriter=None):
            numStatements = 4
            print(self.entity_id,"rdf:type",self.CLASS,file=filewriter)
            print(self.heritability_id,"rdf:type",self.CLASS,file=filewriter)
            print(self.pub_id,"rdf:type",self.IND,file=filewriter)
            print(self.evidence,"rdf:type",self.CLASS,file=filewriter)
            return numStatements

        def printAssociation(self,filewriter=None):
            numStatements = 6
            #note we are making the assumption that all classes will have their
            #labels printed elsewhere
            #this will return the number of statements printed
            print(self.entity_id,self.rel,self.heritability_id,".", sep=" ",file=filewriter)
            #TODO is this "type" correct?
            #TODO do we need the relationship as a predicate here?
            print(self.annot_id, "rdf:type","Annotation:",".", sep=" ",file=filewriter)
            print(self.annot_id, ":hasSubject", self.entity_id, sep=" ",file=filewriter)
            print(self.annot_id, ":hasObject", self.heritability_id, sep=" ",file=filewriter)
            print(self.annot_id, "dc:source", self.pub_id, sep=" ",file=filewriter)
            print(self.annot_id, "dc:evidence", self.evidence, sep=" ", file=filewriter)
            return numStatements

        def addAssociationNodeToGraph(self,g):
            namespaces = self.namespaces
            n = Namespace(namespaces['MONARCH'])
            node = n[self.annot_id]
            s = n[self.entity_id.replace(':','_')]
            p = n[self.rel]
            o = Namespace(namespaces['HP'])[self.heritability_id.replace('HP:','')]

            if (re.compile('http').match(self.pub_id)):
                source = URIRef(urllib.parse.quote_plus(self.pub_id))
            elif (re.compile('PMID').match(self.pub_id)):
                source = Namespace(namespaces['PMID'])[self.pub_id.replace('PMID:','')]
            else:
                source = Literal(self.pub_id)
            evidence = Namespace(namespaces['ECO'])[self.evidence.replace('ECO:','')]

            g.add((s,RDF['type'],OWL['class']))
            g.add((o,RDF['type'],OWL['class']))
            g.add((node, RDF['type'],OWL['annotation']))
            g.add((node, OWL['hasSubject'], s))
            g.add((node, OWL['hasObject'], o))
            g.add((node, DC['source'], source))
            g.add((node, DC['evidence'], evidence))
            g.add((s,p,o))

            return g