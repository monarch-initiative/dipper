__author__ = 'nicole'

import re
import urllib

from rdflib import RDF, Namespace, URIRef
from rdflib.namespace import DC

from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil
from conf import curie_map


#The first of many associations
#This one is specific for making a disease-to-phenotype
class DispositionAssoc(Assoc):
        #FIXME init might not require all these elements
        has_disposition="http://purl.obolibrary.org/obo/GENO_0000208"

        relationships={'has_disposition','http://purl.obolibrary.org/obo/GENO_0000208'}

        def __init__(self,annot_id,entity_id,heritability_id,pub,evidence_code):
            self.annot_id = annot_id
            self.entity_id = entity_id
            self.heritability_id = heritability_id
            self.pub_id = pub
            self.evidence = evidence_code
            self.rel = self.has_disposition  #default to has_disposition
            self.cu = CurieUtil(curie_map.get())
            self.pub_list = None

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
            namespaces = curie_map.get()
            n = Namespace(namespaces['MONARCH'])
            node = n[self.annot_id]

            s = URIRef(self.cu.get_uri(self.entity_id))
            p = URIRef(self.rel)
            o = Namespace(namespaces['HP'])[self.heritability_id.replace('HP:','')]

            if (re.compile('http').match(self.pub_id)):
                source = URIRef(urllib.parse.quote_plus(self.pub_id))
            elif (re.compile('PMID').match(self.pub_id)):
#                source = Namespace(namespaces['PMID'])[self.pub_id.replace('PMID:','')]
                source = URIRef(self.cu.get_uri(self.pub_id))
            else:
                u = self.cu.get_uri(self.pub_id)
                if (u is not None):
                    source = URIRef(u)
                else:
                    source = None

            #evidence = Namespace(namespaces['ECO'])[self.evidence.replace('ECO:','')]
            evidence = URIRef(self.cu.get_uri(self.evidence))

            g.add((s,RDF['type'],self.OWLCLASS))
            g.add((o,RDF['type'],self.OWLCLASS))
            g.add((s,p,o))

            g.add((node, RDF['type'],URIRef(self.cu.get_uri('Annotation:'))))
            g.add((node, self.BASE['hasSubject'], s))
            g.add((node, self.BASE['hasObject'], o))
            if (self.pub_id is None or self.pub_id.strip() == ''):
                print("WARN:",self.entity_id,'+',self.heritability_id,'has no source information for the association (',self.evidence,')')
            else:
                g.add((node, DC['source'], source))
                g.add((source, RDF['type'], self.OWLIND))

            if (self.evidence.strip() == ''):
                print("WARN:",self.entity_id,'+',self.heritability_id,'has no evidence code')
            else:
                g.add((node, DC['evidence'], evidence))


            return g