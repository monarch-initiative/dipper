__author__ = 'nicole'

from rdflib.namespace import OWL, RDF, DC
from rdflib import Namespace, URIRef, Literal
import re
import urllib
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil

#The first of many associations
#This one is specific for making a disease-to-phenotype
class D2PAssoc(Assoc):
        #FIXME init might not require all these elements
        has_phenotype="http://purl.obolibrary.org/obo/RO_0002200" #has_phenotype

        def __init__(self,assoc_id,entity_id,phenotype_id,onset,frequency,pub,evidence_code,curie_map):
            self.annot_id=assoc_id
            self.entity_id = entity_id
            self.phenotype_id = phenotype_id
            self.onset = onset
            self.frequency = frequency
            self.pub_id = pub
            self.evidence = evidence_code
            self.rel = self.has_phenotype  #default to has_phenotype
            self.curie_map = curie_map
            self.cu = CurieUtil(self.curie_map)

            return

        def set_relationship(self,rel):
            self.rel = rel
            return

        def addAssociationNodeToGraph(self,g):
            namespaces = self.curie_map
            n = Namespace(namespaces['MONARCH'])
            node = n[self.annot_id]
#            s = n[self.entity_id.replace(':','_')]
            s = URIRef(self.cu.get_uri(self.entity_id))
            p = URIRef(self.rel)
            #TODO generalize mapping curie to namespace
#            o = Namespace(namespaces['HP'])[self.phenotype_id.replace('HP:','')]
            o = URIRef(self.cu.get_uri(self.phenotype_id))

            #scrub the occasional case errors in the pubs
            self.pub_id = self.pub_id.replace('pmid','PMID')
            if (re.compile('http').match(self.pub_id)):
                source = URIRef(self.pub_id)
#                source = URIRef(urllib.parse.urlencode(self.pub_id))
                #source = URIRef(urllib.parse.quote_plus(self.pub_id))
            elif (re.compile('PMID').match(self.pub_id)):
#                source = Namespace(namespaces['PMID'])[self.pub_id.replace('PMID:','')]
                source = URIRef(self.cu.get_uri(self.pub_id))
            else:
                source = URIRef(self.cu.get_uri(self.pub_id))
#            print("SOURCE: ",source)
            evidence = URIRef(self.cu.get_uri(self.evidence))
            frequency = onset = None

            g.add((s,RDF['type'],self.OWLCLASS))
            g.add((o,RDF['type'],self.OWLCLASS))
            g.add((s,p,o))

            g.add((node, RDF['type'],URIRef(self.cu.get_uri('Annotation:'))))
            g.add((node, self.BASE['hasSubject'], s))
            g.add((node, self.BASE['hasObject'], o))

            if (self.pub_id.strip() != ''):
                if (source != URIRef('[]')):
                    g.add((node, DC['source'], source))
                    g.add((source, RDF['type'], self.OWLIND))
                else:
                    print("WARN: source as a literal -- is this ok?")
                    g.add((node, DC['source'], Literal(self.pub_id)))
#            else:
#                print("WARN:",self.entity_id,'+',self.phenotype_id,'has no source information for the association (',self.evidence,')')

            if (self.evidence.strip() == ''):
                print("WARN:",self.entity_id,'+',self.phenotype_id,'has no evidence code')
            else:
                g.add((node, DC['evidence'], evidence))

            if (self.frequency is not None and self.frequency != ''):
                #FIXME what is the real predicate here?
                g.add((node, self.BASE['frequencyOfPhenotype'],Literal(self.frequency)))
            if (self.onset is not None and self.onset != ''):
                #FIXME what is the real predicate here?
                g.add((node,self.BASE['onset'],n[self.onset]))

#            print(g.serialize(format="turtle").decode())
            return g

