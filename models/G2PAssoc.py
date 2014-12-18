__author__ = 'nicole'

from rdflib.namespace import OWL, RDF, DC
from rdflib import Namespace, URIRef, Literal
import re
import urllib
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil

#The first of many associations
#This one is specific for making a disease-to-phenotype
class G2PAssoc(Assoc):

        def __init__(self,assoc_id,entity_id,phenotype_id,pub,evidence_code,curie_map):
            self.annot_id=assoc_id
            self.entity_id = entity_id
            self.phenotype_id = phenotype_id
            self.pub_id = pub
            self.evidence = evidence_code
            self.rel = self.relationships['has_phenotype']  #default to has_phenotype
            self.curie_map = curie_map
            self.cu = CurieUtil(self.curie_map)
            return

        def set_relationship(self,rel):
            self.rel = rel
            return

        def set_stage(self,start_stage_id,end_stage_id):
            self.start_stage_id = start_stage_id
            self.end_stage_id = end_stage_id
            return

        def addAssociationNodeToGraph(self,g):
            namespaces = self.curie_map
            n = Namespace(namespaces['MONARCH'])
            node = n[self.annot_id]
            s = URIRef(self.cu.get_uri(self.entity_id))
            p = URIRef(self.rel)
            #TODO generalize mapping curie to namespace
            o = URIRef(self.cu.get_uri(self.phenotype_id))

            if (re.compile('http').match(self.pub_id)):
                source = URIRef(self.pub_id)
            else:
                source = URIRef(self.cu.get_uri(self.pub_id))
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

            #TODO add staging information here
            #if (self.start_stage_id is not None):
            #    g.add((node, self.BASE['during'], URIRef(self.cu.get_uri(self.start_stage_id))))

            return g

