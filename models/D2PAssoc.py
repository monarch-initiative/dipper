__author__ = 'nicole'

from rdflib.namespace import OWL, RDF, DC
from rdflib import Namespace, URIRef, Literal
import re
import urllib
from models.Assoc import Assoc
from utils.CurieUtil import CurieUtil

# This one is specific for making a disease-to-phenotype
class D2PAssoc(Assoc):
    '''
    A specific association class for defining Disease-to-Phenotype relationships
    This assumes that a graph is created outside of this class, and nodes get added.
    By default, an association will assume the "has_phenotype" relationship, unless
    otherwise specified.


    '''

    def __init__(self, assoc_id, entity_id, phenotype_id, onset, frequency, pub, evidence_code, curie_map):
        self.annot_id = assoc_id
        self.entity_id = entity_id
        self.phenotype_id = phenotype_id
        self.onset = onset
        self.frequency = frequency
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.relationships['has_phenotype']
        self.curie_map = curie_map
        self.cu = CurieUtil(self.curie_map)

        return

    def set_relationship(self, rel):
        self.rel = rel
        return

    def addAssociationNodeToGraph(self, g):
        '''
        The reified relationship between a disease and a phenotype is decorated with some provenance information.
        This makes the assumption that both the disease and phenotype are classes.

        currently hardcoded to map the annotation to the monarch namespace
        :param g:
        :return:
        '''
        namespaces = self.curie_map
        #FIXME generalize for outside the monarch namespace?
        n = Namespace(namespaces['MONARCH'])
        node = n[self.annot_id]
        s = URIRef(self.cu.get_uri(self.entity_id))
        p = URIRef(self.rel)
        o = URIRef(self.cu.get_uri(self.phenotype_id))

        if (re.compile('http').match(self.pub_id)):
            source = URIRef(self.pub_id)
        else: (re.compile('PMID').match(self.pub_id))

        evidence = URIRef(self.cu.get_uri(self.evidence))
        frequency = onset = None

        g.add((s, RDF['type'], self.OWLCLASS))
        g.add((o, RDF['type'], self.OWLCLASS))
        g.add((s, p, o))

        g.add((node, RDF['type'], URIRef(self.cu.get_uri('Annotation:'))))
        g.add((node, self.BASE['hasSubject'], s))
        g.add((node, self.BASE['hasObject'], o))

        #this is handling the occasional messy pubs that are sometimes literals
        if (self.pub_id.strip() != ''):
            if (source != URIRef('[]')):
                g.add((node, DC['source'], source))
                g.add((source, RDF['type'], self.OWLIND))
            else:
                print("WARN: source as a literal -- is this ok?")
                g.add((node, DC['source'], Literal(self.pub_id)))
            #            else:
            #                print("WARN:",self.entity_id,'+',self.phenotype_id,'has no source information for the association (',self.evidence,')')

        if (self.evidence is None or self.evidence.strip() == ''):
            print("WARN:", self.entity_id, '+', self.phenotype_id, 'has no evidence code')
        else:
            g.add((node, DC['evidence'], evidence))

        if (self.frequency is not None and self.frequency != ''):
            #FIXME what is the real predicate here?
            g.add((node, self.BASE['frequencyOfPhenotype'], Literal(self.frequency)))
        if (self.onset is not None and self.onset != ''):
            #FIXME what is the real predicate here?
            g.add((node, self.BASE['onset'], n[self.onset]))

        return g

