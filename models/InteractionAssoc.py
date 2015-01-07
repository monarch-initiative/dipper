__author__ = 'nicole'

from rdflib import URIRef, Namespace
from rdflib.namespace import RDF, DC
from utils.CurieUtil import CurieUtil
import re
from models.Assoc import Assoc

class InteractionAssoc(Assoc):

    relationships = {
        'genetically_interacts_with' : 'RO:0002435',
        'interacts_with' : 'RO:0002434',  #use this for directly interacts with.  better choice? psi-mi:"MI:0407"(direct interaction)
        'molecularly_interacts_with' : 'RO:0002436',  #should we use this instead for direct interaction?
        'colocalizes_with' : 'RO:0002325', #psi-mi:"MI:0403"(colocalization)
        'ubiquitinates' : 'RO:0002480'
    }

    def __init__(self,assoc_id, gene1, gene2, pub, evidence_code, curie_map):
        if (self.curie_map is None):
            self.curie_map = curie_map
        else:
            self.curie_map.update(curie_map)
        self.cu = CurieUtil(curie_map)
        self.annot_id = assoc_id
        self.gene1 = gene1
        self.gene2 = gene2
        self.pub_id = pub
        self.evidence = evidence_code
        self.rel = self.relationships['interacts_with']  # default
        self.curie_map = curie_map
        self.cu = CurieUtil(self.curie_map)

        self.setSubject(gene1)
        self.setObject(gene2)

        return

    def addInteractionAssociationToGraph(self,g):

        self.addAssociationToGraph(g)

        #todo add some other stuff related to the experimental methods?

        return