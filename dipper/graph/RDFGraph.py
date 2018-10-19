import re
import logging
import sys
import yaml

from rdflib import ConjunctiveGraph, Literal, URIRef, BNode, Namespace
from dipper.graph.Graph import Graph as DipperGraph
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map

LOG = logging.getLogger(__name__)

# have yet to notice a feature of conjunctive graph used
# perhaps it was more aspirational
# I think dropping to a single graph per ingest
# and leaving more complicated processes to dedicated graph engines
# would be possible here. TEC


class RDFGraph(ConjunctiveGraph, DipperGraph):
    """
    Extends RDFLibs ConjunctiveGraph
    The goal of this class is wrap the creation
    of triples and manage creation of URIRef,
    Bnodes, and literals from an input curie
    """

    curie_map = curie_map.get()
    curie_util = CurieUtil(curie_map)

    # make global translation table available outside the ingest
    with open('translationtable/GLOBAL_TERMS.yaml') as fhandle:
        globaltt = yaml.safe_load(fhandle)
        globaltcid = {v: k for k, v in globaltt.items()}

    def __init__(self, are_bnodes_skized=True, identifier=None):
        # print("in RDFGraph  with id: ", identifier)
        super().__init__('IOMemory', identifier)
        self.are_bnodes_skized = are_bnodes_skized

        # Can be removed when this is resolved
        # https://github.com/RDFLib/rdflib/issues/632
        obo_map = curie_map.get()['OBO']
        self.bind('OBO', Namespace(obo_map))

        # try adding them all
        # self.bind_all_namespaces()  # too much

    def addTriple(
            self, subject_id, predicate_id, obj, object_is_literal=False,
            literal_type=None):

        if object_is_literal is True:
            if literal_type is not None and obj is not None:
                literal_type_iri = self._getnode(literal_type)
                self.add(
                    (self._getnode(subject_id), self._getnode(predicate_id),
                     Literal(obj, datatype=literal_type_iri)))
            elif obj is not None:
                self.add((
                    self._getnode(subject_id), self._getnode(predicate_id),
                    Literal(obj)))
            else:
                LOG.warning(
                    "None as literal object for subj: %s and pred: %s",
                    subject_id, predicate_id)

                # magic number here is "steps up the stack"
                LOG.warning('\tfrom: %s', sys._getframe(2).f_code.co_name)
                LOG.warning('\t\tfrom: %s', sys._getframe(1).f_code.co_name)
                LOG.warning('\t\t\tfrom: %s', sys._getframe(0).f_code.co_name)

        elif obj is not None and obj != '':  # object is a resourse
            self.add((
                self._getnode(subject_id),
                self._getnode(predicate_id),
                self._getnode(obj)))
        else:
            LOG.warning(
                "None/empty object IRI for subj: %s and pred: %s",
                subject_id, predicate_id)
        return

    def skolemizeBlankNode(self, curie):
        stripped_id = re.sub(r'^_:|^_', '', curie, 1)
        node = BNode(stripped_id).skolemize(self.curie_util.get_base())
        node = re.sub(r'rdflib/', '', node)  # remove string added by rdflib
        return URIRef(node)

    def _getnode(self, curie):  # convention is lowercase names
        """
        This is a wrapper for creating a URIRef or Bnode object
        with a given a curie or iri as a string.

        If an id starts with an underscore, it assigns it to a BNode, otherwise
        it creates it with a standard URIRef.
        Alternatively, self.skolemize_blank_node is True,
        it will skolemize the blank node

        :param curie: str identifier formatted as curie or iri
        :return: node: RDFLib URIRef or BNode object
        """
        node = None
        if curie[0] == '_':
            if self.are_bnodes_skized is True:
                node = self.skolemizeBlankNode(curie)
            else:  # delete the leading underscore to make it cleaner
                node = BNode(re.sub(r'^_:|^_', '', curie, 1))

        # Check if curie string is actually an IRI
        elif curie[:4] == 'http' or curie[:3] == 'ftp':
            node = URIRef(curie)
        else:
            iri = RDFGraph.curie_util.get_uri(curie)
            if iri is not None:
                node = URIRef(RDFGraph.curie_util.get_uri(curie))
                # Bind prefix map to graph
                prefix = curie.split(':')[0]
                if prefix not in self.namespace_manager.namespaces():
                    mapped_iri = curie_map.get()[prefix]
                    self.bind(prefix, Namespace(mapped_iri))
            else:
                LOG.error("couldn't make URI for %s", curie)
        return node

    def bind_all_namespaces(self):
        '''
            Results in the RDF @prefix directives for every ingest
            being added to this ingest.

        '''
        for prefix in curie_map.get().keys():
            iri = curie_map.get()[prefix]
            self.bind(prefix, Namespace(iri))
        return

    # serialize() conflicts between rdflib & Graph.serialize abstractmethod
    # GraphUtils expects the former.  (too bad there is no multiple dispatch)
    #
    # def serialize(
    #        self, subject_iri, predicate_iri, obj, object_is_literal, literal_type):
    #    '''
    #        abstract in parent class. yet to be implemented here
    #
    #    '''
    #    ConjunctiveGraph.serialize(
    #        subject_iri, predicate_iri, obj, object_is_literal, literal_type)
    #    # raise NotImplementedError
