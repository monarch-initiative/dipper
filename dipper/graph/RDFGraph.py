from rdflib import ConjunctiveGraph, Literal, URIRef, BNode, Namespace
from dipper.graph.Graph import Graph as DipperGraph
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map
import re
import logging

logger = logging.getLogger(__name__)


class RDFGraph(ConjunctiveGraph, DipperGraph):
    """
    Extends RDFLibs ConjunctiveGraph
    The goal of this class is wrap the creation
    of triples and manage creation of URIRef,
    Bnodes, and literals from an input curie
    """

    curie_util = CurieUtil(curie_map.get())
    curie_map = curie_map

    def __init__(self, are_bnodes_skized=True):
        super().__init__()
        self.are_bnodes_skized = are_bnodes_skized

        # Can be removed when this is resolved
        # https://github.com/RDFLib/rdflib/issues/632
        obo_map = curie_map.get()['OBO']
        self.bind('OBO', Namespace(obo_map))

    def addTriple(self, subject_id, predicate_id, obj,
                  object_is_literal=False, literal_type=None):

        if object_is_literal is True:
            if literal_type is not None and obj is not None:
                literal_type_iri = self._getNode(literal_type)
                self.add(
                    (self._getNode(subject_id), self._getNode(predicate_id),
                     Literal(obj, datatype=literal_type_iri)))
            elif obj is not None:
                self.add(
                    (self._getNode(subject_id), self._getNode(predicate_id),
                     Literal(obj)))
            else:
                logger.warn(
                    "None as literal object for subj: %s and pred: %s",
                    subject_id, predicate_id)
        elif obj is not None and obj != '':
            self.add(
                (self._getNode(subject_id), self._getNode(predicate_id),
                 self._getNode(obj)))
        else:
            logger.warn(
                "None/empty object IRI for subj: %s and pred: %s",
                subject_id, predicate_id)
        return

    def skolemizeBlankNode(self, curie):
        stripped_id = re.sub(r'^_:|^_', '', curie, 1)
        node = BNode(stripped_id).skolemize(self.curie_map.get_base())
        node = re.sub(r'rdflib/', '', node)
        return URIRef(node)

    def _getNode(self, curie):
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
        if re.match(r'^_', curie):
            if self.are_bnodes_skized is True:
                node = self.skolemizeBlankNode(curie)
            else:  # replace the leading underscore to make it cleaner
                node = BNode(re.sub(r'^_:|^_', '', curie, 1))
        # Check if curie actually an IRI
        elif re.match(r'^http|^ftp', curie):
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
                logger.error("couldn't make URI for %s", curie)
        return node

    def bind_all_namespaces(self):
        for prefix in curie_map.get().keys():
            iri = curie_map.get()[prefix]
            self.bind(prefix, Namespace(iri))
