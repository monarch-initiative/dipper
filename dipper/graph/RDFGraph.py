import re
import logging
import sys
import os

import yaml
from rdflib import ConjunctiveGraph, Literal, URIRef, BNode, Namespace

from dipper.graph.Graph import Graph as DipperGraph
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map as curie_map_class
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class RDFGraph(DipperGraph, ConjunctiveGraph):
    """
    Extends RDFLibs ConjunctiveGraph
    The goal of this class is wrap the creation
    of triples and manage creation of URIRef,
    Bnodes, and literals from an input curie
    """

    curie_map = curie_map_class.get()
    curie_util = CurieUtil(curie_map)

    # make global translation table available outside the ingest
    with open(
        os.path.join(
            os.path.dirname(__file__),
            '../../translationtable/GLOBAL_TERMS.yaml')) as fhandle:
        globaltt = yaml.safe_load(fhandle)
        globaltcid = {v: k for k, v in globaltt.items()}

    def __init__(self, are_bnodes_skized=True, identifier=None):
        # print("in RDFGraph  with id: ", identifier)
        super().__init__('IOMemory', identifier)
        self.are_bnodes_skized = are_bnodes_skized
        self.prefixes = set()

        # Can be removed when this is resolved
        # https://github.com/RDFLib/rdflib/issues/632
        for pfx in ('OBO',):  # , 'ORPHA'):
            self.bind(pfx, Namespace(self.curie_map[pfx]))

    def _make_category_triple(
            self, subject, category, predicate=blv.terms['category']
    ):
        """
        add a triple to capture subject or object category (in CURIE form) that was
        passed to addTriple()
        """
        try:
            self.add((
                self._getnode(subject),
                self._getnode(predicate),
                self._getnode(category)))
        except:
            LOG.warning(
                "Problem adding triple in _makeCategoryTriple for " + \
                "subj: %s pred: %s obj(category): %s",
                subject, predicate, category)
                
    def _is_literal(self, thing):
        """
        make inference on type (literal or CURIE)

        return: logical
        """
        if self.curie_regexp.match(thing) is not None or\
           thing.split(':')[0].lower() in ('http', 'https', 'ftp'):
            object_is_literal = False
        else:
            object_is_literal = True

        return object_is_literal

    def addTriple(
            self,
            subject_id,
            predicate_id,
            obj,
            object_is_literal=None,
            literal_type=None,
            subject_category=None,
            object_category=None
    ):

        if object_is_literal is None:
            object_is_literal = self._is_literal(obj)

        # add triples for subject category info
        if subject_category is not None:
            self._make_category_triple(subject_id, subject_category)

        # add triples for obj category info, if obj is not a literal
        if not object_is_literal:
            if object_category is not None:
                self._make_category_triple(obj, object_category)
        else: # emit warning if object category is given for a literal
            if object_category is not None:
                LOG.warning("I was given a category %s for obj: %s, " +
                            "which seems to be a literal!",
                            object_category, obj)
            
        if object_is_literal is True:
            if isinstance(obj, str):
                re.sub(r'[\t\n\r\f\v]+', ' ', obj)  # reduce any ws to a space
            if literal_type is not None and obj is not None and obj not in ("", " "):
                literal_type_iri = self._getnode(literal_type)

                self.add(
                    (self._getnode(subject_id), self._getnode(predicate_id),
                     Literal(obj, datatype=literal_type_iri)))
            elif obj is not None:
                # could attempt to infer a type here but there is no use case
                self.add((
                    self._getnode(subject_id), self._getnode(predicate_id),
                    Literal(obj)))
            else:
                LOG.warning(
                    "None as literal object for subj: %s and pred: %s",
                    subject_id, predicate_id)
                # get a sense of where the None is comming from
                # magic number here is "steps up the call stack"
                # TODO there may be easier/ideomatic ways to do this now
                for call in range(2, 0, -1):
                    LOG.warning(
                        '\t%sfrom: %s', '\t' * call, sys._getframe(call).f_code.co_name)

        elif obj is not None and obj != '':  # object is a resource
            self.add((
                self._getnode(subject_id),
                self._getnode(predicate_id),
                self._getnode(obj)))
        else:
            LOG.warning(
                "None/empty object IRI for subj: %s and pred: %s",
                subject_id, predicate_id)

    def skolemizeBlankNode(self, curie):
        stripped_id = re.sub(r'^_:|^_', '', curie, 1)
        return URIRef(stripped_id + self.curie_map['BNODE'])

    def _getnode(self, curie):
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
            if self.are_bnodes_skized:
                node = self.skolemizeBlankNode(curie)
            else:  # delete the leading underscore to make it cleaner
                node = BNode(re.sub(r'^_:|^_', '', curie, 1))

        # Check if curie string is actually an IRI
        elif curie[:4] == 'http' or curie[:3] == 'ftp' or curie[:4] == 'jdbc':
            node = URIRef(curie)
        else:
            iri = RDFGraph.curie_util.get_uri(curie)
            if iri is not None:
                node = URIRef(iri)
                # Bind prefix map to graph
                prefix = curie.split(':')[0]
                self.prefixes.add(prefix)
            else:
                LOG.error("couldn't make URI for %s", curie)
                # get a sense of where the CURIE-ish? thing is comming from
                # magic number here is "steps up the call stack"
                for call in range(3, 0, -1):
                    LOG.warning(
                        '\t%sfrom: %s', '\t' * call, sys._getframe(call).f_code.co_name)
        return node

    def bind_all_namespaces(self):
        """
            Results in the RDF @prefix directives for every ingest
            being added to this ingest.
        """
        for prefix in self.curie_map.keys():
            iri = self.curie_map[prefix]
            self.bind(prefix, Namespace(iri))

    # serialize() conflicts between rdflib & Graph.serialize abstractmethod
    # GraphUtils expects the former.  (too bad there is no multiple dispatch)
    # rdflib version
    def serialize(
            self, destination=None, format='turtle', base=None, encoding=None
    ):
        for prefix in self.prefixes:
            mapped_iri = self.curie_map[prefix]
            self.bind(prefix, Namespace(mapped_iri))
        return ConjunctiveGraph.serialize(self, destination, format)
