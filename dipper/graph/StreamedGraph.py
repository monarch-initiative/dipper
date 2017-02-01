from dipper.graph.Graph import Graph as DipperGraph
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map
import logging
import re

logger = logging.getLogger(__name__)


class StreamedGraph(DipperGraph):
    """
    Stream rdf triples to file or stdout
    Assumes a downstream process will sort then uniquify triples

    Theoretically could support both ntriple, rdfxml formats, for now
    just support nt
    """

    curie_util = CurieUtil(curie_map.get())
    curie_map = curie_map

    def __init__(self, are_bnodes_skized=False, file_handle=None, fmt='nt'):
        self.are_bnodes_skized = are_bnodes_skized
        self.fmt = fmt
        self.file_handle = file_handle

    def addTriple(self, subject_id, predicate_id, object_id,
                  object_is_literal=False, literal_type=None):
        subject_iri = self._getNode(subject_id)
        predicate_iri = self._getNode(predicate_id)
        if not object_is_literal:
            obj = self._getNode(object_id)
        else:
            obj = object_id

        if literal_type is not None:
            lit_type = self._getNode(literal_type)

        self.serialize(subject_iri, predicate_iri, obj,
                       object_is_literal, lit_type)
        return

    def skolemizeBlankNode(self, curie):
        base_iri = StreamedGraph.curie_map.get_base()
        curie_id = curie.split(':')[1]
        skolem_iri = "{0}.wellknown/genid/{1}".format(base_iri, curie_id)
        return skolem_iri

    def serialize(self, subject_iri, predicate_iri, obj,
                  object_is_literal=False, literal_type=None):
        pass

    def _getNode(self, curie):
        """
        Returns IRI, or blank node curie/iri depending on
        self.skolemize_blank_node setting

        :param curie: str id as curie or iri
        :return:
        """
        if re.match(r'^_:', curie):
            if self.are_bnodes_skized is True:
                node = self.skolemizeBlankNode(curie)
            else:
                node = curie
        elif re.match(r'^http|^ftp', curie):
            node = curie
        elif len(curie.split(':')) == 2:
            node = StreamedGraph.curie_util.get_uri(curie)
        else:
            logger.error("Cannot process curie {}".format(curie))
        return node