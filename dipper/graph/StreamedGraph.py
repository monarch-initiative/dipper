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

    def __init__(self, are_bnodes_skized=True, file_handle=None, fmt='nt'):
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
            literal_type = self._getNode(literal_type)

        if object_id is not None:
            self.serialize(subject_iri, predicate_iri, obj,
                           object_is_literal, literal_type)
        else:
            logger.warn("Null value passed as object")
        return

    def skolemizeBlankNode(self, curie):
        base_iri = StreamedGraph.curie_map.get_base()
        curie_id = curie.split(':')[1]
        skolem_iri = "{0}.wellknown/genid/{1}".format(base_iri, curie_id)
        return skolem_iri

    def serialize(self, subject_iri, predicate_iri, obj,
                  object_is_literal=False, literal_type=None):
        if not object_is_literal:
            triple = "<{}> <{}> <{}> .".format(subject_iri, predicate_iri, obj)
        elif literal_type is not None:
            triple = '<{}> <{}> {}^^<{}> .'.format(
                subject_iri, predicate_iri,
                self._quote_encode(str(obj)), literal_type)
        else:
            if isinstance(obj, str):
                triple = '<{}> <{}> {} .'.format(
                    subject_iri, predicate_iri, self._quote_encode(obj))
            else:
                lit_type = self._getLiteralXSDType(obj)
                if type is not None:
                    triple = '<{}> <{}> "{}"^^<{}> .'.format(
                        subject_iri, predicate_iri, obj, lit_type)
                else:
                    raise TypeError("Cannot determine type of {}".format(obj))

        if self.file_handle is None:
            print(triple)
        else:
            self.file_handle.write("{}\n".format(triple))

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
            raise TypeError("Cannot process curie {}".format(curie))
        return node

    def _getLiteralXSDType(self, literal):
        """
        This could be much more nuanced, but for now
        if a literal is not a str, determine if it's
        a xsd int or double
        :param literal:
        :return: str - xsd full iri
        """
        if isinstance(literal, int):
            return self._getNode("xsd:integer")
        if isinstance(literal, float):
            return self._getNode("xsd:double")

    @staticmethod
    def _quote_encode(literal):
        """
        Copy of code in rdflib here:
        https://github.com/RDFLib/rdflib/blob/776b90be/
        rdflib/plugins/serializers/nt.py#L76
        :param literal:
        :return:
        """
        return '"%s"' % literal.replace('\\', '\\\\')\
            .replace('\n', '\\n')\
            .replace('"', '\\"')\
            .replace('\r', '\\r')
