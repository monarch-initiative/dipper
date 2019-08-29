import logging
import re
import yaml
import os

from dipper.graph.Graph import Graph as DipperGraph
from dipper.utils.CurieUtil import CurieUtil
from dipper import curie_map as curimap

LOG = logging.getLogger(__name__)


class StreamedGraph(DipperGraph):
    """
    Stream rdf triples to file or stdout
    Assumes a downstream process will sort then uniquify triples

    Theoretically could support both ntriple, rdfxml formats, for now
    just support nt
    """

    curie_map = curimap.get()
    curie_util = CurieUtil(curie_map)

    with open(
            os.path.join(os.path.dirname(__file__),
                         '../../translationtable/GLOBAL_TERMS.yaml')) as fhandle:
        globaltt = yaml.safe_load(fhandle).copy()
        globaltcid = {v: k for k, v in globaltt.items()}

    def __init__(
            self, are_bnodes_skized=True, identifier=None, file_handle=None, fmt='nt'):
        self.are_bnodes_skized = are_bnodes_skized
        self.fmt = fmt
        self.file_handle = file_handle
        self.identifier = identifier

    def addTriple(
            self, subject_id, predicate_id, obj, object_is_literal=None,
            literal_type=None, subject_category=None, object_category=None):
        # trying making inference on type of object if none is supplied
        if object_is_literal is None:
            if self.curie_regexp.match(obj) or\
                    obj.split(':')[0].lower() in ('http', 'https', 'ftp'):
                object_is_literal = False
        else:
            object_is_literal = True

        subject_iri = self._getnode(subject_id)
        predicate_iri = self._getnode(predicate_id)
        if not object_is_literal:
            obj = self._getnode(obj)

        if literal_type is not None:
            literal_type = self._getnode(literal_type)

        if subject_category is not None:
            subject_category_iri = self._getnode(subject_category)

        if object_category is not None:
            object_category_iri = self._getnode(object_category)

        if obj is not None:
            self.serialize(
                subject_iri, predicate_iri, obj, object_is_literal, literal_type,
                subject_category_iri=subject_category_iri,
                object_category_iri=object_category_iri)
        else:
            LOG.warning("Null value passed as object")
        return

    def skolemizeBlankNode(self, curie):
        base_iri = StreamedGraph.curie_map.get_base()
        curie_id = curie.split(':')[1]
        skolem_iri = "{0}.wellknown/genid/{1}".format(base_iri, curie_id)
        return skolem_iri

    def serialize(self, subject_iri, predicate_iri, obj,
                  object_is_literal=False, literal_type=None,
                  subject_category_iri=None,
                  predicate_category_iri="biolink:category",
                  object_category_iri=None):
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

        biolink_category_triples = []
        if subject_category_iri is not None:
            biolink_category_triples += \
                "<{}> <{}> <{}> .".format(subject_iri, predicate_category_iri,
                                          subject_category_iri)
        if object_category_iri is not None:
            if object_is_literal or literal_type is not None:
                LOG.warning("can't write biolink category triple for literal!")
            else:
                biolink_category_triples += \
                    "<{}> <{}> <{}> .".format(obj, predicate_category_iri,
                                              object_category_iri)
        all_triples = "\n".join(triple, biolink_category_triples)

        if self.file_handle is None:
            print(all_triples)
        else:
            self.file_handle.write("{}\n".format(all_triples))

    def _make_biolink_category_triples(self):
        return biolink_category_triples

    def _getnode(self, curie):
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
            return self._getnode("xsd:integer")
        if isinstance(literal, float):
            return self._getnode("xsd:double")

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
