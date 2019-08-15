from abc import ABCMeta, abstractmethod
import re


class Graph(metaclass=ABCMeta):

    # regular expression pattern to recognize strings which could be well formed curies
    # I am still not subscribing to the whole unicode realm
    # while 100% of the identifiers I actually get are 7-bit ascii

    # https://www.w3.org/TR/curie/
    # https://www.w3.org/TR/xml/
    # https://www.w3.org/TR/rdfa-core/#s_curies
    # TEC:
    #   adding a hyphen as a non-first-last char in the local_id portion of curie
    #   because that is reality on the street (i.e Reactome)

    curie_regexp = re.compile(
        r'^[a-zA-Z_]*[a-zA-Z_0-9-]*:[A-Za-z0-9_][A-Za-z0-9_.-]*[A-Za-z0-9_]*$')

    @abstractmethod
    def addTriple(
            self,
            subject_id,
            predicate_id,
            object_id,
            object_is_literal=False,
            literal_type=None):
        pass

    @abstractmethod
    def skolemizeBlankNode(self, curie):
        pass

    @abstractmethod
    def serialize(self, **kwargs):
        pass
