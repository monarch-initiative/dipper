from abc import ABCMeta, abstractmethod


class Graph(metaclass=ABCMeta):

    @abstractmethod
    def addTriple(self, subject_id, predicate_id, object_id,
                  object_is_literal, literal_type):
        pass

    @abstractmethod
    def skolemizeBlankNode(self, curie):
        pass

    @abstractmethod
    def serialize(self, subject_iri, predicate_iri, obj,
                  object_is_literal, literal_type):
        pass
