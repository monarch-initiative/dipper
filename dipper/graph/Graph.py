from abc import ABCMeta, abstractmethod


class Graph(metaclass=ABCMeta):

    @abstractmethod
    def addTriple(self):
        pass

    @abstractmethod
    def skolemizeBlankNode(self):
        pass

    @abstractmethod
    def serialize(self):
        pass
