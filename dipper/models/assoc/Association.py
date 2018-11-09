
import logging
from dipper.models.Model import Model
from dipper.graph.Graph import Graph
from dipper.utils.GraphUtils import GraphUtils

__author__ = 'nlw'

LOG = logging.getLogger(__name__)
# note: currently no log issued


class Assoc:
    """
    A base class for OBAN (Monarch)-style associations,
    to enable attribution of source and evidence
    on statements.

    """

    def __init__(self, graph, definedby, sub=None, obj=None, pred=None):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".format(graph))
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        # core parts of the association
        self.definedby = definedby
        self.sub = sub
        self.obj = obj
        self.rel = pred
        self.assoc_id = None

        self.description = None
        self.source = []
        self.evidence = []
        self.date = []

        # this is going to be used for the refactored evidence/provenance
        self.provenance = []
        self.score = None
        self.score_type = None
        self.score_unit = None

        return

    def _is_valid(self):
        # check if sub/obj/rel are none...raise error
        if self.sub is None:
            raise ValueError(
                'No subject set for this association <%s> <%s> <%s>',
                self.sub, self.rel, self.obj
            )
        if self.obj is None:
            raise ValueError(
                'No object set for this association <%s> <%s> <%s>',
                self.sub, self.rel, self.obj
            )
        if self.rel is None:
            raise ValueError(
                'No predicate set for this association <%s> <%s> <%s>',
                self.sub, self.rel, self.obj
            )
        # Are subject & predicate, either a curie or IRI
        pfx = self.sub.split(':')[0]
        if pfx not in self.curie_map.keys() and pfx not in ['_', 'http', 'https', 'ftp']:
            raise ValueError(
                'Invalid Subject for this association <%s> <%s> <%s>',
                self.sub, self.rel, self.obj
            )
        pfx = self.rel.split(':')[0]
        if pfx not in self.curie_map.keys() and pfx not in ['_', 'http', 'https', 'ftp']:
            raise ValueError(
                'Invalid Predicate for this association <%s> <%s> <%s>',
                self.sub, self.rel, self.obj
            )

        return True

    def add_association_to_graph(self):

        if not self._is_valid():
            return

        self.graph.addTriple(self.sub, self.rel, self.obj)

        if self.assoc_id is None:
            self.set_association_id()

        assert self.assoc_id is not None

        self.model.addType(self.assoc_id, self.model.globaltt['association'])

        self.graph.addTriple(
            self.assoc_id, self.globaltt['association has subject'], self.sub)
        self.graph.addTriple(
            self.assoc_id, self.globaltt['association has object'], self.obj)
        self.graph.addTriple(
            self.assoc_id, self.globaltt['association has predicate'], self.rel)

        if self.description is not None:
            self.model.addDescription(self.assoc_id, self.description)

        if self.evidence is not None and len(self.evidence) > 0:
            for evi in self.evidence:
                self.graph.addTriple(
                    self.assoc_id, self.globaltt['has evidence'], evi)

        if self.source is not None and len(self.source) > 0:
            for src in self.source:
                object_is_literal = False
                if src[:4] == 'http':   # not a curie
                    object_is_literal = True
                    # TODO assume that the source is a publication?
                    # use Reference class
                self.graph.addTriple(
                    self.assoc_id, self.globaltt['source'], src, object_is_literal)

        if self.provenance is not None and len(self.provenance) > 0:
            for prov in self.provenance:
                self.graph.addTriple(
                    self.assoc_id, self.globaltt['has_provenance'], prov)

        if self.date is not None and len(self.date) > 0:
            for dat in self.date:
                self.graph.addTriple(
                    object_is_literal=True, subject_id=self.assoc_id,
                    predicate_id=self.globaltt['created_on'], obj=dat)

        if self.score is not None:
            self.graph.addTriple(
                self.assoc_id, self.globaltt['has measurement value'], self.score,
                True, 'xsd:float')
            # TODO
            # update with some kind of instance of scoring object
            # that has a unit and type

        return

    def add_predicate_object(
            self, predicate, object_node, object_type=None, datatype=None):

        if object_type == 'Literal':
            if datatype is not None:
                self.graph.addTriple(
                    self.assoc_id, predicate, object_node, True, datatype)
            else:
                self.graph.addTriple(self.assoc_id, predicate, object_node, True)
        else:
            self.graph.addTriple(self.assoc_id, predicate, object_node, False)

        return

    # This isn't java, but predecessors favored the use of property decorators
    # and CamelCase and ...
    def set_subject(self, identifier):
        self.sub = identifier
        return

    def set_object(self, identifier):
        self.obj = identifier
        return

    def set_relationship(self, identifier):
        self.rel = identifier
        return

    def set_association_id(self, assoc_id=None):
        """
        This will set the association ID based on the internal parts
        of the association.
        To be used in cases where an external association identifier
        should be used.

        :param assoc_id:

        :return:

        """
        if assoc_id is None:
            self.assoc_id = self.make_association_id(
                self.definedby, self.sub, self.rel, self.obj)
        else:
            self.assoc_id = assoc_id

        return self.assoc_id

    def get_association_id(self):
        if self.assoc_id is None:
            self.set_association_id()

        return self.assoc_id

    def set_description(self, description):
        self.description = description

        return

    def set_score(self, score, unit=None, score_type=None):

        self.score = score
        self.score_unit = unit
        self.score_type = score_type

        return

    def add_evidence(self, identifier):
        """
        Add an evidence code to the association object (maintained as a list)
        :param identifier:

        :return:

        """

        if identifier is not None and identifier.strip() != '':
            self.evidence += [identifier]

        return

    def add_source(self, identifier):
        """
        Add a source identifier (such as publication id)
        to the association object (maintained as a list)
        TODO we need to greatly expand this function!

        :param identifier:

        :return:

        """

        if identifier is not None and identifier.strip() != '':
            self.source += [identifier]

        return

    def add_date(self, date):
        if date is not None and date.strip() != '':
            self.date += [date]

        return

    def add_provenance(self, identifier):

        if identifier is not None and identifier.strip() != '':
            self.provenance += [identifier]

        return

    @staticmethod
    def make_association_id(definedby, sub, pred, obj, attributes=None):
        """
        A method to create unique identifiers for OBAN-style associations,
        based on all the parts of the association
        If any of the items is empty or None, it will convert it to blank.
        It effectively digests the  string of concatonated values.
        Subclasses of Assoc can submit an additional array of attributes
        that will be appeded to the ID.

        Note this is equivalent to a RDF blank node

        :param definedby: The (data) resource that provided the annotation
        :param subject:
        :param predicate:
        :param object:
        :param attributes:

        :return:

        """

        items_to_hash = [definedby, sub, pred, obj]
        if attributes is not None and len(attributes) > 0:
            items_to_hash += attributes

        items_to_hash = [x for x in items_to_hash if x is not None]

        assoc_id = ':'.join(('MONARCH', GraphUtils.digest_id('+'.join(items_to_hash))))
        assert assoc_id is not None
        return assoc_id
