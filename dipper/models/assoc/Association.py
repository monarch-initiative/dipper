import re
import logging
import hashlib
from dipper.models.Model import Model
from dipper.graph.Graph import Graph

__author__ = 'nlw'

logger = logging.getLogger(__name__)


class Assoc:
    """
    A base class for OBAN (Monarch)-style associations,
    to enable attribution of source and evidence
    on statements.

    """

    assoc_types = {
        'association': 'OBAN:association'
    }

    annotation_properties = {
        'replaced_by': 'IAO:0100001',
        'consider': 'OIO:consider',
        'hasExactSynonym': 'OIO:hasExactSynonym',
        'hasRelatedSynonym': 'OIO:hasRelatedSynonym',
        'definition': 'IAO:0000115',
        'has_xref': 'OIO:hasDbXref',
        'inchi_key': 'CHEBI:InChIKey',
        'probabalistic_quantifier': 'GENO:0000867'
    }

    object_properties = {
        'has disposition': 'RO:0000091',
        'has_phenotype': 'RO:0002200',
        'expressed_in': 'RO:0002206',
        'in_taxon': 'RO:0002162',
        'has_quality': 'RO:0000086',
        'towards': 'RO:0002503',
        'has_subject': 'OBAN:association_has_subject',
        'has_object': 'OBAN:association_has_object',
        'has_predicate': 'OBAN:association_has_predicate',
        'is_about': 'IAO:0000136',
        'has_evidence': 'RO:0002558',
        'has_source': 'dc:source',
        'has_provenance': 'OBAN:has_provenance',
        'causes_or_contributes': 'RO:0003302'
    }

    datatype_properties = {
        'position': 'faldo:position',
        'has_measurement': 'IAO:0000004',
        'has_quantifier': 'GENO:0000866',
        'created_on': 'pav:createdOn'
    }

    properties = annotation_properties.copy()
    properties.update(object_properties)
    properties.update(datatype_properties)

    def __init__(self, graph, definedby, sub=None, obj=None, pred=None):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)

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

    def get_properties(self):
        return self.properties

    def _is_valid(self):

        # check if sub/obj/rel are none...throw error
        if self.sub is None:
            raise ValueError('No subject set for this association')
        if self.obj is None:
            raise ValueError('No object set for this association')
        if self.rel is None:
            raise ValueError('No relation set for this association')

        return True

    def _add_basic_association_to_graph(self):

        if not self._is_valid():
            return

        self.graph.addTriple(self.sub, self.rel, self.obj)

        if self.assoc_id is None:
            self.set_association_id()

        self.model.addType(self.assoc_id, self.assoc_types['association'])

        self.graph.addTriple(
            self.assoc_id, self.object_properties['has_subject'], self.sub)
        self.graph.addTriple(
            self.assoc_id, self.object_properties['has_object'], self.obj)
        self.graph.addTriple(
            self.assoc_id, self.object_properties['has_predicate'], self.rel)

        if self.description is not None:
            self.model.addDescription(self.assoc_id, self.description)

        if self.evidence is not None and len(self.evidence) > 0:
            for e in self.evidence:
                self.graph.addTriple(
                    self.assoc_id, self.object_properties['has_evidence'], e)

        if self.source is not None and len(self.source) > 0:
            for s in self.source:
                if re.match('http', s):
                    # TODO assume that the source is a publication?
                    # use Reference class here
                    self.graph.addTriple(
                        self.assoc_id, self.object_properties['has_source'],
                        s, True)
                else:
                    self.graph.addTriple(
                        self.assoc_id, self.object_properties['has_source'], s)

        if self.provenance is not None and len(self.provenance) > 0:
            for p in self.provenance:
                self.graph.addTriple(
                    self.assoc_id, self.object_properties['has_provenance'], p)

        if self.date is not None and len(self.date) > 0:
            for d in self.date:
                self.graph.addTriple(
                    object_is_literal=True,
                    subject_id=self.assoc_id,
                    predicate_id=self.datatype_properties['created_on'],
                    obj=d)

        if self.score is not None:
            self.graph.addTriple(
                self.assoc_id, self.properties['has_measurement'],
                self.score, True, 'xsd:float')
            # TODO
            # update with some kind of instance of scoring object
            # that has a unit and type

        return

    def add_association_to_graph(self):

        self._add_basic_association_to_graph()

        return

    def add_predicate_object(self, predicate, object_node,
                             object_type=None, datatype=None):

        if object_type == 'Literal':
            if datatype is not None:
                self.graph.addTriple(self.assoc_id, predicate,
                                     object_node, True, datatype)
            else:
                self.graph.addTriple(self.assoc_id, predicate,
                                     object_node, True)
        else:
            self.graph.addTriple(self.assoc_id, predicate,
                                 object_node, False)

        return

    # This isn't java, but if we must,
    # prefer use of property decorator
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
            self.assoc_id = self.make_association_id(self.definedby, self.sub,
                                                     self.rel, self.obj)
        else:
            self.assoc_id = assoc_id

        return

    def get_association_id(self):

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
    def make_association_id(definedby, subject, predicate, object,
                            attributes=None):
        """
        A method to create unique identifiers for OBAN-style associations,
        based on all the parts of the association
        If any of the items is empty or None, it will convert it to blank.
        It effectively digests the  string of concatonated values.
        Subclasses of Assoc can submit an additional array of attributes
        that will be appeded to the ID.

        :param definedby: The (data) resource that provided the annotation
        :param subject:
        :param predicate:
        :param object:
        :param attributes:

        :return:

        """

        # note others available:
        #   md5(), sha1(), sha224(), sha256(), sha384(), and sha512()
        # putting definedby first,
        # as this will usually be the datasource providing the annotation
        # this will end up making the first few parts of the id
        # be the same for all annotations in that resource
        # (although the point of a digest is to render such details moot).

        items_to_hash = [definedby, subject, predicate, object]
        if attributes is not None:
            items_to_hash += attributes

        for i, val in enumerate(items_to_hash):
            if val is None:
                items_to_hash[i] = ''

        byte_string = '+'.join(items_to_hash).encode("utf-8")

        # TODO put this in a util?
        return ':'.join(('MONARCH', hashlib.sha1(byte_string).hexdigest()[0:16]))
