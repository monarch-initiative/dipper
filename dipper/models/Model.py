import logging
import re
import yaml
from dipper.graph.Graph import Graph

logger = logging.getLogger(__name__)


class Model():
    """
    Utility class to add common triples to a graph
    (subClassOf, type, label, sameAs)
    """
    # shared ontology_label to ontology_id mapping
    globaltt = {}

    # reversed to find labels for items that only came with IDs; such as some taxons
    globaltcid = {}

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
            self.globaltt = self.graph.globaltt

        else:
            raise ValueError("{} is not a graph".graph)

    def addTriple(
        self, subject_id, predicate_id, obj, object_is_literal=False, literal_type=None
    ):
        self.graph.addTriple(
            subject_id, predicate_id, obj, object_is_literal, literal_type)

    def addType(self, subject_id, subject_type):
        self.graph.addTriple(
            subject_id, self.globaltt['type'], subject_type)
        return

    def addLabel(self, subject_id, label):
        self.graph.addTriple(
            subject_id, self.globaltt['label'], label, object_is_literal=True)
        return

    def addClassToGraph(
        self, class_id, label=None, class_type=None, description=None
    ):
        """
        Any node added to the graph will get at least 3 triples:
        *(node, type, owl:Class) and
        *(node, label, literal(label))
        *if a type is added,
            then the node will be an OWL:subclassOf that the type
        *if a description is provided,
            it will also get added as a dc:description
        :param class_id:
        :param label:
        :param class_type:
        :param description:
        :return:

        """

        self.graph.addTriple(
            class_id, self.globaltt['type'], self.globaltt['class'])
        if label is not None:
            self.graph.addTriple(
                class_id, self.globaltt['label'], label, object_is_literal=True)

        if class_type is not None:
            self.graph.addTriple(
                class_id, self.globaltt['subclass_of'], class_type)
        if description is not None:
            self.graph.addTriple(
                class_id, self.globaltt['description'], description,
                object_is_literal=True)
        return

    def addIndividualToGraph(self, ind_id, label, ind_type=None, description=None):

        if label is not None:
            self.graph.addTriple(
                ind_id, self.globaltt['label'], label,
                object_is_literal=True
            )
        if ind_type is not None:
            self.graph.addTriple(
                ind_id, self.globaltt['type'], ind_type
            )
        else:
            self.graph.addTriple(
                ind_id, self.globaltt['type'],
                self.globaltt['named_individual']
            )
        if description is not None:
            self.graph.addTriple(
                ind_id, self.globaltt['description'], description,
                object_is_literal=True)
        return

    def addEquivalentClass(self, sub, obj):
        self.graph.addTriple(
            sub, self.globaltt['equivalent_class'], obj)
        return

    def addSameIndividual(self, sub, obj):
        self.graph.addTriple(sub, self.globaltt['same_as'], obj)

        return

    def addOWLPropertyClassRestriction(self, class_id, property_id, property_value):

        # make a blank node to hold the property restrictions
        # scrub the colons, they will make the ttl parsers choke
        bnode = '_:'+re.sub(
            r':', '', property_id)+re.sub(r':', '', property_value)

        self.graph.addTriple(
            bnode, self.globaltt['type'], self.globaltt['restriction'])
        self.graph.addTriple(
            bnode, self.globaltt['on_property'], property_id)
        self.graph.addTriple(
            bnode, self.globaltt['some_values_from'], property_value)
        self.graph.addTriple(
            class_id, self.globaltt['subclass_of'], bnode)

        return

    def addPerson(self, person_id, person_label=None):
        self.graph.addTriple(
            person_id, self.globaltt['type'], self.globaltt['person'])
        if person_label is not None:
            self.graph.addTriple(
                person_id, self.globaltt['label'], person_label, object_is_literal=True)
        return

    def addDeprecatedClass(self, old_id, new_ids=None):
        """
        Will mark the oldid as a deprecated class.
        if one newid is supplied, it will mark it as replaced by.
        if >1 newid is supplied, it will mark it with consider properties
        :param old_id: str - the class id to deprecate
        :param new_ids: list - the class list that is
                       the replacement(s) of the old class.  Not required.
        :return: None

        """
        self.graph.addTriple(
            old_id, self.globaltt['type'], self.globaltt['class'])

        self._addReplacementIds(old_id, new_ids)

        return

    def _addReplacementIds(self, old_id, new_ids):

        self.graph.addTriple(
            old_id, self.globaltt['deprecated'], True, object_is_literal=True,
            literal_type='xsd:boolean')

        if new_ids is not None:
            if isinstance(new_ids, str):
                self.graph.addTriple(old_id, self.globaltt['term replaced by'], new_ids)
            elif len(new_ids) == 1:
                self.graph.addTriple(
                    old_id, self.globaltt['term replaced by'], new_ids[0])
            elif len(new_ids) > 0:
                for new_id in new_ids:
                    self.graph.addTriple(old_id, self.globaltt['consider'], new_id)
        return

    def addDeprecatedIndividual(self, old_id, new_ids=None):
        """
        Will mark the oldid as a deprecated individual.
        if one newid is supplied, it will mark it as replaced by.
        if >1 newid is supplied, it will mark it with consider properties
        :param g:
        :param oldid: the individual id to deprecate
        :param newids: the individual idlist that is the replacement(s) of
                       the old individual.  Not required.
        :return:

        """
        self.graph.addTriple(
            old_id, self.globaltt['type'], self.globaltt['named_individual'])

        self._addReplacementIds(old_id, new_ids)

        return

    def addSubClass(self, child_id, parent_id):
        self.graph.addTriple(
            child_id, self.globaltt['subclass_of'], parent_id)
        return

    def addSynonym(
            self, class_id, synonym, synonym_type=None):
        """
        Add the synonym as a property of the class cid.
        Assume it is an exact synonym, unless otherwise specified
        :param g:
        :param cid: class id
        :param synonym: the literal synonym label
        :param synonym_type: the CURIE of the synonym type (not the URI)
        :return:

        """

        if synonym_type is None:
            synonym_type = self.globaltt['hasExactSynonym']

        if synonym is not None:
            self.graph.addTriple(
                class_id, synonym_type, synonym, object_is_literal=True)
        return

    def addDefinition(self, class_id, definition):
        self.graph.addTriple(
            class_id, self.globaltt['definition'],  definition, object_is_literal=True)
        return

    def addXref(self, class_id, xref_id, xref_as_literal=False):
        self.graph.addTriple(
            class_id, self.globaltt['has_dbxref'], xref_id,
            object_is_literal=xref_as_literal)
        return

    def addDepiction(self, subject_id, image_url):
        self.graph.addTriple(
            subject_id, self.globaltt['depiction'], image_url, object_is_literal=True)
        return

    def addComment(self, subject_id, comment):
        self.graph.addTriple(
            subject_id, self.globaltt['comment'], comment.strip(),
            object_is_literal=True)
        return

    def addDescription(self, subject_id, description):
        self.graph.addTriple(
            subject_id, self.globaltt['description'], description.strip(),
            object_is_literal=True)
        return

    def addOntologyDeclaration(self, ontology_id):
        self.graph.addTriple(
            ontology_id, self.globaltt['type'], self.globaltt['ontology'])
        return

    def addOWLVersionIRI(self, ontology_id, version_iri):
        self.graph.addTriple(ontology_id, self.globaltt['version_iri'], version_iri)

        return

    def addOWLVersionInfo(self, ontology_id, version_info):
        self.graph.addTriple(
            ontology_id, self.globaltt['version_info'],
            version_info, object_is_literal=True)
        return

    def makeLeader(self, node_id):
        """
        Add an annotation property to the given ```node_id```
        to be the clique_leader.
        This is a monarchism.
        :param node_id:
        :return:
        """
        self.graph.addTriple(
            node_id, self.globaltt['clique_leader'], True, object_is_literal=True,
            literal_type='xsd:boolean')
        return

    def addBlankNodeAnnotation(self, node_id):
        """
        Add an annotation property to the given ```node_id```
        to be a pseudo blank node.
        This is a monarchism.
        :param node_id:
        :return:
        """
        self.graph.addTriple(
            node_id, self.globaltt['is_anonymous'], True, object_is_literal=True,
            literal_type='xsd:boolean')
        return

    def _addSexSpecificity(self, subject_id, sex):
        """
        Add sex specificity to a subject (eg association node)

        In our modeling we use this to add a qualifier to a triple
        for example, this genotype to phenotype association
        is specific to this sex (see MGI, IMPC)

        This expects the client to define the ontology term
        for sex (eg PATO)

        Note this class is probably not the right place for this
        method, but putting here until a better home is found
        :param subject_id:
        :param sex:
        :return:
        """
        self.graph.addTriple(subject_id, self.globaltt['has_sex_specificty'], sex)
        return
