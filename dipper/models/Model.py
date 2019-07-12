import logging
import re
from dipper.graph.Graph import Graph
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)
# note: currently no log issued


class Model():
    """
    Utility class to add common triples to a graph
    (subClassOf, type, label, sameAs)
    """

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
            self.globaltt = self.graph.globaltt
            self.globaltcid = self.graph.globaltcid
            self.curie_map = self.graph.curie_map

        else:
            raise ValueError("{} is not a graph".format(graph))

    def addTriple(
            self, subject_id, predicate_id, obj, object_is_literal=False,
            literal_type=None,
            subject_category=None, object_category=None):
        self.graph.addTriple(
            subject_id, predicate_id, obj, object_is_literal, literal_type,
            subject_category=subject_category, object_category=object_category)

    def addType(self, subject_id, subject_type,
                subject_category=None,
                subject_type_category=None):
        self.graph.addTriple(
            subject_id, self.globaltt['type'], subject_type,
            subject_category=subject_category,
            object_category=subject_type_category)

    def addLabel(self, subject_id, label, **args):
        self.graph.addTriple(
            subject_id, self.globaltt['label'], label, object_is_literal=True, **args)

    def addClassToGraph(
            self, class_id, label=None, class_type=None, description=None,
            class_category=None, class_type_category=None):
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
        :param class_category: a biolink category CURIE for class
        :param class_type_category: a biolink category CURIE for class type
        :return:

        """
        if class_id is None:
            raise ValueError("class_id is None")

        self.graph.addTriple(
            class_id, self.globaltt['type'], self.globaltt['class'],
            subject_category=class_category)
        if label is not None:
            self.graph.addTriple(
                class_id, self.globaltt['label'], label, object_is_literal=True)

        if class_type is not None:
            self.graph.addTriple(class_id, self.globaltt['subclass_of'], class_type,
                                 object_category=class_type_category)
        if description is not None:
            self.graph.addTriple(
                class_id, self.globaltt['description'], description,
                object_is_literal=True)

    def addIndividualToGraph(self, ind_id, label, ind_type=None, description=None,
                             ind_category=None,  # blv category for ind_id
                             ind_type_category=None):
        if label is not None:
            self.graph.addTriple(
                ind_id, self.globaltt['label'], label, object_is_literal=True)
        if ind_type is not None:
            self.graph.addTriple(
                ind_id, self.globaltt['type'], ind_type, object_is_literal=False,
                subject_category=ind_category, object_category=ind_type_category)
        else:
            self.graph.addTriple(
                ind_id, self.globaltt['type'], self.globaltt['named_individual'],
                subject_category=ind_category)
        if description is not None:
            self.graph.addTriple(
                ind_id, self.globaltt['description'], description,
                object_is_literal=True)

    def addEquivalentClass(self, sub, obj, subject_category=None, object_category=None,
                           **args):
        self.graph.addTriple(
            sub, self.globaltt['equivalent_class'], obj,
            subject_category=subject_category,
            object_category=object_category,
            **args)

    def addSameIndividual(self, sub, obj, subject_category=None, object_category=None,
                          **args):
        self.graph.addTriple(sub, self.globaltt['same_as'], obj,
                             subject_category=subject_category,
                             object_category=object_category,
                             **args)

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
            person_id, self.globaltt['type'], self.globaltt['person'],
            subject_category=blv.Case.value)
        if person_label is not None:
            self.graph.addTriple(
                person_id, self.globaltt['label'], person_label, object_is_literal=True)

    def addDeprecatedClass(self, old_id, new_ids=None,
                           old_id_category=None, new_ids_category=None):
        """
        Will mark the oldid as a deprecated class.
        if one newid is supplied, it will mark it as replaced by.
        if >1 newid is supplied, it will mark it with consider properties
        :param old_id: str - the class id to deprecate
        :param new_ids: list - the class list that is
                       the replacement(s) of the old class.  Not required.
        :param old_id_category - a biolink category CURIE for old id
        :param new_ids_category - a biolink category CURIE for new ids
        :return: None

        """
        self.graph.addTriple(
            old_id, self.globaltt['type'], self.globaltt['class'],
            subject_category=old_id_category)

        self._addReplacementIds(old_id, new_ids, new_ids_category=new_ids_category)

    def _addReplacementIds(self, old_id, new_ids, new_ids_category=None):

        self.graph.addTriple(
            old_id, self.globaltt['deprecated'], True, object_is_literal=True,
            literal_type='xsd:boolean')

        if new_ids is not None:
            if isinstance(new_ids, str):
                self.graph.addTriple(old_id, self.globaltt['term replaced by'], new_ids)
            elif len(new_ids) == 1:
                self.graph.addTriple(
                    old_id, self.globaltt['term replaced by'], new_ids[0],
                    object_category=new_ids_category)
            elif new_ids:
                for new_id in new_ids:
                    self.graph.addTriple(old_id, self.globaltt['consider'], new_id,
                                         object_category=new_ids_category)

    def addDeprecatedIndividual(self, old_id, new_ids=None,
                                old_id_category=None, new_ids_category=None):
        """
        Will mark the oldid as a deprecated individual.
        if one newid is supplied, it will mark it as replaced by.
        if >1 newid is supplied, it will mark it with consider properties
        :param g:
        :param oldid: the individual id to deprecate
        :param newids: the individual idlist that is the replacement(s) of
                       the old individual.  Not required.
        :param old_id_category - a biolink category CURIE for old id
        :param new_ids_category - a biolink category CURIE for new ids
        :return:

        """
        self.graph.addTriple(
            old_id, self.globaltt['type'], self.globaltt['named_individual'],
            subject_category=old_id_category)

        self._addReplacementIds(old_id, new_ids, new_ids_category=new_ids_category)

    def addSubClass(self, child_id, parent_id, child_category=None,
                    parent_category=None):
        self.graph.addTriple(child_id, self.globaltt['subclass_of'], parent_id,
                             subject_category=child_category,
                             object_category=parent_category)
        
    def addSynonym(
            self, class_id, synonym, synonym_type=None,
            class_category=None, synonym_category=None):
        """
        Add the synonym as a property of the class cid.
        Assume it is an exact synonym, unless otherwise specified
        :param self:
        :param class_id: class id
        :param synonym: the literal synonym label
        :param synonym_type: the CURIE of the synonym type (not the URI)
        :param class_category: biolink category CURIE for class_id
        :param synonym_category = biolink category CURIE for synonym
        :return:

        """
        if synonym_type is None:
            synonym_type = self.globaltt['has_exact_synonym']

        if synonym is not None:
            self.graph.addTriple(
                class_id, synonym_type, synonym, object_is_literal=True,
                subject_category=class_category,
                object_category=synonym_category)

    def addDefinition(self, class_id, definition):
        self.graph.addTriple(
            class_id, self.globaltt['definition'], definition, object_is_literal=True)

    def addXref(self, class_id, xref_id, xref_as_literal=False,
                class_category=None, xref_category=None):
        self.graph.addTriple(
            class_id, self.globaltt['database_cross_reference'], xref_id,
            subject_category=class_category, object_category=xref_category,
            object_is_literal=xref_as_literal)

    def addDepiction(self, subject_id, image_url, subject_category=None,
                     object_category=blv.InformationContentEntity.value):
        self.graph.addTriple(
            subject_id, self.globaltt['depiction'], image_url, object_is_literal=False,
            subject_category=subject_category, object_category=object_category)

    def addComment(self, subject_id, comment, subject_category=None):
        self.graph.addTriple(
            subject_id, self.globaltt['comment'], comment.strip(),
            object_is_literal=True, subject_category=subject_category
        )

    def addDescription(self, subject_id, description, subject_category=None):
        self.graph.addTriple(
            subject_id, self.globaltt['description'], description.strip(),
            object_is_literal=True, subject_category=subject_category)

    def addOntologyDeclaration(self, ontology_id):
        self.graph.addTriple(
            ontology_id, self.globaltt['type'], self.globaltt['ontology'],
            subject_category=blv.OntologyClass.value)

    def addOWLVersionIRI(self, ontology_id, version_iri):
        self.graph.addTriple(ontology_id, self.globaltt['version_iri'], version_iri,
                             subject_category=blv.OntologyClass.value)

    def addOWLVersionInfo(self, ontology_id, version_info):
        self.graph.addTriple(
            ontology_id, self.globaltt['version_info'],
            version_info, object_is_literal=True)

    def makeLeader(self, node_id, node_category=None):
        """
        Add an annotation property to the given ```node_id```
        to be the clique_leader.
        This is a monarchism.
        :param node_id:
        :param node_category: a biolink category CURIE for node_id
        :return:
        """
        self.graph.addTriple(
            node_id, self.globaltt['clique_leader'], True, object_is_literal=True,
            literal_type='xsd:boolean', subject_category=node_category)

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

    def _addSexSpecificity(self, subject_id, sex,
                           subject_category=None):
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
        :param subject_category: a biolink category CURIE for subject_id
        :param sex:
        :return:
        """
        self.graph.addTriple(subject_id, self.globaltt['has_sex_specificty'], sex,
                             subject_category=subject_category,
                             object_category=blv.BiologicalSex.value)
