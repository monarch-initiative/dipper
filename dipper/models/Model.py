import logging
import re
from dipper.graph.Graph import Graph

logger = logging.getLogger(__name__)


class Model():
    """
    Utility class to add common triples to a graph
    (subClassOf, type, label, sameAs)
    """

    types = {
        'ontology': 'owl:Ontology',
        'class': 'owl:Class',
        'named_individual': 'owl:NamedIndividual',
        'object_property': 'owl:ObjectProperty',
        'annotation_property': 'owl:AnnotationProperty',
        'datatype_property': 'owl:DatatypeProperty',
        'restriction': 'owl:Restriction',
        'deprecated': 'owl:deprecated',
        'person': 'foaf:Person'
    }

    annotation_properties = {
        'label': 'rdfs:label',
        'description': 'dc:description',
        'replaced_by': 'IAO:0100001',
        'consider': 'OIO:consider',
        'depiction': 'foaf:depiction',
        'comment': 'dc:comment',
        'version_info': 'owl:versionInfo',
        'hasExactSynonym': 'OIO:hasExactSynonym',
        'hasRelatedSynonym': 'OIO:hasRelatedSynonym',
        'definition': 'IAO:0000115',
        'has_xref': 'OIO:hasDbXref',
        'clique_leader': 'MONARCH:cliqueLeader',
        'inchi_key': 'CHEBI:InChIKey',
        'is_anonymous': 'MONARCH:anonymous'
    }

    object_properties = {
        'type': 'rdf:type',
        'same_as': 'owl:sameAs',
        'equivalent_class': 'owl:equivalentClass',
        'subclass_of': 'rdfs:subClassOf',
        'on_property': 'owl:onProperty',
        'some_values_from': 'owl:someValuesFrom',
        'version_iri': 'owl:versionIRI',
        'has disposition': 'RO:0000091',
        'has_phenotype': 'RO:0002200',
        'in_taxon': 'RO:0002162',
        'has_quality': 'RO:0000086',
        'has_qualifier': 'GENO:0000580',
        'towards': 'RO:0002503',
        'has_subject': ':hasSubject',
        'has_object': ':hasObject',
        'has_predicate': ':hasPredicate',
        'is_about': 'IAO:0000136',
        'involved_in': 'RO:0002331',
        'enables': 'RO:0002327',
        'derives_from': 'RO:0001000',
        'part_of': 'BFO:0000050',
        'has_part': 'BFO:0000051',
        'mentions': 'IAO:0000142',
        'model_of': 'RO:0003301',
        'has_gene_product': 'RO:0002205',
        'existence_starts_at': 'UBERON:existence_starts_at',
        'existence_starts_during': 'RO:0002488',
        'existence_ends_at': 'UBERON:existence_ends_at',
        'existence_ends_during': 'RO:0002492',
        'starts_with': 'RO:0002224',
        'starts_during': 'RO:0002091',
        'ends_during': 'RO:0002093',
        'ends_with': 'RO:0002230',
        'occurs_in': 'BFO:0000066',
        'has_environment_qualifier': 'GENO:0000580',
        'has_begin_stage_qualifier': 'GENO:0000630',
        'has_end_stage_qualifier': 'GENO:0000631',
        'correlates_with': 'RO:0002610',
        'substance_that_treats': 'RO:0002606',
        'is_marker_for': 'RO:0002607',
        'contributes_to': 'RO:0002326',
        'has_origin': 'GENO:0000643',
        'has_author': 'ERO:0000232',
        'dc:source': 'dc:source',
        'dc:evidence': 'dc:evidence',
        'has_evidence': 'RO:0002558',
        'causally_upstream_of_or_within': 'RO:0002418',
        'causally_influences': 'RO:0002566'
    }

    datatype_properties = {
        'position': 'faldo:position',
        'has_measurement': 'IAO:0000004'
    }

    def __init__(self, graph):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)

    def addTriple(self, subject_id, predicate_id, obj,
                  object_is_literal=False, literal_type=None):
        self.graph.addTriple(subject_id, predicate_id, obj,
                             object_is_literal, literal_type)

    def addType(self, subject_id, subject_type):
        self.graph.addTriple(subject_id, self.object_properties['type'],
                             subject_type)
        return

    def addLabel(self, subject_id, label):
        self.graph.addTriple(subject_id, self.annotation_properties['label'],
                             label, object_is_literal=True)
        return

    def addClassToGraph(self, class_id, label=None,
                        class_type=None, description=None):
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

        self.graph.addTriple(class_id, self.object_properties['type'],
                             self.types['class'])
        if label is not None:
            self.graph.addTriple(
                class_id, self.annotation_properties['label'], label,
                object_is_literal=True)

        if class_type is not None:
            self.graph.addTriple(class_id, self.object_properties['subclass_of'], class_type)
        if description is not None:
            self.graph.addTriple(
                class_id, self.annotation_properties['description'], description,
                object_is_literal=True)
        return

    def addIndividualToGraph(self, ind_id, label, ind_type=None, description=None):

        if label is not None:
            self.graph.addTriple(
                ind_id, self.annotation_properties['label'], label,
                object_is_literal=True
            )
        if ind_type is not None:
            self.graph.addTriple(
                ind_id, self.object_properties['type'], ind_type
            )
        else:
            self.graph.addTriple(
                ind_id, self.object_properties['type'],
                self.types['named_individual']
            )
        if description is not None:
            self.graph.addTriple(
                ind_id, self.annotation_properties['description'], description,
                object_is_literal=True)
        return

    def addEquivalentClass(self, sub, obj):
        self.graph.addTriple(sub, self.object_properties['equivalent_class'],
                             obj)
        return

    def addSameIndividual(self, sub, obj):
        self.graph.addTriple(sub, self.object_properties['same_as'], obj)

        return

    def addOWLPropertyClassRestriction(
            self, class_id, property_id, property_value):

        # make a blank node to hold the property restrictions
        # scrub the colons, they will make the ttl parsers choke
        bnode = \
            '_:'+re.sub(r':', '', property_id)+re.sub(r':', '', property_value)

        self.graph.addTriple(bnode, self.object_properties['type'],
                             self.types['restriction'])
        self.graph.addTriple(bnode, self.object_properties['on_property'],
                             property_id)
        self.graph.addTriple(bnode, self.object_properties['some_values_from'],
                             property_value)
        self.graph.addTriple(class_id, self.object_properties['subclass_of'],
                             bnode)

        return

    def addPerson(self, person_id, person_label=None):
        self.graph.addTriple(person_id, self.object_properties['type'],
                             self.types['person'])
        if person_label is not None:
            self.graph.addTriple(
                person_id, self.annotation_properties['label'],
                person_label, object_is_literal=True
            )
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
        self.graph.addTriple(old_id, self.object_properties['type'],
                             self.types['class'])

        self._addReplacementIds(old_id, new_ids)

        return

    def _addReplacementIds(self, old_id, new_ids):
        consider = self.annotation_properties['consider']
        replaced_by = self.annotation_properties['replaced_by']

        self.graph.addTriple(old_id, self.types['deprecated'],
                             True, object_is_literal=True,
                             literal_type='xsd:boolean')

        if new_ids is not None:
            if isinstance(new_ids, str):
                self.graph.addTriple(old_id, replaced_by, new_ids)
            elif len(new_ids) == 1:
                self.graph.addTriple(old_id, replaced_by, new_ids[0])
            elif len(new_ids) > 0:
                for new_id in new_ids:
                    self.graph.addTriple(old_id, consider, new_id)
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
        self.graph.addTriple(old_id, self.object_properties['type'],
                             self.types['named_individual'])

        self._addReplacementIds(old_id, new_ids)

        return

    def addSubClass(self, child_id, parent_id):
        self.graph.addTriple(child_id,
                             self.object_properties['subclass_of'],
                             parent_id)
        return


    def addSynonym(self, class_id, synonym, synonym_type=None):
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
            # default
            synonym_type = self.annotation_properties['hasExactSynonym']

        if synonym is not None:
            self.graph.addTriple(class_id, synonym_type, synonym,
                                 object_is_literal=True)
        return

    def addDefinition(self, class_id, definition):
        self.graph.addTriple(class_id, self.annotation_properties['definition'],
                             definition, object_is_literal=True)
        return

    def addXref(self, class_id, xref_id, xref_as_literal=False):
        self.graph.addTriple(
            class_id, self.annotation_properties['has_xref'], xref_id,
            object_is_literal=xref_as_literal)
        return

    def addDepiction(self, subject_id, image_url):
        self.graph.addTriple(
            subject_id, self.annotation_properties['depiction'],
            image_url, object_is_literal=True)
        return

    def addComment(self, subject_id, comment):
        self.graph.addTriple(
            subject_id, self.annotation_properties['comment'],
            comment.strip(), object_is_literal=True)
        return

    def addDescription(self, subject_id, description):
        self.graph.addTriple(
            subject_id, self.annotation_properties['description'],
            description.strip(), object_is_literal=True)
        return

    def addOntologyDeclaration(self, ontology_id):

        self.graph.addTriple(ontology_id, self.object_properties['type'],
                             self.types['ontology'])
        return

    def addOWLVersionIRI(self, ontology_id, version_iri):
        self.graph.addTriple(
            ontology_id, self.object_properties['version_iri'], version_iri)

        return

    def addOWLVersionInfo(self, ontology_id, version_info):
        self.graph.addTriple(
            ontology_id, self.annotation_properties['version_info'],
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
            node_id, self.annotation_properties['clique_leader'], True,
            object_is_literal=True, literal_type='xsd:boolean')
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
            node_id, self.annotation_properties['is_anonymous'], True,
            object_is_literal=True, literal_type='xsd:boolean')
        return
