import re
import logging
from dipper.graph.Graph import Graph
from dipper.models.Model import Model

__author__ = 'nlw'

logger = logging.getLogger(__name__)


class Reference:
    """
    To model references for associations
        (such as journal articles, books, etc.).

    By default, references will be typed as "documents",
        unless if the type is set otherwise.

    If a short_citation is set, this will be used for the individual's label.
        We may wish to subclass this later.

    """

    ref_types = {
        'person': 'foaf:Person',
        'journal_article': 'IAO:0000013',
        'publication': 'IAO:0000311',  # book
        'document': 'IAO:0000310',  # document???
        'photograph': 'IAO:0000185',
        'webpage': 'SIO:000302',
    }

    annotation_properties = {
        'page': 'foaf:page',
        'title': 'dc:title'
    }

    def __init__(self, graph, ref_id=None, ref_type=None):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.ref_id = ref_id
        self.ref_url = None
        self.title = None
        self.year = None
        self.author_list = None
        self.short_citation = None

        self.model = Model(self.graph)

        if ref_type is None:
            self.ref_type = self.ref_types['document']
        else:
            self.ref_type = ref_type

        if ref_id is not None and re.match(r'http', ref_id):
            self.ref_url = ref_id

        return

    def setTitle(self, title):
        self.title = title
        return

    def setYear(self, year):

        self.year = year

        return

    def setType(self, reference_type):

        self.ref_type = reference_type

        return

    def setAuthorList(self, author_list):
        """

        :param author_list: Array of authors
        :return:
        """

        self.author_list = author_list
        return

    def addAuthor(self, author):

        self.author_list += [author]

        return

    def setShortCitation(self, citation):
        self.short_citation = citation
        return

    def addPage(self, subject_id, page_url):
        self.graph.addTriple(
            subject_id, self.annotation_properties['page'],
            page_url, object_is_literal=True)
        return

    def addTitle(self, subject_id, title):
        self.graph.addTriple(
            subject_id, self.annotation_properties['title'],
            title, object_is_literal=True)
        return

    def addRefToGraph(self):

        n = self.short_citation
        if n is None:
            n = self.title

        if self.ref_url is not None:
            self.addTitle(self.ref_url, self.title)
            self.model.addType(self.ref_url, self.ref_type)
            self.model.addLabel(self.ref_url, n)
        elif self.ref_id is not None:
            self.model.addIndividualToGraph(self.ref_id, n, self.ref_type)
            if self.title is not None:
                self.addTitle(self.ref_id, self.title)
        else:
            # should never be true
            logger.error("You are missing an identifier for a reference.")

        # TODO what is the property here to add the date?
        # if self.year is not None:
        #    gu.addTriple()

        # if self.author_list is not None:
        #    for a in self.author_list:
        #        gu.addTriple(
        #           g, self.ref_id, self.props['has_author'], a, True)
        return
