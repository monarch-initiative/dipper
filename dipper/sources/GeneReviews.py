import re
import os
import csv
import logging
from bs4 import BeautifulSoup

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Model import Model
from dipper.sources.OMIM import OMIM, filter_keep_phenotype_entry_ids
from dipper import config
from dipper.models.Reference import Reference

__author__ = 'nicole'

logger = logging.getLogger(__name__)
GRDL = 'http://ftp.ncbi.nih.gov/pub/GeneReviews'


class GeneReviews(Source):
    """
    Here we process the GeneReviews mappings to OMIM,
    plus inspect the GeneReviews (html) books to pull the clinical descriptions
    in order to populate the definitions of the terms in the ontology.
    We define the GeneReviews items as classes that are either grouping classes
    over OMIM disease ids (gene ids are filtered out),
    or are made as subclasses of DOID:4 (generic disease).

    Note that GeneReviews
    [copyright policy](http://www.ncbi.nlm.nih.gov/books/NBK138602/)
    (as of 2015.11.20) says:

    GeneReviews® chapters are owned by the University of Washington, Seattle,
    © 1993-2015. Permission is hereby granted to reproduce, distribute,
    and translate copies of content materials provided that
    (i) credit for source (www.ncbi.nlm.nih.gov/books/NBK1116/)
    and copyright (University of Washington, Seattle)
    are included with each copy;
    (ii) a link to the original material is provided whenever the material is
    published elsewhere on the Web; and
    (iii) reproducers, distributors, and/or translators comply with this
    copyright notice and the GeneReviews Usage Disclaimer.

    This script doesn't pull the GeneReviews books from the NCBI Bookshelf
    directly; scripting this task is expressly prohibited by
    [NCBIBookshelf policy](http://www.ncbi.nlm.nih.gov/books/NBK45311/).
    However, assuming you have acquired the books (in html format) via
    permissible means, a parser for those books is provided here to extract
    the clinical descriptions to define the NBK identified classes.

    """

    files = {
        'idmap': {'file': 'NBKid_shortname_OMIM.txt',
                  'url': GRDL + '/NBKid_shortname_OMIM.txt'},
        'titles': {'file': 'GRtitle_shortname_NBKid.txt',
                   'url': GRDL + '/GRtitle_shortname_NBKid.txt'}
        }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'genereviews')

        self.dataset = Dataset(
            'genereviews', 'Gene Reviews', 'http://genereviews.org/',
            None, 'http://www.ncbi.nlm.nih.gov/books/NBK138602/')
        self.dataset.set_citation('GeneReviews:NBK1116')

        self.book_ids = set()
        self.all_books = {}

        if 'test_ids' not in config.get_config() or\
                'disease' not in config.get_config()['test_ids']:
            logger.warning("not configured with disease test ids.")
            self.test_ids = list()
        else:
            # select ony those test ids that are omim's.
            self.test_ids = config.get_config()['test_ids']['disease']

        return

    def fetch(self, is_dl_forced=False):
        """
        We fetch GeneReviews id-label map and id-omim mapping files from NCBI.
        :return: None
        """

        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        """
        :return: None
        """

        if self.testOnly:
            self.testMode = True

        self._get_titles(limit)
        self._get_equivids(limit)

        self.create_books()
        self.process_nbk_html(limit)

        # no test subset for now; test == full graph
        self.testgraph = self.graph
        return

    def _get_equivids(self, limit):
        """
        The file processed here is of the format:
        #NBK_id GR_shortname    OMIM
        NBK1103 trimethylaminuria       136132
        NBK1103 trimethylaminuria       602079
        NBK1104 cdls    122470
        Where each of the rows represents a mapping between
        a gr id and an omim id. These are a 1:many relationship,
        and some of the omim ids are genes(not diseases).
        Therefore, we need to create a loose coupling here.
        We make the assumption that these NBKs are generally higher-level
        grouping classes; therefore the OMIM ids are treated as subclasses.
        (This assumption is poor for those omims that are actually genes,
        but we have no way of knowing what those are here...
        we will just have to deal with that for now.)
        :param limit:
        :return:

        """
        raw = '/'.join((self.rawdir, self.files['idmap']['file']))
        model = Model(self.graph)
        line_counter = 0

        # we look some stuff up in OMIM, so initialize here
        omim = OMIM(self.graph_type, self.are_bnodes_skized)
        id_map = {}
        allomimids = set()
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                if line_counter == 1:  # skip header
                    continue
                (nbk_num, shortname, omim_num) = row
                gr_id = 'GeneReviews:'+nbk_num
                omim_id = 'OMIM:'+omim_num
                if not (
                        (self.testMode and
                         len(self.test_ids) > 0 and
                         omim_id in self.test_ids) or not
                        self.testMode):
                    continue

                # sometimes there's bad omim nums
                if len(omim_num) > 6:
                    logger.warning(
                        "OMIM number incorrectly formatted " +
                        "in row %d; skipping:\n%s",
                        line_counter, '\t'.join(row))
                    continue

                # build up a hashmap of the mappings; then process later
                if nbk_num not in id_map:
                    id_map[nbk_num] = set()
                id_map[nbk_num].add(omim_num)

                # add the class along with the shortname
                model.addClassToGraph(gr_id, None)
                model.addSynonym(gr_id, shortname)

                allomimids.add(omim_num)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

            # end looping through file

        # get the omim ids that are not genes
        entries_that_are_phenotypes = \
            omim.process_entries(
                list(allomimids), filter_keep_phenotype_entry_ids,
                None, None, limit)

        logger.info("Filtered out %d/%d entries that are genes or features",
                    len(allomimids)-len(entries_that_are_phenotypes),
                    len(allomimids))

        for nbk_num in self.book_ids:
            gr_id = 'GeneReviews:'+nbk_num
            if nbk_num in id_map:
                omim_ids = id_map.get(nbk_num)
                for omim_num in omim_ids:
                    omim_id = 'OMIM:'+omim_num
                    # add the gene reviews as a superclass to the omim id,
                    # but only if the omim id is not a gene
                    if omim_id in entries_that_are_phenotypes:
                        model.addClassToGraph(omim_id, None)
                        model.addSubClass(omim_id, gr_id)
            # add this as a generic subclass of DOID:4
            model.addSubClass(gr_id, 'DOID:4')

        return

    def _get_titles(self, limit):
        """
        The file processed here is of the format:
        #NBK_id GR_shortname    OMIM
        NBK1103 trimethylaminuria       136132
        NBK1103 trimethylaminuria       602079
        NBK1104 cdls    122470
        Where each of the rows represents a mapping between
        a gr id and an omim id. These are a 1:many relationship,
        and some of the omim ids are genes (not diseases).
        Therefore, we need to create a loose coupling here.
        We make the assumption that these NBKs are generally higher-level
        grouping classes; therefore the OMIM ids are treated as subclasses.
        (This assumption is poor for those omims that are actually genes,
        but we have no way of knowing what those are here...
        we will just have to deal with that for now.)
        :param limit:
        :return:
        """
        raw = '/'.join((self.rawdir, self.files['titles']['file']))
        model = Model(self.graph)
        line_counter = 0
        with open(raw, 'r', encoding='latin-1') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                if line_counter == 1:  # skip header
                    continue
                (shortname, title, nbk_num) = row
                gr_id = 'GeneReviews:'+nbk_num

                self.book_ids.add(nbk_num)  # a global set of the book nums

                if limit is None or line_counter < limit:
                    model.addClassToGraph(gr_id, title)
                    model.addSynonym(gr_id, shortname)

        return

    def create_books(self):

        # note that although we put in the url to the book,
        # NCBI Bookshelf does not allow robots to download content
        book_item = {'file': 'books/',
                     'url': ''}

        for nbk in self.book_ids:
            b = book_item.copy()
            b['file'] = '/'.join(('books', nbk+'.html'))
            b['url'] = 'http://www.ncbi.nlm.nih.gov/books/'+nbk
            self.all_books[nbk] = b

        return

    def process_nbk_html(self, limit):
        """
        Here we process the gene reviews books to fetch
        the clinical descriptions to include in the ontology.
        We only use books that have been acquired manually,
        as NCBI Bookshelf does not permit automated downloads.
        This parser will only process the books that are found in
        the ```raw/genereviews/books``` directory,
        permitting partial completion.

        :param limit:
        :return:
        """
        model = Model(self.graph)
        c = 0
        books_not_found = set()
        for nbk in self.book_ids:
            c += 1
            nbk_id = 'GeneReviews:'+nbk
            book_item = self.all_books.get(nbk)
            url = '/'.join((self.rawdir, book_item['file']))

            # figure out if the book is there; if so, process, otherwise skip
            book_dir = '/'.join((self.rawdir, 'books'))
            book_files = os.listdir(book_dir)
            if ''.join((nbk, '.html')) not in book_files:
                # logger.warning("No book found locally for %s; skipping", nbk)
                books_not_found.add(nbk)
                continue
            logger.info("Processing %s", nbk)

            page = open(url)
            soup = BeautifulSoup(page.read())

            # sec0 == clinical description
            clin_summary = \
                soup.find(
                    'div', id=re.compile(".*Summary.sec0"))
            if clin_summary is not None:
                p = clin_summary.find('p')
                ptext = p.text
                ptext = re.sub(r'\s+', ' ', ptext)

                ul = clin_summary.find('ul')
                if ul is not None:
                    item_text = list()
                    for li in ul.find_all('li'):
                        item_text.append(re.sub(r'\s+', ' ', li.text))
                    ptext += ' '.join(item_text)

                # add in the copyright and citation info to description
                ptext = \
                    ' '.join(
                        (ptext,
                         '[GeneReviews:NBK1116, GeneReviews:NBK138602, ' +
                         nbk_id+']'))

                model.addDefinition(nbk_id, ptext.strip())

            # get the pubs
            pmid_set = set()
            pub_div = soup.find('div', id=re.compile(r".*Literature_Cited"))
            if pub_div is not None:
                ref_list = pub_div.find_all('div', attrs={'class': "bk_ref"})
                for r in ref_list:
                    for a in r.find_all(
                            'a', attrs={'href': re.compile(r"pubmed")}):
                        if re.match(r'PubMed:', a.text):
                            pmnum = re.sub(r'PubMed:\s*', '', a.text)
                        else:
                            pmnum = \
                                re.search(
                                    r'\/pubmed\/(\d+)$', a['href']).group(1)
                        if pmnum is not None:
                            pmid = 'PMID:'+str(pmnum)
                            self.graph.addTriple(
                                pmid,
                                model.object_properties['is_about'],
                                nbk_id)
                            pmid_set.add(pmnum)
                            reference = Reference(
                                self.graph,
                                pmid, Reference.ref_types['journal_article'])
                            reference.addRefToGraph()

            # TODO add author history, copyright, license to dataset

            # TODO get PMID-NBKID equivalence (near foot of page),
            # and make it "is about" link
            # self.gu.addTriple(
            #   self.graph, pmid,
            #   self.gu.object_properties['is_about'], nbk_id)
            # for example: NBK1191 PMID:20301370

            # add the book to the dataset
            self.dataset.setFileAccessUrl(book_item['url'])

            if limit is not None and c > limit:
                break

            # finish looping through books

        l = len(books_not_found)
        if len(books_not_found) > 0:
            if l > 100:
                logger.warning("There were %d books not found.", l)
            else:
                logger.warning(
                    "The following %d books were not found locally: %s",
                    l, str(books_not_found))
        logger.info(
            "Finished processing %d books for clinical descriptions", c-l)

        return

    def getTestSuite(self):
        import unittest
        from tests.test_genereviews import GeneReviewsTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(GeneReviewsTestCase)

        return test_suite
