import re
import os
import csv
import logging

from bs4 import BeautifulSoup
from dipper.sources.OMIMSource import OMIMSource
from dipper.models.Model import Model
from dipper.models.Reference import Reference

__author__ = 'nicole'

LOG = logging.getLogger(__name__)
GRDL = 'http://ftp.ncbi.nih.gov/pub/GeneReviews'


class GeneReviews(OMIMSource):
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
        'idmap': {
            'file': 'NBKid_shortname_OMIM.txt',
            'url': GRDL + '/NBKid_shortname_OMIM.txt'
        },
        'titles': {
            'file': 'GRtitle_shortname_NBKid.txt',
            'url': GRDL + '/GRtitle_shortname_NBKid.txt'
        },
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skolemized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='genereviews',
            ingest_title='Gene Reviews',
            ingest_url='http://genereviews.org/',
            ingest_logo='source-genereviews.png',
            license_url=None,
            data_rights='http://www.ncbi.nlm.nih.gov/books/NBK138602/',
            # file_handle=None
        )

        self.dataset.set_citation('GeneReviews:NBK1116')

        self.book_ids = set()
        self.all_books = {}

        if 'disease' not in self.all_test_ids:
            LOG.warning("not configured with disease test ids.")
            self.test_ids = list()
        else:
            # select ony those test ids that are omim's.
            self.test_ids = self.all_test_ids['disease']

    def fetch(self, is_dl_forced=False):
        """
        We fetch GeneReviews id-label map and id-omim mapping files from NCBI.
        :return: None
        """
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        """
        :return: None
        """
        if self.test_only:
            self.test_mode = True

        self._get_titles(limit)
        self._get_equivids(limit)

        self.create_books()
        self.process_nbk_html(limit)

        # no test subset for now; test == full graph
        self.testgraph = self.graph

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

        :param limit:

        """
        raw = '/'.join((self.rawdir, self.files['idmap']['file']))
        model = Model(self.graph)
        LOG.info('Looping over %s', raw)
        # we look some stuff up in OMIM, so initialize here
        # omim = OMIM(self.graph_type, self.are_bnodes_skized)
        id_map = {}
        allomimids = set()
        col = ['NBK_id', 'GR_shortname', 'OMIM']

        with open(raw, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row = next(reader)
            row[0] = row[0][1:]
            if not self.check_fileheader(col, row):
                pass

            for row in reader:

                nbk_num = row[col.index('NBK_id')]
                shortname = row[col.index('GR_shortname')]
                omim_num = row[col.index('OMIM')]
                gr_id = 'GeneReviews:' + nbk_num
                omim_id = 'OMIM:' + omim_num
                if not (
                        (self.test_mode and
                         len(self.test_ids) > 0 and
                         omim_id in self.test_ids) or not
                        self.test_mode):
                    continue

                # sometimes there's bad omim nums
                omim_num = omim_num.strip()
                if len(omim_num) != 6:
                    LOG.warning(
                        "OMIM number incorrectly formatted in row %i; skipping:\n%s",
                        reader.line_num, '\t'.join(row))
                    continue

                # build up a hashmap of the mappings; then process later
                if nbk_num not in id_map:
                    id_map[nbk_num] = set()
                id_map[nbk_num].add(omim_num)

                # add the class along with the shortname
                model.addClassToGraph(gr_id, None)
                model.addSynonym(gr_id, shortname)

                allomimids.add(omim_num)

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break

            # end looping through file

        # given all_omim_ids from GR,
        # we want to update any which are changed or removed
        # before deciding which are disease / phenotypes
        replaced = allomimids & self.omim_replaced.keys()
        if replaced is not None and len(replaced) > 0:
            LOG.warning("These OMIM ID's are past their pull date: %s", str(replaced))
            for oid in replaced:
                allomimids.remove(oid)
                replacements = self.omim_replaced[oid]
                for rep in replacements:
                    allomimids.update(rep)
        # guard against omim identifiers which have been removed
        obsolete = [
            o for o in self.omim_type
            if self.omim_type[o] == self.globaltt['obsolete']]
        removed = allomimids & set(obsolete)
        if removed is not None and len(removed) > 0:
            LOG.warning("These OMIM ID's are gone: %s", str(removed))
            for oid in removed:
                allomimids.remove(oid)
        # filter for disease /phenotype types (we can argue about what is included)
        omim_phenotypes = set([
            omim for omim in self.omim_type if self.omim_type[omim] in (
                self.globaltt['phenotype'],
                self.globaltt['has_affected_feature'],  # both a gene and a phenotype
                self.globaltt['heritable_phenotypic_marker'])])  # probable phenotype
        LOG.info(
            "Have %i omim_ids globally typed as phenotypes from OMIM",
            len(omim_phenotypes))

        entries_that_are_phenotypes = allomimids & omim_phenotypes
        LOG.info(
            "Filtered out %d/%d entries that are genes or features",
            len(allomimids - entries_that_are_phenotypes), len(allomimids))

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
            # add this as a generic subclass  -- TEC: this is the job of inference
            model.addSubClass(gr_id, self.globaltt['disease'])

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
        col = ['GR_shortname', 'GR_Title', 'NBK_id', 'PMID']
        with open(raw, 'r', encoding='latin-1') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row = next(reader)
            row[0] = row[0][1:]
            colcount = len(col)
            if not self.check_fileheader(col, row):
                pass
            for row in reader:
                if len(row) != colcount:
                    LOG.error("Unexpected row. got: %s", row)
                    LOG.error("Expected data for: %s", col)
                    exit(-1)
                nbk_num = row[col.index('NBK_id')]
                gr_id = 'GeneReviews:' + nbk_num
                self.book_ids.add(nbk_num)  # a global set of the book nums
                if limit is None or reader.line_num < limit:
                    model.addClassToGraph(gr_id, row[col.index('GR_Title')])
                    model.addSynonym(gr_id, row[col.index('GR_shortname')])
                # TODO include the new PMID?

    def create_books(self):

        # note that although we put in the url to the book,
        # NCBI Bookshelf does not allow robots to download content
        book_item = {
            'file': 'books/',
            'url': ''
        }

        for nbk in self.book_ids:
            nbki = book_item.copy()
            nbki['file'] = '/'.join(('books', nbk + '.html'))
            nbki['url'] = 'http://www.ncbi.nlm.nih.gov/books/' + nbk
            self.all_books[nbk] = nbki

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
        cnt = 0
        books_not_found = set()
        clin_des_regx = re.compile(r".*Summary.sec0")
        lit_cite_regex = re.compile(r".*Literature_Cited")
        pubmed_regex = re.compile(r"pubmed")  # ??? for a static string?
        for nbk in self.book_ids:
            cnt += 1
            nbk_id = 'GeneReviews:'+nbk
            book_item = self.all_books.get(nbk)
            url = '/'.join((self.rawdir, book_item['file']))

            # figure out if the book is there; if so, process, otherwise skip
            book_dir = '/'.join((self.rawdir, 'books'))
            book_files = os.listdir(book_dir)
            if ''.join((nbk, '.html')) not in book_files:
                # LOG.warning("No book found locally for %s; skipping", nbk)
                books_not_found.add(nbk)
                continue
            LOG.info("Processing %s", nbk)

            page = open(url)
            soup = BeautifulSoup(page.read())

            # sec0 == clinical description
            clin_summary = soup.find('div', id=clin_des_regx)
            if clin_summary is not None:
                ptext = clin_summary.find('p').text
                ptext = re.sub(r'\s+', ' ', ptext)

                unlst = clin_summary.find('ul')
                if unlst is not None:
                    item_text = list()
                    for lst_itm in unlst.find_all('li'):
                        item_text.append(re.sub(r'\s+', ' ', lst_itm.text))
                    ptext += ' '.join(item_text)

                # add in the copyright and citation info to description
                ptext = ' '.join((
                    ptext, '[GeneReviews:NBK1116, GeneReviews:NBK138602, ' +
                    nbk_id + ']'))

                model.addDefinition(nbk_id, ptext.strip())

            # get the pubs
            pmid_set = set()
            pub_div = soup.find('div', id=lit_cite_regex)
            if pub_div is not None:
                ref_list = pub_div.find_all('div', attrs={'class': "bk_ref"})
                for ref in ref_list:
                    for anchor in ref.find_all(
                            'a', attrs={'href': pubmed_regex}):
                        if re.match(r'PubMed:', anchor.text):
                            pmnum = re.sub(r'PubMed:\s*', '', anchor.text)
                        else:
                            pmnum = re.search(
                                r'\/pubmed\/(\d+)$', anchor['href']).group(1)
                        if pmnum is not None:
                            pmid = 'PMID:'+str(pmnum)
                            self.graph.addTriple(
                                pmid, self.globaltt['is_about'], nbk_id)
                            pmid_set.add(pmnum)
                            reference = Reference(
                                self.graph, pmid, self.globaltt['journal article'])
                            reference.addRefToGraph()

            # TODO add author history, copyright, license to dataset

            # TODO get PMID-NBKID equivalence (near foot of page),
            # and make it "is about" link
            # self.gu.addTriple(
            #   self.graph, pmid,
            #   self.globaltt['is_about'], nbk_id)
            # for example: NBK1191 PMID:20301370

            # add the book to the dataset
            self.dataset.set_ingest_source(book_item['url'])

            if limit is not None and cnt > limit:
                break

            # finish looping through books

        bknfd = len(books_not_found)
        if len(books_not_found) > 0:
            if bknfd > 100:
                LOG.warning("There were %d books not found.", bknfd)
            else:
                LOG.warning(
                    "The following %d books were not found locally: %s", bknfd,
                    str(books_not_found))
        LOG.info("Finished processing %d books for clinical descriptions", cnt - bknfd)

    def getTestSuite(self):
        import unittest
        from tests.test_genereviews import GeneReviewsTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(GeneReviewsTestCase)

        return test_suite
