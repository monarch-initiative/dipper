import csv
import re
import logging
from os import listdir
from os.path import isfile, join

from dipper.sources.Source import Source
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.Dataset import Dataset
from dipper.models.Model import Model

logger = logging.getLogger(__name__)


class Monarch(Source):
    """
    This is the parser for data curated by the
    [Monarch Initiative](https://monarchinitiative.org).
    Data is currently maintained in a private repository, soon to be released.

    """

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'monarch')

        self.dataset = Dataset(
            'monarch', 'MonarchInitiative', 'https://monarchinitiative.org',
            None, 'https://creativecommons.org/licenses/by/4.0/', None)

        return

    def fetch(self, is_dl_forced=False):

        # fetch all the files
        # self.get_files(is_dl_forced)

        logger.info(
            "Temporarily using local files until they move to public git")

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows of each file", limit)
        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        self.process_omia_phenotypes(limit)
        logger.info("Finished parsing.")

        return

    def process_omia_phenotypes(self, limit):

        # process the whole directory
        # TODO get the file listing
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        model = Model(g)

        logger.info(
            "Processing Monarch OMIA Animal disease-phenotype associations")

        # get file listing
        mypath = '/'.join((self.rawdir, 'OMIA-disease-phenotype'))
        file_list = [
            f for f in listdir(mypath)
            if isfile(join(mypath, f)) and re.search(r'.txt$', f)]

        for f in file_list:
            logger.info("Processing %s", f)
            print(f)
            line_counter = 0
            count_missing = 0
            bad_rows = list()
            fname = '/'.join((mypath, f))
            with open(fname, 'r') as csvfile:
                filereader = csv.reader(
                    csvfile, delimiter='\t', quotechar='\"')
                for row in filereader:
                    line_counter += 1
                    if line_counter <= 1:
                        continue  # skip header
                    if len(row) != 22:
                        logger.info("Not enough cols (%d) in %s - please fix",
                                    len(row), f)
                        continue
                    (disease_num, species_id, breed_name, variant, inheritance,
                     phenotype_id, phenotype_name, entity_id, entity_name,
                     quality_id, quality_name, related_entity_id,
                     related_entity_name, abnormal_id, abnormal_name,
                     phenotype_description, assay, frequency, pubmed_id,
                     pub_description, curator_notes, date_created) = row

                    if phenotype_id == '':
                        # logger.warning('Missing phenotype in row:\n%s', row)
                        count_missing += 1
                        bad_rows.append(row)
                        continue
                    if len(str(disease_num)) < 6:
                        disease_num = str(disease_num).zfill(6)
                    disease_id = 'OMIA:'+disease_num.strip()
                    species_id = species_id.strip()
                    if species_id != '':
                        disease_id = '-'.join((disease_id, species_id))
                    assoc = D2PAssoc(g, self.name, disease_id, phenotype_id)
                    if pubmed_id != '':
                        for p in re.split(r'[,;]', pubmed_id):
                            pmid = 'PMID:'+p.strip()
                            assoc.add_source(pmid)
                    else:
                        assoc.add_source(
                            '/'.join(('http://omia.angis.org.au/OMIA' +
                                      disease_num.strip(),
                                      species_id.strip())))
                    assoc.add_association_to_graph()
                    aid = assoc.get_association_id()
                    if phenotype_description != '':
                        model.addDescription(aid, phenotype_description)
                    if breed_name != '':
                        model.addDescription(
                            aid, breed_name.strip()+' [observed in]')
                    if assay != '':
                        model.addDescription(aid, assay.strip()+' [assay]')
                    if curator_notes != '':
                        model.addComment(aid, curator_notes.strip())

                    if entity_id != '' or quality_id != '':
                        logger.info("EQ not empty for %s: %s + %s", disease_id,
                                    entity_name, quality_name)
            if count_missing > 0:
                logger.warning(
                    "You are missing %d/%d D2P annotations from id %s",
                    count_missing, line_counter-1, f)
                # TODO PYLINT Used builtin function 'map'.
                # Using a list comprehension can be clearer.
                logger.warning("Bad rows:\n"+"\n".join(map(str, bad_rows)))
            # finish loop through all files

        return

    # def get_files_from_git(self):
    #
    #     base_url = 'https://raw.githubusercontent.com'
    #     path = \
    #       'monarch-initiative/data-boutique/master/OMIA-disease-phenotype/000060.txt?'
    #     params = {
    #         'token': None
    #     }
    #
    #     from git import Repo
    #     import os
    #     join = os.path.join
    #     git_url = 'https://github.com/'
    #     repo_dir = 'monarch-initiative/data-boutique/'
    #     repo = Repo.clone_from(git_url, repo_dir)
    #
    #
    #
    #     return
