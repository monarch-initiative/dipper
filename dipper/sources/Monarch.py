import csv
import re
import logging
from os import listdir
from os.path import isfile, join

from dipper.sources.Source import Source
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class Monarch(Source):
    """
    This is the parser for data curated by the
    [Monarch Initiative](https://monarchinitiative.org).
    Data is currently maintained in a private repository, soon to be released.

    """

    files = {
        'omia_d2p': {
            'file': r'[0-9]{6}.txt',
            'columns': [
                'Disease ID',
                'Species ID',
                'Breed Name',
                'Variant',
                'Inheritance',
                'Phenotype ID',
                'Phenotype Name',
                'Entity ID',
                'Entity Name',
                'Quality ID',
                'Quality Name',
                'Related Entity ID',
                'Related Entity Name',
                'Abnormal ID',
                'Abnormal Name',
                'Phenotype Desc',
                'Assay',
                'Frequency',
                'Pubmed ID',
                'Pub Desc',
                'Curator Notes',
                'Date Created',
            ]
        }
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='monarch',
            ingest_title='The Monarch Initiative',
            ingest_url='https://monarchinitiative.org',
            ingest_logo='source-monarch.png',
            license_url='https://creativecommons.org/licenses/by/4.0/'
            # data_rights=None,
            # file_handle=None
        )
        return

    def fetch(self, is_dl_forced=False):

        # fetch all the files
        # self.get_files(is_dl_forced)

        LOG.info(
            "Temporarily using local files until they move to public git")

        return

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)
        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self.process_omia_phenotypes(limit)
        LOG.info("Finished parsing.")

        return

    def process_omia_phenotypes(self, limit):

        # process the whole directory
        # TODO get the file listing
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        model = Model(graph)

        LOG.info("Processing Monarch OMIA Animal disease-phenotype associations")

        src_key = 'omia_d2p'

        # get file listing
        mypath = '/'.join((self.rawdir, 'OMIA-disease-phenotype'))
        file_list = [
            f for f in listdir(mypath)
            if isfile(join(mypath, f)) and re.search(r'.txt$', f)]

        col = self.files[src_key]['columns']
        # reusable initial code generator
        # for c in col:
        #   print(
        #    '# '+str.lower(c.replace(" ",""))+" = row[col.index('"+c+"')].strip()")

        for filename in file_list:
            LOG.info("Processing %s", filename)
            count_missing = 0
            bad_rows = list()
            fname = '/'.join((mypath, filename))
            with open(fname, 'r') as csvfile:
                filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
                fileheader = next(filereader)
                if fileheader != col:
                    LOG.error('Expected  %s to have columns: %s', fname, col)
                    LOG.error('But Found %s to have columns: %s', fname, fileheader)
                    raise AssertionError('Incoming data headers have changed.')

                for row in filereader:
                    if len(row) != len(col):
                        LOG.info(
                            "Not enough cols %d in %s - please fix", len(row), filename)
                        continue

                    disease_num = row[col.index('Disease ID')].strip()
                    species_id = row[col.index('Species ID')].strip()
                    breed_name = row[col.index('Breed Name')].strip()
                    # variant = row[col.index('Variant')]
                    # inheritance = row[col.index('Inheritance')]
                    phenotype_id = row[col.index('Phenotype ID')].strip()
                    # phenotype_name = row[col.index('Phenotype Name')]
                    entity_id = row[col.index('Entity ID')].strip()
                    entity_name = row[col.index('Entity Name')]
                    quality_id = row[col.index('Quality ID')].strip()
                    quality_name = row[col.index('Quality Name')]
                    # related_entity_id = row[col.index('Related Entity ID')]
                    # related_entity_name = row[col.index('Related Entity Name')]
                    # abnormal_id = row[col.index('Abnormal ID')]
                    # abnormal_name = row[col.index('Abnormal Name')]
                    # phenotype_desc = row[col.index('Phenotype Desc')]
                    assay = row[col.index('Assay')].strip()
                    # frequency = row[col.index('Frequency')]
                    pubmed_id = row[col.index('Pubmed ID')].strip()
                    phenotype_description = row[col.index('Pub Desc')].strip()
                    curator_notes = row[col.index('Curator Notes')].strip()
                    # date_created = row[col.index('Date Created')]

                    if phenotype_id == '':
                        # LOG.warning('Missing phenotype in row:\n%s', row)
                        count_missing += 1
                        bad_rows.append(row)
                        continue
                    if len(str(disease_num)) < 6:
                        disease_num = str(disease_num).zfill(6)
                    disease_id = 'OMIA:' + disease_num
                    if species_id != '':
                        disease_id = '-'.join((disease_id, species_id))
                    assoc = D2PAssoc(graph, self.name, disease_id, phenotype_id)
                    if pubmed_id != '':
                        for pnum in re.split(r'[,;]', pubmed_id):
                            pnum = re.sub(r'[^0-9]', '', pnum)
                            pmid = 'PMID:' + pnum
                            assoc.add_source(pmid)
                    else:
                        assoc.add_source(
                            '/'.join((
                                self.curie_map['OMIA'] + disease_num, species_id)))
                    assoc.add_association_to_graph()
                    aid = assoc.get_association_id()
                    if phenotype_description != '':
                        model.addDescription(aid, phenotype_description,
                                             subject_category=blv.terms.Association.value)
                    if breed_name != '':
                        model.addDescription(aid, breed_name + ' [observed in]',
                                             subject_category=blv.terms.Association.value)
                    if assay != '':
                        model.addDescription(aid, assay + ' [assay]',
                                             subject_category=blv.terms.Association.value)
                    if curator_notes != '':
                        model.addComment(aid, curator_notes,
                                         subject_category=blv.terms.Association.value)

                    if entity_id != '' or quality_id != '':
                        LOG.info(
                            "EQ not empty for %s: %s + %s",
                            disease_id, entity_name, quality_name)
            if count_missing > 0:
                LOG.warning(
                    "We are missing %d of %d D2P annotations from id %s",
                    count_missing, filereader.line_num-1, filename)
                LOG.warning("Bad rows:\n%s", '\n'.join([str(x) for x in bad_rows]))
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
