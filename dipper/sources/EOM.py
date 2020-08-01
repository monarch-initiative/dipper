import os
from datetime import datetime
from stat import ST_CTIME
import re
import logging
import csv

from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.Model import Model
from dipper import config
from dipper.models.Reference import Reference
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class EOM(PostgreSQLSource):
    """
    Elements of Morphology is a resource from NHGRI that has definitions of
    morphological abnormalities, together with image depictions.
    We pull those relationships, as well as our local mapping of equivalences
    between EOM and HP terminologies.

    The website is crawled monthly by NIF's DISCO crawler system,
    which we utilize here.
    Be sure to have pg user/password connection details in your conf.yaml file,
    like:
    dbauth : {'disco' : {'user' : '<username>', 'password' : '<password>'}}

    Monarch-curated data for the HP to EOM mapping is stored at
    https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/src/mappings/hp-to-eom-mapping.tsv

    Since this resource is so small, the entirety of it is the "test" set.

    """

    # we are using the production view here; should we be using services?
    tables = [
        'dvp.pr_nlx_157874_1'
    ]

    GHRAW = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology'
    files = {
        'map': {  # in 2019 this repo has beem untouched in four years
            'file': 'hp-to-eom-mapping.tsv',
            'url': GHRAW + '/master/src/mappings/hp-to-eom-mapping.tsv',
            'columns': [
                'morphology_term_id',
                'morphology_term_label',
                'HP ID',
                'HP Label',
                'Notes'
            ]
        }
    }
    resources = {
        'tables': {
            'url': 'nif-db.crbs.ucsd.edu:5432',
            'file': 'dvp.pr_nlx_157874_1',
            'columns': [
                # head -1 dvp.pr_nlx_157874_1 | tr '\t' '\n' |sed "s|\(.*\)|'\1',|g"
                'morphology_term_id',
                'morphology_term_num',
                'morphology_term_label',
                'morphology_term_url',
                'terminology_category_label',
                'terminology_category_url',
                'subcategory',
                'objective_definition',
                'subjective_definition',
                'comments',
                'synonyms',
                'replaces',
                'small_figure_url',
                'large_figure_url',
                'e_uid',
                'v_uid',
                'v_uuid',
                'v_lastmodified',
                'v_status',
                'v_lastmodified_epoch'
            ],
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
            name='eom',
            ingest_title='Elements of Morphology',
            ingest_url='http://elementsofmorphology.nih.gov',
            ingest_logo='source-eom.png',
            data_rights='http://www.genome.gov/copyright.cfm',
            license_url='https://creativecommons.org/publicdomain/mark/1.0/'
            # file_handle=None
        )

    def fetch(self, is_dl_forced=False):
        '''connection details for DISCO'''
        cxn = {}
        cxn['host'] = 'nif-db.crbs.ucsd.edu'
        cxn['database'] = 'disco_crawler'
        cxn['port'] = '5432'
        cxn['user'] = config.get_config()['user']['disco']
        cxn['password'] = config.get_config()['keys'][cxn['user']]

        pg_iri = 'jdbc:postgresql://'+cxn['host']+':'+cxn['port']+'/'+cxn['database']
        self.dataset.set_ingest_source(pg_iri)

        # process the tables
        # self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables, cxn)

        self.get_files(is_dl_forced)

        fstat = os.stat('/'.join((self.rawdir, 'dvp.pr_nlx_157874_1')))
        filedate = datetime.utcfromtimestamp(fstat[ST_CTIME]).strftime("%Y-%m-%d")

        self.dataset.set_ingest_source_file_version_date(pg_iri, filedate)

    def parse(self, limit=None):
        '''
            Over ride Source.parse inherited via PostgreSQLSource
        '''

        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)

        if self.test_only:
            self.test_mode = True

        LOG.info("Parsing files...")

        self._process_nlx_157874_1_view(
            '/'.join((self.rawdir, 'dvp.pr_nlx_157874_1')), limit)
        self._map_eom_terms(
            '/'.join((self.rawdir, self.files['map']['file'])), limit)

        LOG.info("Finished parsing.")

        # since it's so small,
        # we default to copying the entire graph to the test set
        self.testgraph = self.graph

    def _process_nlx_157874_1_view(self, raw, limit=None):
        """
        This table contains the Elements of Morphology data .
        Note that foaf:depiction is inverse of foaf:depicts relationship.

        Since it is bad form to have two definitions,
        we concatenate the two into one string.

        Turtle:
            <eom id> a owl:Class
                rdfs:label Literal(eom label)
                oboInOwl:has_related_synonym Literal(synonym list)
                IAO:definition Literal(objective_def. subjective def)
                foaf:depiction Literal(small_image_url),
                               Literal(large_image_url)
                foaf:page Literal(page_url)
                rdfs:comment Literal(long commented text)

        TEC_note: URL are not literals.


        :param raw:
        :param limit:
        :return:
        """

        src_key = 'tables'
        model = Model(self.graph)
        col = self.resources[src_key]['columns']
        with open(raw, 'r') as rawread:
            reader = csv.reader(rawread, delimiter='\t', quotechar='\"')
            row = next(reader)
            if not self.check_fileheader(col, row):
                pass

            for row in reader:
                # head -1 dvp.pr_nlx_157874_1|tr '\t' '\n'|
                # sed "s|\(.*\)|# \1 = row[col.index('\1')]|g"

                morphology_term_id = row[col.index('morphology_term_id')].strip()
                # morphology_term_num = row[col.index('morphology_term_num')]
                morphology_term_label = row[col.index('morphology_term_label')].strip()
                morphology_term_url = row[col.index('morphology_term_url')].strip()
                # terminology_category_label = row[
                #   col.index('terminology_category_label')]
                # terminology_category_url = row[col.index('terminology_category_url')]
                # subcategory = row[col.index('subcategory')]
                objective_definition = row[col.index('objective_definition')].strip()
                subjective_definition = row[col.index('subjective_definition')].strip()
                comments = row[col.index('comments')].strip()
                synonyms = row[col.index('synonyms')].strip()
                replaces = row[col.index('replaces')].strip()
                small_figure_url = row[col.index('small_figure_url')].strip()
                large_figure_url = row[col.index('large_figure_url')].strip()
                # e_uid = row[col.index('e_uid')]
                # v_uid = row[col.index('v_uid')]
                # v_uuid = row[col.index('v_uuid')]
                # v_lastmodified = row[col.index('v_lastmodified')]
                # v_status = row[col.index('v_status')]
                # v_lastmodified_epoch = row[col.index('v_lastmodified_epoch')]

                # Add morphology term to graph as a class
                # with label, type, and description.
                model.addClassToGraph(
                    morphology_term_id,
                    morphology_term_label,
                    blv.terms['PhenotypicFeature']
                )

                # Assemble the description text

                if subjective_definition != '' and not (
                        re.match(r'.+\.$', subjective_definition)):
                    # add a trailing period.
                    subjective_definition = subjective_definition + '.'
                if objective_definition != '' and not (
                        re.match(r'.+\.$', objective_definition)):
                    # add a trailing period.
                    objective_definition = objective_definition + '.'

                definition = '  '.join(
                    (objective_definition, subjective_definition))

                model.addDefinition(morphology_term_id, definition,
                                    class_category=blv.terms['PhenotypicFeature'])

                # <term id> FOAF:depicted_by literal url
                # <url> type foaf:depiction

                # do we want both images?
                # morphology_term_id has depiction small_figure_url
                if small_figure_url != '':
                    model.addDepiction(morphology_term_id, small_figure_url)

                # morphology_term_id has depiction large_figure_url
                if large_figure_url != '':
                    model.addDepiction(morphology_term_id, large_figure_url)

                # morphology_term_id has comment comments
                if comments != '':
                    model.addComment(morphology_term_id, comments)

                for syn in synonyms.split(';'):
                    model.addSynonym(
                        morphology_term_id,
                        syn.strip(),
                        self.globaltt['has_exact_synonym']
                    )

                # morphology_term_id has_related_synonym replaces (; delimited)
                if replaces not in ['', synonyms]:
                    for syn in replaces.split(';'):
                        model.addSynonym(
                            morphology_term_id,
                            syn.strip(),
                            self.globaltt['has_related_synonym']
                        )

                # <morphology_term_id> <foaf:page> morphology_term_url
                if morphology_term_id is not None:
                    reference = Reference(
                        self.graph, morphology_term_id, self.globaltt['web page'])

                    # TEC 201905:
                    # Not so sure we need explicit   <eom_uri> <webpage> <eom_url>.
                    # since <eom_uri> IS the <eom_url>.

                    reference.addPage(morphology_term_id, morphology_term_url)

                if limit is not None and reader.line_num > limit:
                    break

    def _map_eom_terms(self, raw, limit=None):
        """
        This table contains the HP ID mappings from the local tsv file.
        Triples:
            <eom id> owl:equivalentClass <hp id>
        :param raw:
        :param limit:
        :return:
        """

        model = Model(self.graph)
        line_counter = 0
        col = self.files['map']['columns']
        with open(raw, 'r') as reader:
            line = reader.readline().strip()
            line = line.strip('/n')
            if self.check_fileheader(col, line.split('\t')):
                pass
            for line in reader:
                line_counter += 1
                row = line.strip('\n').split('\t')

                morphology_term_id = row[col.index('morphology_term_id')].strip()
                # morphology_term_label = row[col.index('morphology_term_label')]
                hp_id = row[col.index('HP ID')].strip()
                # hp_label = row[col.index('HP_Label')]
                # notes = row[col.index('Notes')]

                # Sub out the underscores for colons.
                hp_id = re.sub('_', ':', hp_id)
                if re.match(".*HP:.*", hp_id):
                    # add the HP term as a class
                    model.addClassToGraph(hp_id, None)
                    # Add the HP ID as an equivalent class
                    model.addEquivalentClass(morphology_term_id, hp_id)
                else:
                    LOG.warning('No matching HP term for %s', morphology_term_id)

                if limit is not None and line_counter > limit:
                    break

    def getTestSuite(self):
        import unittest
        from tests.test_eom import EOMTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(EOMTestCase)

        return test_suite
