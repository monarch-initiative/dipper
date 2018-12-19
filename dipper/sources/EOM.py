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
        'map': {
            'file': 'hp-to-eom-mapping.tsv',
            'url': GHRAW + '/master/src/mappings/hp-to-eom-mapping.tsv'
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'eom',
            ingest_title='Elements of Morphology',
            ingest_url='http://elementsofmorphology.nih.gov',
            data_rights='http://www.genome.gov/copyright.cfm',
            license_url='https://creativecommons.org/publicdomain/mark/1.0/'
            # file_handle=None
        )

        return

    def fetch(self, is_dl_forced=False):
        '''connection details for DISCO'''
        cxn = {}
        cxn['host'] = 'nif-db.crbs.ucsd.edu'
        cxn['database'] = 'disco_crawler'
        cxn['port'] = '5432'
        cxn['user'] = config.get_config()['user']['disco']
        cxn['password'] = config.get_config()['keys'][cxn['user']]

        self.dataset.setFileAccessUrl(
            'jdbc:postgresql://'+cxn['host']+':'+cxn['port']+'/'+cxn['database'],
            is_object_literal=True)

        # process the tables
        # self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables, cxn)

        self.get_files(is_dl_forced)

        # FIXME: Everything needed for data provenance?
        fstat = os.stat('/'.join((self.rawdir, 'dvp.pr_nlx_157874_1')))
        filedate = datetime.utcfromtimestamp(fstat[ST_CTIME]).strftime("%Y-%m-%d")
        self.dataset.setVersion(filedate)

        return

    def parse(self, limit=None):
        '''
            Over ride Source.parse inherited via PostgreSQLSource
        '''

        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)

        if self.testOnly:
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

        return

    def _process_nlx_157874_1_view(self, raw, limit=None):
        """
        This table contains the Elements of Morphology data that has been
        screen-scraped into DISCO.
        Note that foaf:depiction is inverse of foaf:depicts relationship.

        Since it is bad form to have two definitions,
        we concatenate the two into one string.

        Turtle:
            <eom id> a owl:Class
                rdf:label Literal(eom label)
                OIO:hasRelatedSynonym Literal(synonym list)
                IAO:definition Literal(objective_def. subjective def)
                foaf:depiction Literal(small_image_url),
                               Literal(large_image_url)
                foaf:page Literal(page_url)
                rdfs:comment Literal(long commented text)


        :param raw:
        :param limit:
        :return:
        """

        model = Model(self.graph)
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            reader = csv.reader(f1, delimiter='\t', quotechar='\"')
            for line in reader:
                (morphology_term_id, morphology_term_num,
                 morphology_term_label, morphology_term_url,
                 terminology_category_label, terminology_category_url,
                 subcategory, objective_definition, subjective_definition,
                 comments, synonyms, replaces, small_figure_url,
                 large_figure_url, e_uid, v_uid, v_uuid,
                 v_last_modified, v_status, v_lastmodified_epoch) = line

                # note:
                # e_uid v_uuid v_last_modified terminology_category_url
                # subcategory v_uid morphology_term_num
                # terminology_category_label hp_label notes
                # are currently unused.

                # Add morphology term to graph as a class
                # with label, type, and description.
                model.addClassToGraph(morphology_term_id, morphology_term_label)

                # Assemble the description text

                if subjective_definition != '' and not (
                        re.match(r'.+\.$', subjective_definition)):
                    # add a trailing period.
                    subjective_definition = subjective_definition.strip() + '.'
                if objective_definition != '' and not (
                        re.match(r'.+\.$', objective_definition)):
                    # add a trailing period.
                    objective_definition = objective_definition.strip() + '.'

                definition = '  '.join(
                    (objective_definition, subjective_definition)).strip()

                model.addDefinition(morphology_term_id, definition)

                # <term id> FOAF:depicted_by literal url
                # <url> type foaf:depiction

                # do we want both images?
                # morphology_term_id has depiction small_figure_url
                if small_figure_url != '':
                    model.addDepiction(morphology_term_id,
                                       small_figure_url)

                # morphology_term_id has depiction large_figure_url
                if large_figure_url != '':
                    model.addDepiction(morphology_term_id,
                                       large_figure_url)

                # morphology_term_id has comment comments
                if comments != '':
                    model.addComment(morphology_term_id, comments.strip())

                for syn in synonyms.split(';'):
                    model.addSynonym(
                        morphology_term_id, syn.strip(),
                        self.globaltt['hasExactSynonym'])

                # morphology_term_id hasRelatedSynonym replaces (; delimited)
                if replaces != '' and replaces != synonyms:
                    for syn in replaces.split(';'):
                        model.addSynonym(
                            morphology_term_id, syn.strip(),
                            self.globaltt['hasRelatedSynonym'])

                # <morphology_term_id> <foaf:page> morphology_term_url
                if morphology_term_id is not None and morphology_term_url is not None:
                    reference = Reference(
                        self.graph, morphology_term_url, self.globaltt['web page'])
                    reference.addPage(morphology_term_id, morphology_term_url)

                if limit is not None and reader.line_num > limit:
                    break
        return

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
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1
                row = line.split('\t')
                (
                    morphology_term_id, morphology_term_label, hp_id, hp_label,
                    notes) = row

                # Sub out the underscores for colons.
                hp_id = re.sub('_', ':', hp_id)
                if re.match(".*HP:.*", hp_id):
                    # add the HP term as a class
                    model.addClassToGraph(hp_id, None)
                    # Add the HP ID as an equivalent class
                    model.addEquivalentClass(morphology_term_id, hp_id)
                else:
                    LOG.warning('No matching HP term for %s', morphology_term_label)

                if limit is not None and line_counter > limit:
                    break

        return

    def getTestSuite(self):
        import unittest
        from tests.test_eom import EOMTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(EOMTestCase)

        return test_suite
