
import re
import csv
import logging

from dipper.sources.Source import Source, USER_AGENT
from dipper import config

LOG = logging.getLogger(__name__)

class OmimType(Source):
    '''
    '''

    OMIMURL = 'https://data.omim.org/downloads/'
    OMIMFTP = OMIMURL + config.get_config()['keys']['omim']

    files = {
        'mimtitles': {
            'file': 'mimTitles.txt',
            'url':  OMIMFTP + '/mimTitles.txt',
            'headers': {'User-Agent': USER_AGENT},
            'clean': OMIMURL,
            'columns': (  # expected
                'Prefix',
                'Mim Number',
                'Preferred Title; symbol',
                'Alternative Title(s); symbol(s)',
                'Included Title(s); symbols',
            ),
        },
    }

##################################################

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'omimtype',
            ingest_title='OMIM Types',
            ingest_url='https://data.omim.org',
            license_url=None,
            data_rights='http://omim.org/help/agreement',
            # file_handle=None
        )

        self.omim_replaced = {}  # id_num to SET of id nums
        self.omim_type = {}      # id_num to onto_term

        return

##################################################

    def fetch(self, is_dl_forced=False):
        """

        """
        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        """
        :return: None
        """
        self.find_omim_type()

        return

    def write(self):

        for num in omim_type:
            stmt = " ".join('OMIM:'+num, self.globaltt['domain'], omim_type[num],)+'.'

        return

    def find_omim_type(self):
        '''
        Use OMIM's discription of their identifiers
        to heuristically partition them into genes | phenotypes-diseases
        type could be
            - `obsolete`  Check `omim_replaced`  populated as side effect
            - 'Suspected' (phenotype)  Ignoring thus far
            - 'gene'
            - 'Phenotype'
            - 'heritable_phenotypic_marker'   Probable phenotype
            - 'has_affected_feature'  Use as both a gene and a phenotype

        :return hash of omim_number to ontology_curie
        '''
        myfile = '/'.join((self.rawdir, self.files['mimtitles']['file']))
        omim_type = {}
        line_counter = 1
        with open(myfile, 'r') as fh:
            reader = csv.reader(fh, delimiter='\t')
            for row in reader:
                line_counter += 1
                if row[0][0] == '#':     # skip comments
                    continue
                elif row[0] == 'Caret':  # moved|removed|split -> moved twice
                    # populating a dict from an omim to a set of omims
                    # here as a side effect which is less than ideal
                    (prefix, omim_id, destination, empty, empty) = row
                    omim_type[omim_id] = self.globaltt['obsolete']
                    if row[2][:9] == 'MOVED TO ':
                        token = row[2].split(' ')
                        rep = token[2]
                        if not re.match(r'^[0-9]{6}$', rep):
                            LOG.error('Report malformed omim replacement %s', rep)
                            # clean up ones I know about
                            if rep[0] == '{' and rep[7] == '}':
                                rep = rep[1:6]
                            if len(rep) == 7 and rep[6] == ',':
                                rep = rep[:5]
                        # asuming splits are typically to both gene & phenotype
                        if len(token) > 3:
                            self.omim_replaced[omim_id] = {rep, token[4]}
                        else:
                            self.omim_replaced[omim_id] = {rep}

                elif row[0] == 'Asterisk':  # declared as gene
                    (prefix, omim_id, pref_label, alt_label, inc_label) = row
                    omim_type[omim_id] = self.globaltt['gene']
                elif row[0] == 'NULL':
                    #  potential model of disease?
                    (prefix, omim_id, pref_label, alt_label, inc_label) = row
                    #
                    omim_type[omim_id] = self.globaltt['Suspected']   # NCIT:C71458
                elif row[0] == 'Number Sign':
                    (prefix, omim_id, pref_label, alt_label, inc_label) = row
                    omim_type[omim_id] = self.globaltt['Phenotype']
                elif row[0] == 'Percent':
                    (prefix, omim_id, pref_label, alt_label, inc_label) = row
                    omim_type[omim_id] = self.globaltt['heritable_phenotypic_marker']
                elif row[0] == 'Plus':
                    (prefix, omim_id, pref_label, alt_label, inc_label) = row
                    # to be interperted as  a gene and/or a phenotype
                    omim_type[omim_id] = self.globaltt['has_affected_feature']
                else:
                    LOG.error('Unknown OMIM type line %i :', line_counter)
        return omim_type
