import logging
import re
import csv
from dipper.sources.Source import Source
from dipper import config

LOG = logging.getLogger(__name__)


# omimftp key EXPIRES
# get a new one here: https://omim.org/help/api

OMIMURL = 'https://data.omim.org/downloads/'
OMIMFTP = OMIMURL + config.get_config()['keys']['omim']
USER_AGENT = "The Monarch Initiative (https://monarchinitiative.org/; " \
             "info@monarchinitiative.org)"


class OMIMSource(Source):
    '''
    Intended as a shim between Source.py and an ingest to make omim information
    available in a consistant way, without requireing other ingests to import config,py.

    populates two dicts
        one for omims which are replaced
        one for our interpertation of an omims type as an ontology term.

    These dicts are to be used to guide the treatment of omims found in ingests.

    note: If an omim is not a key in either table it is presumed removed.
    '''

    mimfiles = {  # do not conflict with subclasses 'files' dict
        'mimtitles': {
            'file': 'mimTitles.txt',
            'url':  OMIMFTP + '/mimTitles.txt',
            'clean': OMIMURL,
            'headers': {'User-Agent': USER_AGENT},
            'columns': [  # expected
                'Prefix',
                'Mim Number',
                'Preferred Title; symbol',
                'Alternative Title(s); symbol(s)',
                'Included Title(s); symbols',
            ]
        },
    }

    def __init__(
            self,
            graph_type,
            are_bnodes_skolemized,
            data_release_version=None,
            name=None,
            ingest_title=None,
            ingest_url=None,
            ingest_logo=None,
            license_url=None,
            data_rights=None,
            file_handle=None
    ):

        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name=name,
            ingest_title=ingest_title,
            ingest_url=ingest_url,
            ingest_logo=ingest_logo,
            license_url=license_url,
            data_rights=data_rights,
            file_handle=file_handle)

        self.omim_type = {}
        self.omim_replaced = {}
        self.populate_omim_type()

    # abstract
    def fetch(self, is_dl_forced=False):
        """
        abstract method to fetch all data from an external resource.
        this should be overridden by subclasses
        :return: None

        """
        raise NotImplementedError

    def parse(self, limit):
        """
        abstract method to parse all data from an external resource,
        that was fetched in fetch() this should be overridden by subclasses
        :return: None

        """
        raise NotImplementedError

    def populate_omim_type(self):
        '''
        This utility may still need to be rehomed. have considered:
            - make it its own ingest, output rdf
            - make it a dipper utility
            - inherited method from the super class
            - make a intermediate subclass between Source & ingest
        trying the last out as path of least [resistance, kludgy, questions]

        Use OMIM's discription of their identifiers
        to heuristically partition them into genes | phenotypes-diseases
        and shades between.

        populate omim_type map from an omim number to an ontology term
        the ontology terms's labels as:

        -  'gene'
            Asterisk (*)  Gene

        -   'has_affected_feature'
            Plus (+)  Gene and phenotype, combined

        -   'Phenotype'
            Number Sign (#)  Phenotype, molecular basis known

        -   'heritable_phenotypic_marker'
            Percent (%)  Phenotype or locus, molecular basis unknown

        -   'obsolete'
            Caret (^)  Entry has been removed from the database
            or moved to another entry.
                further processed into Removed or Moved & Split
                (`omim_replaced`  populated where Moved or Split )

        -   'Suspected'
            NULL (<null>)  Other, mainly phenotypes with suspected mendelian basis

        Populates dict of omim_number to ontology_curie
        Populates dict of omim_number to list to replacments

        note:
            If an omim id is neither in omim_replaced nor omim_types
            Then it was removed.

        '''

        src_key = 'mimtitles'

        self.get_files(is_dl_forced=False, files=self.mimfiles)

        myfile = '/'.join((self.rawdir, self.mimfiles[src_key]['file']))

        col = self.mimfiles[src_key]['columns']
        with open(myfile, 'r') as readfile:
            reader = csv.reader(readfile, delimiter='\t')
            row = next(reader)  # copyright
            row = next(reader)  # date  generated
            row = next(reader)  # column header
            row[0] = row[0][2:]  # remove octothorp '# '
            if not self.check_fileheader(col, row):
                pass

            for row in reader:
                if row[0][0] == '#':  # they have comments at the end
                    continue

                declared = row[col.index('Prefix')].strip()
                omim_id = row[col.index('Mim Number')].strip()
                pref_label = row[col.index('Preferred Title; symbol')].strip()
                # alt_label =row[col.index('Alternative Title(s); symbol(s)')].strip()
                # inc_label = row[col.index('Included Title(s); symbols')].strip()

                if declared == 'Caret':  # moved|removed|split -> moved twice
                    # populating a dict from an omim to a set of omims
                    self.omim_type[omim_id] = self.globaltt['obsolete']
                    self.omim_replaced[omim_id] = []
                    if pref_label[:9] == 'MOVED TO ':
                        token = pref_label.split(' ')
                        rep = token[2]
                        if not re.match(r'^[0-9]{6}$', rep):
                            LOG.error('Report malformed omim replacement %s', rep)
                            # clean up ones I know about
                            if rep[0] == '{' and rep[7] == '}':
                                rep = rep[1:7]
                                LOG.info('Repaired malformed omim replacement %s', rep)
                            if len(rep) == 7 and rep[6] == ',':
                                rep = rep[:6]
                                LOG.info('Repaired malformed omim replacement %s', rep)
                        if len(token) > 3:
                            self.omim_replaced[omim_id] = [rep, token[4]]
                        else:
                            self.omim_replaced[omim_id] = [rep]

                elif declared == 'Asterisk':
                    self.omim_type[omim_id] = self.globaltt['gene']
                elif declared == 'NULL':
                    self.omim_type[omim_id] = self.globaltt['Suspected']   # NCIT:C71458
                elif declared == 'Number Sign':
                    self.omim_type[omim_id] = self.globaltt['phenotype']
                elif declared == 'Percent':
                    self.omim_type[omim_id] = self.globaltt[
                        'heritable_phenotypic_marker']
                elif declared == 'Plus':
                    self.omim_type[omim_id] = self.globaltt['has_affected_feature']
                else:
                    LOG.error('Unknown OMIM type line %s', reader.line_num)
        LOG.info('Have %i obsoleted OMIMS ids', len(self.omim_replaced))
        LOG.info('Have %i OMIM typed', len(self.omim_type))
'''
cut -f2 raw/omim/mim2gene.txt | grep -v '^#' | dist
  16068 gene
   7095 phenotype                   ->      Number Sign AND   Percent
   1754 predominantly phenotypes
   1268 moved/removed
     45 gene/phenotype

cut -f1 raw/omim/mimTitles.txt | grep -v  '^#' | dist
  16068 Asterisk        ->  mim2gene.gene
   5531 Number Sign     ->  mim2gene.phenotype*
   1754 NULL            ->  mim2gene.predominantly phenotypes
   1564 Percent         ->  mim2gene.phenotype*
   1268 Caret           ->  mim2gene.moved/removed
     45 Plus            ->  mim2gene.gene/phenotype


    gene-ish        (Asterisk, Plus)
    phenotype-ish   ( Number Sign, NULL, Percent, Plus )

26230 prod calls
24962 dev calls   (no obs)

'''
