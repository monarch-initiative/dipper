import csv
import os
from datetime import datetime
from stat import ST_CTIME
import logging
import re
import shutil
from git import Repo
from git import GitCommandError

from dipper.utils import pysed
from dipper.sources.Source import Source
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.assoc.DispositionAssoc import DispositionAssoc
from dipper.models.Dataset import Dataset
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper import config

logger = logging.getLogger(__name__)

HPOADL = 'http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc'


class HPOAnnotations(Source):
    """
    The [Human Phenotype Ontology](http://human-phenotype-ontology.org) group
    curates and assembles over 115,000 annotations to hereditary diseases
    using the HPO ontology. Here we create OBAN-style associations
    between diseases and phenotypic features, together with their evidence,
    and age of onset and frequency (if known).
    The parser currently only processes the "abnormal" annotations.
    Association to "remarkable normality" will be added in the near future.

    We create additional associations from text mining.  See info at
    http://pubmed-browser.human-phenotype-ontology.org/.

    Also, you can read about these annotations in
    [PMID:26119816](http://www.ncbi.nlm.nih.gov/pubmed/26119816).

    In order to properly test this class,
    you should have a conf.json file configured with some test ids, in
    the structure of:
    # as examples.  put your favorite ids in the config.
    <pre>
    test_ids: {"disease" : ["OMIM:119600", "OMIM:120160"]}
    </pre>

    """

    files = {
        'annot': {
            'file': 'phenotype_annotation.tab',
            'url': HPOADL + '/phenotype_annotation.tab'},
        'version': {
            'file': 'data_version.txt',
            'url': HPOADL + '/data_version.txt'},
        # 'neg_annot': {
        #   'file': 'phenotype_annotation.tab',
        #    'url': HPOADL + '/negative_phenotype_annotation.tab'},
        'doid': {
            'file': 'doid.owl',
            'url': 'http://purl.obolibrary.org/obo/doid.owl'
        }
    }

    # note, two of these codes are awaiting term requests.  see #114 and
    # https://code.google.com/p/evidenceontology/issues/detail?id=32
    # TODO TEC see if the GC issue translated into a GH issue
    eco_dict = {
        # FIXME currently using "curator inference used in manual assertion"
        "ICE": "ECO:0000305",
        # Inferred from Electronic Annotation
        "IEA": "ECO:0000501",
        # FIXME currently is"experimental evidence used in manual assertion"
        "PCS": "ECO:0000269",
        # Traceable Author Statement
        "TAS": "ECO:0000304",
        # FIXME currently using computational combinatorial evidence
        # in automatic assertion
        "ITM": "ECO:0000246",
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'hpoa')

        self.dataset = Dataset(
            'hpoa', 'Human Phenotype Ontology',
            'http://www.human-phenotype-ontology.org', None,
            'http://www.human-phenotype-ontology.org/contao/index.php/legal-issues.html')

        self.replaced_id_count = 0

        if 'test_ids' not in config.get_config()\
                or 'disease' not in config.get_config()['test_ids']:
            logger.warning("not configured with disease test ids.")
            self.test_ids = []
        else:
            self.test_ids = config.get_config()['test_ids']['disease']

        # data-source specific warnings to be removed when issues are cleared
        logger.warning(
            "note that some ECO classes are missing for ICE, PCS, and ITM;" +
            " using temporary mappings.")

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)

        self.scrub()

        # get the latest build from jenkins

        # use the files['version'] file as the version
        fname = '/'.join((self.rawdir, self.files['version']['file']))

        with open(fname, 'r', encoding="utf8") as f:
            # 2015-04-23 13:01
            v = f.readline()  # read the first line (the only line, really)
            d = datetime.strptime(
                v.strip(), '%Y-%m-%d %H:%M').strftime("%Y-%m-%d-%H-%M")
        f.close()

        st = os.stat(fname)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        # this will cause two dates to be attached to the dataset
        # (one from the filedate, and the other from here)
        # TODO when #112 is implemented,
        # this will result in only the whole dataset being versioned
        self.dataset.setVersion(filedate, d)

        self.get_common_files()

        return

    def scrub(self):
        """
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:
        * revise errors in identifiers for some OMIM and PMIDs

        :return: None

        """

        # scrub file of the oddities...lots of publication rewriting
        f = '/'.join((self.rawdir, self.files['annot']['file']))
        logger.info('scrubbing PubMed:12345 --> PMID:12345')
        pysed.replace(r'PubMed:', 'PMID:', f)

        logger.info('scrubbing pmid:12345 --> PMID:12345')
        pysed.replace(r'pmid:', 'PMID:', f)

        logger.info('scrubbing PMID:    12345 --> PMID:12345')
        pysed.replace(r'PMID:  *', 'PMID:', f)

        logger.info('scrubbing PMID12345 --> PMID:12345')
        pysed.replace(r'PMID([0-9][0-9]*)', r'PMID:\1', f)

        logger.info('scrubbing MIM12345 --> OMIM:12345')
        pysed.replace(r'MIM([0-9][0-9]*)', r'OMIM:\1', f)

        logger.info('scrubbing MIM:12345 --> OMIM:12345')
        pysed.replace(r";MIM", ";OMIM", f)

        logger.info('scrubbing ORPHANET --> Orphanet')
        pysed.replace("ORPHANET", "Orphanet", f)

        logger.info('scrubbing ORPHA --> Orphanet')
        pysed.replace("ORPHA", "Orphanet", f)
        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows", limit)

        self.add_common_files_to_file_list()

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        # rare disease-phenotype associations
        self._process_phenotype_tab('/'.join((self.rawdir,
                                              self.files['annot']['file'])),
                                    limit)

        # TODO add negative phenotype statements #113
        # self._process_negative_phenotype_tab(self.rawfile,self.outfile,limit)

        # common disease-phenotype associations from text mining work
        self.process_all_common_disease_files(limit)

        logger.info("Finished parsing.")

        return

    def _map_evidence_to_codes(self, code_string):
        """
        A simple mapping of the code_string to it's ECO class
        using the dictionary defined here
        Currently includes ICE, IEA, PCS, TAS
        :param code_string:
        :return:

        """
        return self.eco_dict.get(code_string)

    def _process_phenotype_tab(self, raw, limit):
        """
        see info on format here:
        http://www.human-phenotype-ontology.org/contao/index.php/annotation-guide.html

        :param raw:
        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                row = [str(col).strip() for col in row]
                (db, num, name, qual, pheno_id, publist, eco, onset, freq, w,
                 asp, syn, date, curator) = row
                disease_id = db + ":" + num

                if self.testMode:
                    try:
                        id_list = self.test_ids
                        if id_list is None \
                                or disease_id not in id_list:
                            continue
                    except AttributeError:
                        continue

                # logger.info('adding %s', disease_id)

                model.addClassToGraph(disease_id, None)
                model.addClassToGraph(pheno_id, None)
                eco_id = self._map_evidence_to_codes(eco)
                model.addClassToGraph(eco_id, None)
                if onset is not None and onset != '':
                    model.addClassToGraph(onset, None)

                # we want to do things differently depending on
                # the aspect of the annotation
                # TODO PYLINT Redefinition of assoc type from
                #   dipper.models.assoc.D2PAssoc.D2PAssoc to
                #   dipper.models.assoc.DispositionAssoc.DispositionAssoc
                if asp == 'O' or asp == 'M':  # organ abnormality or mortality
                    assoc = D2PAssoc(
                        g, self.name, disease_id, pheno_id, onset, freq)
                elif asp == 'I':  # inheritance patterns for the whole disease
                    assoc = DispositionAssoc(
                        g, self.name, disease_id, pheno_id)
                elif asp == 'C':  # clinical course / onset
                    assoc = DispositionAssoc(
                        g, self.name, disease_id, pheno_id)
                else:
                    logger.error("I don't know what this aspect is: %s", asp)

                assoc.add_evidence(eco_id)

                publist = re.split(r'[,;]', publist)
                # blow these apart if there is a list of pubs
                for pub in publist:
                    pub = pub.strip()
                    pubtype = None
                    if pub != '':
                        # if re.match(
                        #       r'http://www.ncbi.nlm.nih.gov/bookshelf/br\.fcgi\?book=gene',
                        #        pub):
                        #     #http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=gene&part=ced
                        #     m = re.search(r'part\=(\w+)', pub)
                        #     pub_id = 'GeneReviews:'+m.group(1)
                        # elif re.search(
                        #        r'http://www.orpha.net/consor/cgi-bin/OC_Exp\.php\?lng\=en\&Expert\=',
                        #        pub):
                        #     m = re.search(r'Expert=(\d+)', pub)
                        #     pub_id = 'Orphanet:'+m.group(1)

                        if re.match(r'(PMID|ISBN-13|ISBN-10|ISBN|HPO)', pub):
                            if re.match(r'PMID', pub):
                                pubtype = \
                                    Reference.ref_types['journal_article']
                            elif re.match(r'HPO', pub):
                                pubtype = Reference.ref_types['person']
                            else:
                                pubtype = Reference.ref_types['publication']
                            r = Reference(g, pub, pubtype)
                            r.addRefToGraph()
                        elif re.match(r'(OMIM|Orphanet|DECIPHER)', pub):
                            # make the pubs a reference to the website,
                            # instead of the curie
                            if re.match(r'OMIM', pub):
                                omimnum = re.sub(r'OMIM:', '', pub)
                                omimurl = '/'.join(('http://omim.org/entry',
                                                    str(omimnum).strip()))
                                pub = omimurl
                            elif re.match(r'Orphanet:', pub):
                                orphanetnum = re.sub(r'Orphanet:', '', pub)
                                orphaneturl = \
                                    ''.join((
                                        'http://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=',
                                        str(orphanetnum)))
                                pub = orphaneturl
                            elif re.match(r'DECIPHER:', pub):
                                deciphernum = re.sub(r'DECIPHER:', '', pub)
                                decipherurl = '/'.join(
                                    ('https://decipher.sanger.ac.uk/syndrome',
                                     deciphernum))
                                pub = decipherurl
                            pubtype = Reference.ref_types['webpage']
                        elif re.match(r'http', pub):
                            pass
                        else:
                            logger.error('Unknown pub type for %s: %s',
                                         disease_id, pub)
                            print(disease_id, 'pubs:', str(publist))
                            continue

                        if pub is not None:
                            assoc.add_source(pub)

                        # TODO add curator

                assoc.add_association_to_graph()

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def get_common_files(self):
        """
        Fetch the raw hpo-annotation-data by cloning/pulling the
        [repository](https://github.com/monarch-initiative/hpo-annotation-data.git)
        These files get added to the files object,
        and iterated over separately.
        :return:

        """

        repo_dir = '/'.join((self.rawdir, 'git'))
        REMOTE_URL = \
            "git@github.com:monarch-initiative/hpo-annotation-data.git"
        HTTPS_URL = \
            "https://github.com/monarch-initiative/hpo-annotation-data.git"

        # TODO if repo doesn't exist, then clone otherwise pull
        if os.path.isdir(repo_dir):
            shutil.rmtree(repo_dir)

        logger.info("Cloning common disease files from %s", REMOTE_URL)
        try:
            Repo.clone_from(REMOTE_URL, repo_dir)
        except GitCommandError:
            # Try with https and if this doesn't work fail
            Repo.clone_from(HTTPS_URL, repo_dir)

        return

    def add_common_files_to_file_list(self):
        repo_dir = '/'.join((self.rawdir, 'git'))
        common_disease_dir = '/'.join((repo_dir, 'common-diseases'))

        # add the files to the self.files object
        filelist = os.listdir(common_disease_dir)
        fcount = 0
        for f in filelist:
            if not re.search(r'\.tab', f):
                continue
            fcount += 1
            self.files['common'+str(fcount).zfill(7)] = {
                'file': '/'.join((common_disease_dir, f)),
                # TODO add url to reference the file?
                # need to get git version somehow?
            }
            # TODO add this to the dataset
        logger.info("Found %d common disease files", fcount)

        return

    def process_all_common_disease_files(self, limit=None):
        """
        Loop through all of the files that we previously fetched from git,
        creating the disease-phenotype assoc.
        :param limit:
        :return:

        """

        self.replaced_id_count = 0
        unpadded_doids = self.get_doid_ids_for_unpadding()
        total_processed = 0
        logger.info("Iterating over all common disease files")
        common_file_count = 0
        for f in self.files:
            if not re.match(r'common', f):
                continue
            common_file_count += 1
            raw = self.files[f]['file']
            total_processed += self.process_common_disease_file(
                raw, unpadded_doids, limit)
            if not self.testMode \
                    and limit is not None and total_processed > limit:
                break
        logger.info("Finished iterating over all common disease files.")
        logger.info("Fixed %d/%d incorrectly zero-padded ids",
                    self.replaced_id_count, common_file_count)
        return

    def get_doid_ids_for_unpadding(self):
        """
        Here, we fetch the doid owl file, and get all the doids.
        We figure out which are not zero-padded, so we can map the DOID
        to the correct identifier when processing the common annotation files.

        This may become obsolete when
        https://github.com/monarch-initiative/hpo-annotation-data/issues/84
        is addressed.

        :return:

        """

        logger.info("Building list of non-zero-padded DOIDs")
        raw_file = '/'.join((self.rawdir, self.files['doid']['file']))
        doids = set()
        # scan the file and get all doids
        with open(raw_file, 'r', encoding="utf8") as f:
            for line in f:
                matches = re.search(r'(DOID_\d+)', line)
                if matches is not None:
                    for m in matches.groups():
                        doids.add(re.sub(r'_', ':', m))

        nopad_doids = set()
        for d in doids:
            num = re.sub(r'DOID[:_]', '', d)
            # look for things not starting with zero
            if not re.match(r'^0', str(num)):
                nopad_doids.add(num)

        logger.info("Found %d/%d DOIDs are not zero-padded",
                    len(nopad_doids), len(doids))

        return nopad_doids

    def process_common_disease_file(self, raw, unpadded_doids, limit=None):
        """
        Make disaese-phenotype associations.
        Some identifiers need clean up:
        * DOIDs are listed as DOID-DOID: --> DOID:
        * DOIDs may be unnecessarily zero-padded.
        these are remapped to their non-padded equivalent.

        :param raw:
        :param unpadded_doids:
        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        assoc_count = 0
        replace_id_flag = False

        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            header = csvfile.readline()  # skip the header row
            logger.info("HEADER: %s", header)
            disease_id = None
            for row in filereader:

                if 21 == len(row):
                    (did, dname, gid, gene_name, genotype, gene_symbols,
                     phenotype_id, phenotype_name, age_of_onset_id,
                     age_of_onset_name, eid, evidence_name, frequency, sex_id,
                     sex_name, negation_id, negation_name, description,
                     pub_ids, assigned_by,
                     date_created) = [str(col).strip() for col in row]
                else:
                    logger.warning(
                        "Wrong number of columns! expected 21, got: %s in: %s",
                        len(row), raw)
                    logger.warning("%s", row)
                    continue
                # b/c "PMID:    17223397"
                pub_ids = re.sub(r'  *', '', pub_ids)

                disease_id = re.sub(r'DO(ID)?[-\:](DOID:)?', 'DOID:', did)
                disease_id = re.sub(r'MESH-', 'MESH:', disease_id)
                if not re.search(r'(DOID\:|MESH\:\w)\d+', disease_id):
                    logger.warning("Invalid id format: %s", disease_id)

                # figure out if the doid should be unpadded,
                # then use the unpadded version instead
                if re.match(r'DOID', disease_id):
                    unpadded_num = re.sub(r'DOID:', '', disease_id)
                    unpadded_num = unpadded_num.lstrip('0')
                    if unpadded_num in unpadded_doids:
                        fixed_id = 'DOID:' + unpadded_num
                        replace_id_flag = True
                        disease_id = fixed_id.strip()

                if self.testMode and disease_id not in self.test_ids:
                    # since these are broken up into disease-by-disease,
                    # just skip the whole file
                    return 0
                else:
                    line_counter += 1

                if negation_id != '':
                    continue  # TODO add negative associations

                if disease_id != '' and phenotype_id != '':
                    assoc = D2PAssoc(
                        g, self.name, disease_id, phenotype_id.strip())
                    if age_of_onset_id != '':
                        assoc.onset = age_of_onset_id
                    if frequency != '':
                        assoc.frequency = frequency
                    eco_id = self._map_evidence_to_codes(eid)
                    if eco_id is None:
                        eco_id = self._map_evidence_to_codes('ITM')
                    assoc.add_evidence(eco_id)
                    # TODO add sex? - not in dataset yet
                    if description != '':
                        assoc.set_description(description)
                    if pub_ids != '':
                        for p in pub_ids.split(';'):
                            p = re.sub(r'  *', '', p)
                            if re.search(r'(DOID|MESH)', p) \
                                    or re.search(r'Disease name contained',
                                                 description):
                                # skip "pubs" that are derived from
                                # the classes themselves
                                continue
                            assoc.add_source(p.strip())
                    # TODO assigned by?

                    assoc.add_association_to_graph()
                    assoc_count += 1

                if not self.testMode and limit is not None\
                        and line_counter > limit:
                    break

            if replace_id_flag:
                logger.info("replaced DOID with unpadded version")
                self.replaced_id_count += 1
            logger.info(
                "Added %d associations for %s.", assoc_count, disease_id)
        return assoc_count

    def getTestSuite(self):
        import unittest
        from tests.test_hpoa import HPOATestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(HPOATestCase)

        return test_suite
