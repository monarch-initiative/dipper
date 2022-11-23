import csv
import os
import logging
import re
import shutil
import tarfile
import glob
from datetime import datetime
from stat import ST_CTIME

import requests
from dipper.sources.Source import Source
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv
from dipper import config

LOG = logging.getLogger(__name__)

# summer 2018 PR mentioned switching to the new format
PNR = 'http://compbio.charite.de/jenkins/job/hpo.annotations.current'
HPOADL2 = PNR + '/lastSuccessfulBuild/artifact/current'
              # '/lastSuccessfulBuild/artifact/misc_2018'   disappeared by 2020 03

# Fall 2018 CM made a mondo version of common disease (but we decided to hold off now)
GITAPI = "https://api.github.com/repos/monarch-initiative"

ORPHANET = 'http://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert='
CONF = config.get_config()


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
    you should have a resources/test_ids.yaml file configured with some test ids,
    in the structure of:
    # as examples.  put your favorite ids in the config.
    <pre>
    test_ids: {"disease" : ["OMIM:119600", "OMIM:120160"]}
    </pre>

    """

    files = {
        'hpoa': {
            'file': 'phenotype.hpoa',
            'url':  HPOADL2 + '/phenotype.hpoa',
            'columns': [
                '#DatabaseID',
                'DiseaseName',
                'Qualifier',
                'HPO_ID',
                'Reference',
                'Evidence',
                'Onset',
                'Frequency',
                'Sex',
                'Modifier',
                'Aspect',
                'Biocuration'
            ]
        },
        'doid': {  # DOID plans to remain all about inconsistant left zero padding
            'file': 'doid.owl',
            'url': 'http://purl.obolibrary.org/obo/doid.owl'
        }
    }

    small_files = {
        # this is a placeholder for the columns in the "common-disease"
        # small file format which were added to the 'files' dict programaticly
        'columns': [
            'Disease ID',
            'Disease Name',
            'Gene ID',
            'Gene Name',
            'Genotype',
            'Gene Symbol(s)',
            'Phenotype ID',
            'Phenotype Name',
            'Age of Onset ID',
            'Age of Onset Name',
            'Evidence ID',
            'Evidence Name',
            'Frequency',
            'Sex ID',
            'Sex Name',
            'Negation ID',
            'Negation Name',
            'Description',
            'Pub',
            'Assigned by',
            'Date Created'
        ]
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='hpoa',
            ingest_title='Human Phenotype Ontology',
            ingest_url='https://hpo.jax.org/app/',
            ingest_logo='source-hpo.png',
            license_url=None,
            data_rights='https://hpo.jax.org/app/license',
            # file_handle=None
        )
        self.dataset.set_citation('https://hpo.jax.org/app/citation')
        self.replaced_id_count = 0
        self.test_ids = self.all_test_ids['disease']

        return

    def fetch(self, is_dl_forced=False):
        # get the latest build from jenkins
        self.get_files(is_dl_forced)
        # and git repo
        # self.get_common_files()
        return

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows", limit)

        # self.add_common_files_to_file_list()
        LOG.info("Parsing files...")
        if self.test_only:
            self.test_mode = True

        self._process_phenotype_hpoa(self.files['hpoa'], limit=limit)

        # TODO add negative phenotype statements #113
        # self._process_negative_phenotype_tab(self.rawfile,self.outfile,limit)

        # common disease-phenotype associations from text mining work
        # self.process_all_common_disease_files(limit)

        LOG.info("Finished parsing.")

        return

    def _process_phenotype_hpoa(self, file_info, limit=None):
        """
        see info on format here:
        http://www.human-phenotype-ontology.org/contao/index.php/annotation-guide.html

        :param raw:
        :param limit:
        :return:

        """
        src_key = 'hpoa'

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)

        raw = '/'.join((self.rawdir, file_info['file']))

        # this will cause two dates to be attached to the dataset
        # (one from the filedate, and the other from here)
        # TODO when #112 is implemented,
        # this will result in only the whole dataset being versioned

        col = self.files[src_key]['columns']
        with open(raw, 'r', encoding="utf8") as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t', quotechar='\"')
            row = next(reader)  # drop Description
            row = str(next(reader))[9:19]
            LOG.info("Ingest from %s", row)
            date = datetime.strptime(
                row.strip(), '%Y-%m-%d').strftime("%Y-%m-%d-%H-%M")

            if file_info.get("url") is not None:
                self.dataset.set_ingest_source_file_version_date(
                    file_info.get("url"), date)

            row = next(reader)  # drop tracker url
            row = next(reader)  # drop release url
            row = next(reader)  # headers
            # row[0] = row[0][1:]  # uncomment;  but not allways needed ?!
            if not self.check_fileheader(col, row):
                pass

            for row in reader:
                row = [str(col).strip() for col in row]

                disease_id = row[col.index('#DatabaseID')]
                # 98246 OMIM
                # 68646 ORPHA
                # 297 DECIPHER

                if self.test_mode:
                    try:
                        id_list = self.test_ids
                        if id_list is None or disease_id not in id_list:
                            continue
                    except AttributeError:
                        continue

                # row[col.index('DiseaseName')]  unused

                if row[col.index('Qualifier')] == 'NOT':
                    continue

                hpo_id = row[col.index('HPO_ID')]
                publist = row[col.index('Reference')]
                eco_id = self.resolve(row[col.index('Evidence')])
                onset = row[col.index('Onset')]
                freq = row[col.index('Frequency')]
                sex = row[col.index('Sex')].lower()
                # row[col.index('Modifier')]   unused
                asp = row[col.index('Aspect')]
                # row[col.index('Biocuration')]  unused

                # LOG.info(
                #    'adding <%s>-to-<%s> because <%s>', disease_id, hpo_id, eco_id)

                model.addClassToGraph(disease_id)
                model.addClassToGraph(eco_id)
                if onset is not None and onset != '':
                    model.addClassToGraph(onset)

                if asp in ('P', 'M'):  # phenotype? abnormality or mortality
                    model.addClassToGraph(hpo_id)
                    assoc = D2PAssoc(  # default rel=self.globaltt['has phenotype']
                        graph, self.name, disease_id, hpo_id, onset, freq
                    )
                elif asp in ('I', 'C'):  # inheritance pattern or clinical course/onset
                    model.addClassToGraph(hpo_id)
                    assoc = D2PAssoc(
                        graph,
                        self.name,
                        disease_id,
                        hpo_id,
                        rel=self.globaltt['has disposition']
                    )
                else:
                    LOG.error("Unknown aspect : %s at line %i", asp, reader.line_num)

                assoc.add_evidence(eco_id)
                if sex is not None and sex != '':
                    self.graph.addTriple(
                        assoc.get_association_id(),
                        self.globaltt['has_sex_specificty'],
                        self.globaltt[sex],
                        object_category=blv.terms['BiologicalSex']
                    )

                # Publication
                # cut -f 5 phenotype.hpoa | grep ";" | tr ';' '\n' | cut -f1 -d ':' |\
                # sort | uniq -c | sort -nr
                # 629 PMID
                # 63 OMIM
                # 42 ISBN-13
                # 36 http

                for pub in publist.split(';'):
                    pub = pub.strip()

                    # there have been several malformed PMIDs
                    if pub[:4] != 'http' and \
                            graph.curie_regexp.fullmatch(pub) is None:
                        LOG.warning(
                            'Record %s has a malformed Reference %s', disease_id, pub)
                        continue

                    pubtype = None

                    if pub[:5] == 'PMID:':
                        pubtype = self.globaltt['journal article']

                    elif pub[:4] == 'ISBN':
                        pubtype = self.globaltt['publication']

                    elif pub[:5] == 'OMIM:':
                        pub = 'http://omim.org/entry/' + pub[5:]
                        pubtype = self.globaltt['web page']

                    elif pub[:9] == 'DECIPHER:':
                        pubtype = self.globaltt['web page']

                    elif pub[:6] == 'ORPHA:':
                        pubtype = self.globaltt['web page']

                    elif pub[:4] == 'http':
                        pubtype = self.globaltt['web page']

                    else:
                        LOG.error(
                            'Unknown pub type for disease %s from "%s"',
                            disease_id, pub)
                        continue

                    if pub is not None:
                        assoc.add_source(pub)
                        if pubtype is not None:
                            ref = Reference(graph, pub, pubtype)
                            # ref.setTitle('');  ref.setYear()

                            ref.addRefToGraph()
                    # TODO add curator

                    # pprint.pprint(assoc)

                    assoc.add_association_to_graph()

                if not self.test_mode and limit is not None and reader.line_num > limit:
                    break
        return

    def get_common_files(self):
        """
        Fetch the hpo-annotation-data
        [repository](https://github.com/monarch-initiative/hpo-annotation-data.git)
        as a tarball

        :return:

        """
        # curl -sLu "username:personal-acess-token" \
        # GITAPI + "/hpo-annotation-data/tarball/master" > hpoa.tgz

        repo_dir = self.rawdir + '/git'
        username = CONF['user']['hpoa']
        response = requests.get(
            GITAPI + '/hpo-annotation-data/tarball/master',
            auth=(username, CONF['keys'][username]))

        with open(self.rawdir + '/hpoa.tgz', 'wb') as tgz:
            tgz.write(response.content)

        if os.path.isdir(repo_dir):  # scoarched earth approach
            shutil.rmtree(repo_dir)
        os.mkdir(repo_dir)

        with tarfile.open('raw/hpoa/hpoa.tgz', 'r') as tarball:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tarball, repo_dir)

        # TO-DO add this to the dataset object
        # hmm ...kay... have git commit-hash in tarball repo name

        repo_hash = glob.glob(
            str(
                '/'.join(
                    (self.rawdir, 'git', 'monarch-initiative-hpo-annotation-data-*')
                )))[-42:]

        print(repo_hash)
        repo_hash = str(repo_hash)
        # (note this makes little sense as it is a private repo)
        self.dataset.set_ingest_source(
            '/'.join((
                'https://github.com/monarch-initiative/hpo-annotation-data/tree',
                repo_hash)))

        return

    def add_common_files_to_file_list(self):
        '''
            The (several thousands) common-disease files from the repo tarball
            are added to the files object.
            try adding the 'common-disease-mondo' files as well?

        '''
        repo_dir = '/'.join((self.rawdir, 'git'))
        common_disease_dir = '/'.join((
            repo_dir,
            'monarch-initiative-hpo-annotation-*', 'common-diseases-mondo/*.tab'))

        # add the files to the self.files object
        filelist = glob.glob(common_disease_dir)
        fcount = 0
        for small_file in filelist:
            if small_file[-4:] == '.tab':
                fcount += 1
                self.files[
                    'common' + str(fcount).zfill(7)] = {
                        'file': '/'.join((common_disease_dir, small_file)),
                    }
        LOG.info("Found %d common disease files", fcount)

        return

    def process_all_common_disease_files(self, limit=None):
        """
        Loop through all of the files that we previously fetched from git,
        creating the disease-phenotype association.
        :param limit:
        :return:

        """
        LOG.info("Iterating over all common disease files")
        common_file_count = 0
        total_processed = ""  # stopgap gill we fix common-disease files
        unpadded_doids = ""   # stopgap gill we fix common-disease files
        for ingest in self.files:
            if ingest[:5] == 'common':
                common_file_count += 1
                raw = self.files[ingest]['file']
                total_processed += self.process_common_disease_file(
                    raw, unpadded_doids, limit)
            if not self.test_mode and limit is not None and total_processed > limit:
                break
        LOG.info("Finished iterating over all common disease files.")
        return

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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        assoc_count = 0
        replace_id_flag = False
        col = self.small_files['columns']

        with open(raw, 'r', encoding="utf8") as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t', quotechar='\"')
            header = tsvfile.readline()
            if header != col:
                LOG.error("HEADER: has changed in %s.", raw)
                raise ValueError(col - header)

            disease_id = None
            for row in reader:
                row = [str(x).strip() for x in row]

                did = row[col.index('Disease ID')]
                # genotype = row[col.index('Genotype')]
                phenotype_id = row[col.index('Phenotype ID')]
                age_of_onset_id = row[col.index('Age of Onset ID')]
                eid = row[col.index('Evidence ID')]
                frequency = row[col.index('Frequency')]
                negation_id = row[col.index('Negation ID')]
                description = row[col.index('Description')]
                pub_ids = row[col.index('Pub')]

                disease_id = re.sub(r'DO(ID)?[-\:](DOID:)?', 'DOID:', did)
                disease_id = re.sub(r'MESH-', 'MESH:', disease_id)

                if not re.search(r'(DOID\:|MESH\:\w)\d+', disease_id):
                    LOG.warning("Invalid id format: %s", disease_id)

                # figure out if the doid should be unpadded,
                # then use the unpadded version instead
                if re.match(r'DOID', disease_id):
                    unpadded_num = re.sub(r'DOID:', '', disease_id)
                    unpadded_num = unpadded_num.lstrip('0')
                    if unpadded_num in unpadded_doids:
                        fixed_id = 'DOID:' + unpadded_num
                        replace_id_flag = True
                        disease_id = fixed_id.strip()

                if self.test_mode and disease_id not in self.test_ids:
                    # since these are broken up into disease-by-disease,
                    # just skip the whole file
                    return 0

                if negation_id != '':
                    continue  # TODO add negative associations

                if disease_id != '' and phenotype_id != '':
                    assoc = D2PAssoc(
                        graph, self.name, disease_id, phenotype_id.strip())
                    if age_of_onset_id != '':
                        assoc.onset = age_of_onset_id
                    if frequency != '':
                        assoc.frequency = frequency
                    eco_id = self.localtt[eid]
                    if eco_id is None:
                        eco_id = self.localtt['ITM']

                    assoc.add_evidence(eco_id)
                    # TODO add sex? - not in dataset yet
                    if description != '':
                        assoc.set_description(description)
                    if pub_ids != '':
                        for pub in pub_ids.split(';'):
                            pub = re.sub(r'  *', '', pub)  # fixed now but just in case

                            # there have been several malformed PMIDs curies
                            if pub[:4] != 'http' and \
                                    graph.curie_regexp.fullmatch(pub) is None:
                                LOG.warning(
                                    'Record %s has a malformed Pub %s', did, pub)
                                continue

                            if re.search(
                                    r'(DOID|MESH)', pub) or re.search(
                                        r'Disease name contained', description):
                                # skip "pubs" that are derived from
                                # the classes themselves
                                continue
                            assoc.add_source(pub.strip())
                    # TODO assigned by?

                    assoc.add_association_to_graph()
                    assoc_count += 1

                if not self.test_mode and limit is not None\
                        and reader.line_num > limit:
                    break

            if replace_id_flag:
                LOG.info("replaced DOID with unpadded version")
                self.replaced_id_count += 1
            LOG.info(
                "Added %d associations for %s.", assoc_count, disease_id)
        return assoc_count

    def getTestSuite(self):
        import unittest
        from tests.test_hpoa import HPOATestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(HPOATestCase)

        return test_suite
