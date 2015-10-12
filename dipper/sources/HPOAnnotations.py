import csv
import os
from datetime import datetime
from stat import *
import logging
import re

from dipper.utils import pysed
from dipper.utils.GraphUtils import GraphUtils
from dipper.sources.Source import Source
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.assoc.DispositionAssoc import DispositionAssoc
from dipper.models.Dataset import Dataset
from dipper.models.assoc.Association import Assoc
from dipper.models.Reference import Reference
from dipper import curie_map
from dipper import config


'''
#see info on format here:http://www.human-phenotype-ontology.org/contao/index.php/annotation-guide.html
#file schema:
#1 	DB 	required 	MIM
#2 	DB_Object_ID 	required 	154700
#3 	DB_Name 	required 	Achondrogenesis, type IB
#4 	Qualifier 	optional 	NOT
#5 	HPO ID 	required 	HP:0002487
#6 	DB:Reference 	required 	OMIM:154700 or PMID:15517394
#7 	Evidence code 	required 	IEA
#8 	Onset modifier  optional	HP:0003577
#9	Frequency modifier	optional	"70%" or "12 of 30" or from the vocabulary show in table 2
#10 With 	optional
#11 Aspect 	required 	O
#12 	Synonym 	optional 	ACG1B|Achondrogenesis, Fraccaro type
#13 	Date 	required 	YYYY.MM.DD
#14 	Assigned by 	required 	HPO
'''

logger = logging.getLogger(__name__)


class HPOAnnotations(Source):
    """
    The [Human Phenotype Ontology](http://human-phenotype-ontology.org) group curates and assembles
    over 115,000 annotations to hereditary diseases using the HPO ontology.
    Here we create OBAN-style associations between diseases and phenotypic features, together with their
    evidence, and age of onset and frequency (if known).
    The parser currently only processes the "abnormal" annotations.  Association to "remarkable normality"
    will be added in the near future.

    In order to properly test this class, you should have a conf.json file configured with some test ids, in
    the structure of:
        <pre>
        test_ids: {
            "disease" : ["OMIM:119600", "OMIM:120160"]  # as examples.  put your favorite ids in the config.
        }
        </pre>
    """

    files = {
        'annot': {'file' : 'phenotype_annotation.tab',
                   'url' : 'http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab'},
        'version': {'file' : 'data_version.txt',
                    'url' : 'http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/data_version.txt'},
#       'neg_annot': {'file' : 'phenotype_annotation.tab',
#                     'url' : 'http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/negative_phenotype_annotation.tab'
#        },
    }

    # note, two of these codes are awaiting term requests.  see #114 and
    # https://code.google.com/p/evidenceontology/issues/detail?id=32
    eco_dict = {
        "ICE": "ECO:0000305",  # FIXME currently using "curator inference used in manual assertion"
        "IEA": "ECO:0000501",  # Inferred from Electronic Annotation
        "PCS": "ECO:0000269",  # FIXME currently using "experimental evidence used in manual assertion"
        "TAS": "ECO:0000304"   # Traceable Author Statement
    }

    def __init__(self):
        Source.__init__(self, 'hpoa')

        self.load_bindings()

        self.dataset = Dataset('hpoa', 'Human Phenotype Ontology',
                               'http://www.human-phenotype-ontology.org', None,
                               'http://www.human-phenotype-ontology.org/contao/index.php/legal-issues.html')

        if 'test_ids' not in config.get_config() or 'disease' not in config.get_config()['test_ids']:
            logger.warn("not configured with disease test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['disease']

        # data-source specific warnings (will be removed when issues are cleared)
        logger.warn("note that some ECO classes are missing for ICE and PCS; using temporary mappings.")

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)

        self.scrub()

        # get the latest build from jenkins
        # NOT DOING THIS ANY MORE - but leaving it in for reference
        # jenkins_info = eval(urllib.request.urlopen('http://compbio.charite.de/hudson/job/hpo.annotations/lastSuccessfulBuild/api/python').read())
        # version = jenkins_info['number']

        # use the files['version'] file as the version
        fname = '/'.join((self.rawdir, self.files['version']['file']))

        with open(fname, 'r', encoding="utf8") as f:
            # 2015-04-23 13:01
            v = f.readline()  # read the first line (the only line, really)
            d = datetime.strptime(v.strip(), '%Y-%m-%d %H:%M').strftime("%Y-%m-%d-%H-%M")
        f.close()

        st = os.stat(fname)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        # this will cause two dates to be attached to the dataset (one from the filedate, and the other from here)
        # TODO when #112 is implemented, this will result in only the whole dataset being versioned
        self.dataset.setVersion(filedate, d)

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
        pysed.replace("PubMed", 'PMID', f)

        logger.info('scrubbing pmid:12345 --> PMID:12345')
        pysed.replace("pmid", 'PMID', f)

        logger.info('scrubbing PMID12345 --> PMID:12345')
        pysed.replace("PMID([0-9][0-9]*)", 'PMID:\\1', f)

        logger.info('scrubbing MIM12345 --> OMIM:12345')
        pysed.replace('MIM([0-9][0-9]*)', 'OMIM:\\1', f)

        logger.info('scrubbing MIM:12345 --> OMIM:12345')
        pysed.replace(";MIM", ";OMIM", f)

        logger.info('scrubbing ORPHANET --> Orphanet')
        pysed.replace("ORPHANET", "Orphanet", f)
        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_phenotype_tab('/'.join((self.rawdir, self.files['annot']['file'])), limit)

        # TODO add negative phenotype statements #113
        # self._process_negative_phenotype_tab(self.rawfile,self.outfile,limit)

        logger.info("Finished parsing.")

        return

    def _map_evidence_to_codes(self, code_string):
        """
        A simple mapping of the code_string to it's ECO class using the dictionary defined here
        Currently includes ICE, IEA, PCS, TAS
        :param code_string:
        :return:
        """
        return self.eco_dict.get(code_string)

    def _process_phenotype_tab(self, raw, limit):
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (db, num, name, qual, pheno_id, publist, eco, onset, freq, w, asp, syn, date, curator) = row
                disease_id = db + ":" + str(num)

                if self.testMode:
                    try:
                        id_list = self.test_ids
                        if id_list is None or disease_id.strip() not in id_list:
                            continue
                    except AttributeError:
                        continue

                # logger.info('adding %s', disease_id)

                gu.addClassToGraph(g, disease_id, None)
                gu.addClassToGraph(g, pheno_id, None)
                eco_id = self._map_evidence_to_codes(eco)
                gu.addClassToGraph(g, eco_id, None)
                if onset is not None and onset.strip() != '':
                    gu.addClassToGraph(g, onset, None)

                # we want to do things differently depending on the aspect of the annotation
                if asp == 'O' or asp == 'M':  # organ abnormality or mortality
                    assoc = D2PAssoc(self.name, disease_id, pheno_id, onset, freq)
                elif asp == 'I':  # inheritance patterns for the whole disease
                    assoc = DispositionAssoc(self.name, disease_id, pheno_id)
                elif asp == 'C':  # clinical course / onset
                    assoc = DispositionAssoc(self.name, disease_id, pheno_id)
                else:
                    logger.error("I don't know what this aspect is:", asp)

                assoc.add_evidence(eco_id)

                publist = publist.split(';')
                # blow these apart if there is a list of pubs
                for pub in publist:
                    pub = pub.strip()
                    if pub != '':
                        # if re.match('http://www.ncbi.nlm.nih.gov/bookshelf/br\.fcgi\?book=gene', pub):
                        #     #http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=gene&part=ced
                        #     m = re.search('part\=(\w+)', pub)
                        #     pub_id = 'GeneReviews:'+m.group(1)
                        # elif re.search('http://www.orpha.net/consor/cgi-bin/OC_Exp\.php\?lng\=en\&Expert\=', pub):
                        #     m = re.search('Expert=(\d+)', pub)
                        #     pub_id = 'Orphanet:'+m.group(1)
                        if not re.match('http', pub):
                            r = Reference(pub)
                            if re.match('PMID', pub):
                                r.setType(Reference.ref_types['journal_article'])
                            r.addRefToGraph(g)
                        # TODO add curator
                        assoc.add_source(pub)

                assoc.add_association_to_graph(g)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

            Assoc(None).load_all_properties(g)

        return

    def getTestSuite(self):
        import unittest
        from tests.test_hpoa import HPOATestCase
        # TODO add D2PAssoc tests

        test_suite = unittest.TestLoader().loadTestsFromTestCase(HPOATestCase)

        return test_suite
