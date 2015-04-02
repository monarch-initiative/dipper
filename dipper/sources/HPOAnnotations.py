import csv
import os
from datetime import datetime
from stat import *
import urllib
import logging

from dipper.utils import pysed
from dipper.utils.GraphUtils import GraphUtils
from dipper.sources.Source import Source
from dipper.models.D2PAssoc import D2PAssoc
from dipper.models.DispositionAssoc import DispositionAssoc
from dipper.models.Dataset import Dataset
from dipper.models.Assoc import Assoc
from dipper import curie_map


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

    files = {
        'annot' : {'file' : 'phenotype_annotation.tab',
                   'url' : 'http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab'
        },
#        'neg_annot' : {'file' : 'phenotype_annotation.tab',
#                   'url' : 'http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/negative_phenotype_annotation.tab'
#        },
    }

    #note, two of these codes are awaiting term requests
    #https://code.google.com/p/evidenceontology/issues/detail?id=32
    eco_dict = {
        "ICE": "ECO:0000305",  #FIXME currently using "curator inference used in manual assertion"
        "IEA": "ECO:0000501",  # Inferred from Electronic Annotation
        "PCS": "ECO:0000269",  #FIXME currently using "experimental evidence used in manual assertion"
        "TAS": "ECO:0000304"   #Traceable Author Statement
    }

    test_ids = ['OMIM:119600','OMIM:120160','OMIM:157140','OMIM:158900',
                'OMIM:166220','OMIM:168600','OMIM:219700','OMIM:253250',
                'OMIM:305900','OMIM:600669','OMIM:601278','OMIM:602421',
                'OMIM:605073','OMIM:607822',  #from coriell
                'Orphanet:99889','Orphanet:99','DECIPHER:14','DECIPHER:1',
                'OMIM:194072','OMIM:100100','Orphanet:99798']  #these ones xref the other disease DBs


    def __init__(self):
        Source.__init__(self, 'hpoa')

        self.load_bindings()

        self.dataset = Dataset('hpoa', 'Human Phenotype Ontology', 'http://www.human-phenotype-ontology.org')

        #data-source specific warnings (will be removed when issues are cleared)
        logger.warn("note that some ECO classes are missing for ICE and PCS; using temporary mappings.")

        return

    def fetch(self, is_dl_forced):
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                ('/').join((self.rawdir,file['file'])),
                                is_dl_forced)
            self.dataset.setFileAccessUrl(file['url'])
            # zfin versions are set by the date of download.
            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        self.scrub()

        #get the latest build from jenkins
        jenkins_info=eval(urllib.request.urlopen('http://compbio.charite.de/hudson/job/hpo.annotations/lastSuccessfulBuild/api/python').read())
        version=jenkins_info['number']

        self.dataset.setVersion(filedate,str(version))

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:
        * revise errors in identifiers for some OMIM and PMIDs
        :return: None
        '''
        # scrub file of the oddities...lots of publication rewriting
        f = ('/').join((self.rawdir,self.files['annot']['file']))
        logger.info('scrubbing PubMed:12345 --> PMID:12345')
        pysed.replace("PubMed", 'PMID', f)

        logger.info('scrubbing pmid:12345 --> PMID:12345')
        pysed.replace("pmid", 'PMID', f)

        logger.info('scrubbing PMID12345 --> PMID:12345')
        pysed.replace("PMID([0-9][0-9]*)", 'PMID:\\1', f)

        logger.info('scrubbing MIM12345 --> OMIM:12345')
        pysed.replace('MIM([0-9][0-9]*)', 'OMIM:\\1', f)

        logger.info('scrubbing MIM:12345 --> OMIM:12345')
        pysed.replace(";MIM",";OMIM", f)

        logger.info('scrubbing ORPHANET --> Orphanet')
        pysed.replace("ORPHANET","Orphanet", f)
        return


    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    #supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %s rows", limit)

        logger.info("Parsing files...")
        for testMode in [True,False]:
            self._process_phenotype_tab(('/').join((self.rawdir,self.files['annot']['file'])),limit,testMode)

        for g in [self.graph,self.testgraph]:
            Assoc().loadAllProperties(g)

        #TODO add negative phenotype statements
        #self._process_negative_phenotype_tab(self.rawfile,self.outfile,limit)

        logger.info("Finished parsing.")

        logger.info("Loaded %d nodes", len(self.graph))
        return


    def _map_evidence_to_codes(self, code_string):
        '''
        A simple mapping of the code_string to it's ECO class using the dictionary defined here
        Currently includes ICE, IEA, PCS, TAS
        :param code_string:
        :return:
        '''
        return self.eco_dict.get(code_string)

    def _process_phenotype_tab(self, raw, limit, testMode):
        if (testMode):
            g = self.testgraph
            logging.info("Processing testset.")
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (db, num, name, qual, pheno_id, publist, eco, onset, freq, w, asp, syn, date, curator) = row
                disease_id = db + ":" + num

                if testMode and (disease_id not in self.test_ids):
                    continue

                #blow these apart if there is a list of pubs
                publist=publist.split(';')
                for pub in publist:
                    pub = pub.strip()
                    assoc_id = self.make_id(db + num + qual + pheno_id + pub + eco + onset + freq + date + curator)
                    assoc = None
                    eco_id = self._map_evidence_to_codes(eco)
                    #make sure to add the disease, phenotype, eco, as classes.
                    #pub as individual is taken care of in the association function
                    gu.addClassToGraph(g,disease_id,None)
                    gu.addClassToGraph(g,pheno_id,None)
                    gu.addClassToGraph(g,eco_id,None)

                    # we want to do things differently depending on the aspect of the annotation
                    if (asp == 'O' or asp == 'M'):  #organ abnormality or mortality
                        assoc = D2PAssoc(assoc_id, disease_id, pheno_id, onset, freq, pub, eco_id)
                        g = assoc.addAssociationNodeToGraph(g)
                    elif (asp == 'I'):  #inheritance patterns for the whole disease
                        assoc = DispositionAssoc(assoc_id, disease_id, pheno_id, pub, eco_id)
                        g = assoc.addAssociationNodeToGraph(g)
                    elif (asp == 'C'):  #clinical course / onset
                        #FIXME is it correct for these to be dispositions?
                        assoc = DispositionAssoc(assoc_id, disease_id, pheno_id, pub, eco_id)
                        g = assoc.addAssociationNodeToGraph(g)
                    else:
                        #TODO throw an error?
                        logger.error("I don't know what this aspect is:", asp)


                if not testMode and (limit is not None and line_counter > limit):
                   break

        return

