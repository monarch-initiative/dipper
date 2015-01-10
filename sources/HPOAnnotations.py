import csv, os, datetime
from datetime import datetime
from stat import *
import urllib
from utils import pysed


from sources.Source import Source

from models.D2PAssoc import D2PAssoc
from models.DispositionAssoc import DispositionAssoc
from rdflib import Namespace
from models.Dataset import Dataset
from models.Assoc import Assoc

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
    disease_prefixes = {
        'OMIM' : 'http://purl.obolibrary.org/obo/OMIM_',
        'DECIPHER' : 'http://purl.obolibrary.org/obo/DECIPHER_',
        'ORPHANET' : 'http://purl.obolibrary.org/obo/ORPHANET_'
    }

    curie_map = {
        'HPO' : 'http://human-phenotype-ontology.org/' #used for HPO-person identifiers, but they don't resolve
    }

    def __init__(self):
        Source.__init__(self, 'hpoa')

        self.curie_map.update(D2PAssoc.curie_map)
        self.curie_map.update(DispositionAssoc.curie_map)
        self.curie_map.update(self.disease_prefixes)

        self.load_bindings()

        self.dataset = Dataset('hpoa', 'Human Phenotype Ontology', 'http://www.human-phenotype-ontology.org')

        #data-source specific warnings (will be removed when issues are cleared)
        print("WARN: note that some ECO classes are missing for ICE and PCS; using temporary mappings.")

        return

    def fetch(self):
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],('/').join((self.rawdir,file['file'])))
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
        print('INFO: scrubbing PubMed:12345 --> PMID:12345')
        pysed.replace("PubMed", 'PMID', f)

        print('INFO: scrubbing pmid:12345 --> PMID:12345')
        pysed.replace("pmid", 'PMID', f)

        print('INFO: scrubbing PMID12345 --> PMID:12345')
        pysed.replace("PMID([0-9][0-9]*)", 'PMID:\\1', f)

        print('INFO: scrubbing MIM12345 --> OMIM:12345')
        pysed.replace('MIM([0-9][0-9]*)', 'OMIM:\\1', f)

        print('INFO: scrubbing MIM:12345 --> OMIM:12345')
        pysed.replace(";MIM",";OMIM", f)
        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.curie_map.keys():
            v=self.curie_map[k]
            self.graph.bind(k, Namespace(v))
        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    #supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows")

        print("Parsing files...")
        self._process_phenotype_tab(('/').join((self.rawdir,self.files['annot']['file'])),self.outfile,self.graph,limit)
        Assoc().loadObjectProperties(self.graph)

        #TODO add negative phenotype statements
        #self._process_negative_phenotype_tab(self.rawfile,self.outfile,limit)

        print("Finished parsing.")

        print("Loaded", len(self.graph), "nodes")
        return


    def _map_evidence_to_codes(self, code_string):
        '''
        A simple mapping of the code_string to it's ECO class using the dictionary defined here
        Currently includes ICE, IEA, PCS, TAS
        :param code_string:
        :return:
        '''
        return self.eco_dict.get(code_string)

    def _process_phenotype_tab(self, raw, out, g, limit=None):
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (db, num, name, qual, pheno_id, publist, eco, onset, freq, w, asp, syn, date, curator) = row
                disease_id = db + ":" + num

                #blow these apart if there is a list of pubs
                publist=publist.split(';')
                for pub in publist:
                    pub = pub.strip()
                    assoc_id = self.make_id(db + num + qual + pheno_id + pub + eco + onset + freq + date + curator)
                    assoc = None
                    # we want to do things differently depending on the aspect of the annotation
                    if (asp == 'O' or asp == 'M'):  #organ abnormality or mortality
                        assoc = D2PAssoc(assoc_id, disease_id, pheno_id, onset, freq, pub, self._map_evidence_to_codes(eco), self.curie_map)
                        g = assoc.addAssociationNodeToGraph(g)
                    elif (asp == 'I'):  #inheritance patterns for the whole disease
                        assoc = DispositionAssoc(assoc_id, disease_id, pheno_id, pub, self._map_evidence_to_codes(eco), self.curie_map)
                        g = assoc.addAssociationNodeToGraph(g)
                    elif (asp == 'C'):  #clinical course / onset
                        #FIXME is it correct for these to be dispositions?
                        assoc = DispositionAssoc(assoc_id, disease_id, pheno_id, pub, self._map_evidence_to_codes(eco), self.curie_map)
                        g = assoc.addAssociationNodeToGraph(g)
                    else:
                        #TODO throw an error?
                        print("I don't know what this is:", asp)

                    #if (assoc is not None):
                    #   self.triple_count += assoc.printTypes(filewriter)
                    #    self.triple_count += assoc.printAssociation(filewriter)

                if (limit is not None and line_counter > limit):
                   break

        return

