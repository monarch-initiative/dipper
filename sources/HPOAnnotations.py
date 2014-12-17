import csv, os, datetime
from datetime import datetime
from stat import *
import urllib


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
    ANNOT_URL = "http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab"
    ANNOT_FILE = "phenotype_annotation.tab"

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


    def __init__(self):
        Source.__init__(self, 'hpoa')
        self.outfile = self.outdir+'/'+self.name + ".ttl"
        self.rawfile = ('/').join((self.rawdir,self.ANNOT_FILE))
        self.datasetfile = self.outdir + '/' + self.name + '_dataset.ttl'

        print("Setting outfile to", self.outfile)
        self.curie_map = D2PAssoc.curie_map.copy()
        self.curie_map.update(DispositionAssoc.curie_map)
        self.curie_map.update(self.disease_prefixes)

        self.load_bindings()

        self.dataset = Dataset('hpoa', 'Human Phenotype Ontology', 'http://www.human-phenotype-ontology.org')

        print("WARN: note that some ECO classes are missing for ICE and PCS; using temporary mappings.")

        return

    def fetch(self):
        self.fetch_from_url(self.ANNOT_URL,self.rawfile)

        self.dataset.setFileAccessUrl(self.ANNOT_URL)

        st = os.stat(self.rawfile)
        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        #get the latest build from jenkins
        jenkins_info=eval(urllib.request.urlopen('http://compbio.charite.de/hudson/job/hpo.annotations/lastSuccessfulBuild/api/python').read())
        version=jenkins_info['number']
        self.dataset.setVersion(filedate,str(version))


        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.curie_map.keys():
            v=self.curie_map[k]
            self.graph.bind(k, Namespace(v))
#            print("bound ", k, " to ", v)
        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    #supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        Source.parse(self)
        if (limit is not None):
            print("Only parsing first", limit, "rows")
        line_counter = 0
        #todo move this into super

        self._process_phenotype_tab(self.rawfile,self.outfile,self.graph,limit)
        Assoc().loadObjectProperties(self.graph)

        #TODO add negative phenotype statements
        # http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/negative_phenotype_annotation.tab
        #self._process_negative_phenotype_tab(self.rawfile,self.outfile,limit)

        filewriter = open(self.outfile, 'w')
        self.load_bindings()
        print("Finished parsing",self.rawfile, ". Writing turtle to",self.outfile)
        print(self.graph.serialize(format="turtle").decode(),file=filewriter)
        filewriter.close()

        filewriter = open(self.datasetfile,'w')
        print(self.dataset.getGraph().serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        print("Wrote", len(self.graph), "nodes")
        return

    def sanity_checks(self):
        #TODO syntactic checking in the file
        #will return False if any of the conditions are not satisfied

        #check if there's any evidence codes not in our dictionary

        #check if there's "aspect" codes not in our dictionary

        return True



    def _map_evidence_to_codes(self, code_string):
        #the following evidence codes are used as literals: ICE, IEA, PCS, TAS
        #here, we map them to the ECO


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

    def verify(self):
        status = True
        self._verify(self.outfile)
        self._verifyowl(self.outfile)
        #TODO verify some kind of relationship that should be in the file

        return status