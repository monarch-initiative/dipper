import csv
import os
from datetime import datetime
from stat import *
import re
import psycopg2


from rdflib import Literal
from rdflib.namespace import RDFS, OWL, RDF, DC, FOAF
from rdflib import Namespace, URIRef, BNode, Graph

from utils import pysed
from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from utils.CurieUtil import CurieUtil
import config
import curie_map
from utils.GraphUtils import GraphUtils


class EOM(Source):
    '''
    Be sure to have pg user/password connection details in your conf.json file, like:
      dbauth : {
        'disco' : {'user' : '<username>', 'password' : '<password>'}
      }
    '''
    tables = [
        'dv.nlx_157874_1'
    ]

    terms = {
        'phenotype': 'MONARCH:phenotype'  # Is this correct? What about GENO:0000348 - phenotype? MONARCH:phenotype

    }

    relationship = {
        'is_mutant_of': 'GENO:0000440',
        'derives_from': 'RO:0001000',
        'has_alternate_part': 'GENO:0000382',
        'has_reference_part': 'GENO:0000385',
        'in_taxon': 'RO:00002162',
        'has_zygosity': 'GENO:0000608',
        'is_sequence_variant_instance_of': 'GENO:0000408',
        'is_reference_instance_of': 'GENO:0000610',
        'hasRelatedSynonym': 'OIO:hasRelatedSynonym',
        'has_disposition': 'GENO:0000208',
        'has_phenotype': 'RO:0002200',
        'has_part': 'BFO:0000051',
        'has_variant_part': 'GENO:0000382'
    }


    def __init__(self):
        Source.__init__(self, 'eom')
        self.namespaces.update(curie_map.get())

        #update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('eom', 'EOM', 'http://elementsofmorphology.nih.gov')

        #check if config exists; if it doesn't, error out and let user know
        if (not (('dbauth' in config.get_config()) and ('disco' in config.get_config()['dbauth']))):
            print("ERROR: not configured with PG user/password.")

        #source-specific warnings.  will be cleared when resolved.
        #print("WARN: we are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self, is_dl_forced):

        #create the connection details for DISCO
        cxn = config.get_config()['dbauth']['disco']
        cxn.update({'host' : 'nif-db.crbs.ucsd.edu', 'database' : 'disco_crawler', 'port' : 5432 })

        self.dataset.setFileAccessUrl(('').join(('jdbc:postgresql://',cxn['host'],':',str(cxn['port']),'/',cxn['database'])))

        #process the tables
        #self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables,cxn)


        #FIXME: Fix
        datestamp=ver=None
        #get the resource version information from table mgi_dbinfo, already fetched above
        #outfile=('/').join((self.rawdir,'mgi_dbinfo'))
        '''
        if os.path.exists(outfile):
            st = os.stat(outfile)
            with open(outfile, 'r') as f:
                f.readline() #read the header row; skip
                info = f.readline()
                cols = info.split('\t')
                ver = cols[0] #col 0 is public_version
                ver = ver.replace('EOM ','')  #MGI 5.20 --> 5.20
                #MGI has a datestamp for the data within the database; use it instead of the download date
                #datestamp in the table: 2014-12-23 00:14:20
                d = cols[7].strip()  #modification date
                datestamp = datetime.strptime(d, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
                f.close()
        '''
        #datestamp = datetime.strptime(d, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
        st = os.stat(('/').join((self.rawdir,'dv.nlx_157874_1')))
        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")
        self.dataset.setVersion(filedate)






        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:

        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        #pysed.replace("\r", '', ('/').join((self.rawdir,'dv.nlx_157874_1')))

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")

        self._process_nlx_157874_1_view(('/').join((self.rawdir,'dv.nlx_157874_1')),limit)
        self._process_eom_terms(('/').join((self.rawdir,'eom_terms.tsv')),limit)

        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Found", len(self.graph), "nodes")
        return



    def _process_nlx_157874_1_view(self, raw, limit=None):
        '''
        This table contains the Elements of Morphology data that has been screen-scraped into DISCO.
        :param raw:
        :param limit:
        :return:
        '''

        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())
        line_counter = 0
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1

                (morphology_term_id, morphology_term_num, morphology_term_label, morphology_term_url,
                 terminology_category_label, terminology_category_url, subcategory, objective_definition,
                 subjective_definition, comments, synonyms, replaces, small_figure_url, large_figure_url,
                 e_uid, v_uid, v_uuid, v_last_modified) = line.split('\t')





                self.cm = curie_map.get()

                self.gu = GraphUtils(self.cm)
                self.cu = CurieUtil(self.cm)

                self.eom = Graph()

                #Add morphology_term_id as a class? An instance of what type? Phenotype, yes?
                self.graph.add((URIRef(cu.get_uri(morphology_term_id)),RDF['type'],URIRef(cu.get_uri(self.terms['phenotype']))))

                #morphology_term_id has label morphology_term_label
                self.graph.add((URIRef(cu.get_uri(morphology_term_id)),RDFS['label'],Literal(morphology_term_label)))

                #morphology_term_id has depiction small_figure_url
                if small_figure_url != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),FOAF['depiction'],Literal(small_figure_url)))
                #The below statement doesn't hang the image on the term id, so I used the above statement instead.
                #self.graph.add((Literal(small_figure_url),FOAF['depicts'],(URIRef(cu.get_uri(morphology_term_id)))))

                #morphology_term_id has depiction large_figure_url
                if large_figure_url != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),FOAF['depiction'],Literal(large_figure_url)))
                #The below statement doesn't hang the image on the term id, so I used the above statement instead.
                #self.graph.add((Literal(large_figure_url),FOAF['depicts'],(URIRef(cu.get_uri(morphology_term_id)))))


                #Do we want the label even if there is only one or the other definition?
                if subjective_definition != '' and objective_definition == '':
                    description = 'Subjective Description: '+subjective_definition
                elif objective_definition != '' and subjective_definition == '':
                    description = 'Objective Description: '+objective_definition
                elif subjective_definition != '' and objective_definition != '':
                    description = 'Objective Description: '+objective_definition+'; Subjective Description: '+subjective_definition
                else:
                    description = None

                if description is not None:
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),DC['description'],Literal(description)))

                #morphology_term_id has comment comments
                if comments != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),DC['comment'],Literal(comments)))


                #TODO: Need to fix handling of '; ' delimited entries.
                #morphology_term_id hasRelatedSynonym synonyms
                if synonyms != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),URIRef(cu.get_uri(self.relationship['hasRelatedSynonym'])),Literal(synonyms)))

                #TODO: Need to fix handling of '; ' delimited entries.
                if replaces != '' and replaces == synonyms:
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),URIRef(cu.get_uri(self.relationship['hasRelatedSynonym'])),Literal(replaces)))


                #Question: Add subject_definition and objective definition as a combined definition?
                #self.graph.add((URIRef(cu.get_uri(morphology_term_id)),RDF['type'],URIRef(cu.get_uri(self.terms['phenotype']))))






                if (limit is not None and line_counter > limit):
                    break
                #self.graph = eom.getGraph().__iadd__(self.graph)
        return

    def _process_eom_terms(self, raw, limit=None):
        '''
        This table contains the Elements of Morphology data that has been screen-scraped into DISCO.
        :param raw:
        :param limit:
        :return:
        '''

        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())
        line_counter = 0
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1

                (morphology_term_id, morphology_term_label, hp_id, hp_label,notes) = line.split('\t')



                if (limit is not None and line_counter > limit):
                    break

        return



    #TODO generalize this to a set of utils
    def _getcols(self,cur,table):
        query=(' ').join(("SELECT * FROM",table,"LIMIT 0"))  #for testing
        #print("COMMAND:",query)
        cur.execute(query)
        colnames = [desc[0] for desc in cur.description]
        print("COLS ("+table+"):",colnames)

        return


    def file_len(self,fname):
        with open(fname) as f:
            return sum(1 for line in f)




    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


