import os
from datetime import datetime
from stat import *
import re
import logging


from rdflib import Literal
from rdflib.namespace import DC, FOAF
from rdflib import URIRef

from sources.Source import Source
from models.Assoc import Assoc
from models.Dataset import Dataset
from utils.CurieUtil import CurieUtil
from conf import config, curie_map
from utils.GraphUtils import GraphUtils

logger = logging.getLogger(__name__)

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
        'hasRelatedSynonym': 'OIO:hasRelatedSynonym',
    }


    def __init__(self):
        Source.__init__(self, 'eom')
        self.namespaces.update(curie_map.get())

        #update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('eom', 'EOM', 'http://elementsofmorphology.nih.gov')

        #check if config exists; if it doesn't, error out and let user know
        if (not (('dbauth' in config.get_config()) and ('disco' in config.get_config()['dbauth']))):
            logger.error("ERROR: not configured with PG user/password.")

        #source-specific warnings.  will be cleared when resolved.

        return


    def fetch(self, is_dl_forced):

        #create the connection details for DISCO
        cxn = config.get_config()['dbauth']['disco']
        cxn.update({'host' : 'nif-db.crbs.ucsd.edu', 'database' : 'disco_crawler', 'port' : 5432 })

        self.dataset.setFileAccessUrl(('').join(('jdbc:postgresql://',cxn['host'],':',str(cxn['port']),'/',cxn['database'])))

        #process the tables
        #self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables,cxn)


        #FIXME: Everything needed for data provenance?
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
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")

        self._process_nlx_157874_1_view(('/').join((self.rawdir,'dv.nlx_157874_1')),limit)
        self._map_eom_terms(('/').join((self.rawdir,'eom_terms.tsv')),limit)

        logger.info("Finished parsing.")


        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        logger.info("Found %s nodes", len(self.graph))
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

                #Assemble the description text
                if subjective_definition != '' and objective_definition == '':
                    description = 'Subjective Description: '+subjective_definition
                elif objective_definition != '' and subjective_definition == '':
                    description = 'Objective Description: '+objective_definition
                elif subjective_definition != '' and objective_definition != '':
                    description = 'Objective Description: '+objective_definition+'; Subjective Description: '+subjective_definition
                else:
                    description = None

                #Add morphology term to graph as a class with label, type, and description.
                gu.addClassToGraph(self.graph,morphology_term_id,morphology_term_label,self.terms['phenotype'],description)

                #morphology_term_id has depiction small_figure_url
                if small_figure_url != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),FOAF['depiction'],Literal(small_figure_url)))

                #morphology_term_id has depiction large_figure_url
                if large_figure_url != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),FOAF['depiction'],Literal(large_figure_url)))

                #morphology_term_id has comment comments
                if comments != '':
                    self.graph.add((URIRef(cu.get_uri(morphology_term_id)),DC['comment'],Literal(comments)))

                #morphology_term_id hasRelatedSynonym synonyms (; delimited)
                if synonyms != '':
                    items = synonyms.split(';')
                    for i in items:
                        self.graph.add((URIRef(cu.get_uri(morphology_term_id)),URIRef(cu.get_uri(self.relationship['hasRelatedSynonym'])),Literal(i)))

                #morphology_term_id hasRelatedSynonym replaces (; delimited)
                if replaces != '' and replaces != synonyms:
                    items = replaces.split(';')
                    for i in items:
                        self.graph.add((URIRef(cu.get_uri(morphology_term_id)),URIRef(cu.get_uri(self.relationship['hasRelatedSynonym'])),Literal(i)))

                #FIXME: What is the proper predicate for a web page link? morphology_term_id hasURL morphology_term_url?
                #morphology_term_id has page morphology_term_url
                self.graph.add((URIRef(cu.get_uri(morphology_term_id)),FOAF['page'],Literal(morphology_term_url)))


                if (limit is not None and line_counter > limit):
                    break
                #self.graph = eom.getGraph().__iadd__(self.graph)
        return

    def _map_eom_terms(self, raw, limit=None):
        '''
        This table contains the HP ID mappings from the local tsv file.
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

                #Sub out the underscores for colons.
                hp_id = re.sub('_', ':', hp_id)
                if re.match(".*HP:.*", hp_id):
                    #Add the HP ID as an equivalent class
                    gu.addEquivalentClass(self.graph,morphology_term_id,hp_id)
                else:
                    logger.warning('No matching HP term for %s',morphology_term_label)

                if (limit is not None and line_counter > limit):
                    break

        return





    #TODO generalize this to a set of utils
    def _getcols(self,cur,table):
        query=(' ').join(("SELECT * FROM",table,"LIMIT 0"))  #for testing
        #print("COMMAND:",query)
        cur.execute(query)
        colnames = [desc[0] for desc in cur.description]
        logger.info("COLS (%s): %s", table, colnames)

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


