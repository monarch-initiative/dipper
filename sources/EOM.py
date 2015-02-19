import os
from datetime import datetime
from stat import *

from sources.Source import Source

from models.Assoc import Assoc
from models.Dataset import Dataset
from utils.CurieUtil import CurieUtil
from conf import config, curie_map
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



                if (limit is not None and line_counter > limit):
                    break

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


