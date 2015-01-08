import csv
import os
from datetime import datetime
from stat import *
import re
import psycopg2


from rdflib import Literal
from rdflib.namespace import RDFS, OWL, RDF
from rdflib import Namespace, URIRef

from utils import pysed
from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from utils.CurieUtil import CurieUtil
import config


class MGI(Source):
    '''
    Be sure to have pg user/password connection details in your conf.json file, like:
      dbauth : {
        'mgi' : {'user' : '<username>', 'password' : '<password>'}
      }
    '''
# tables in existing interop
# mgi_organism_acc_view mgi_organism_view gxd_genotype_view gxd_allelepair_view mrk_marker_view
# mgi_reference_allele_view bib_acc_view voc_annot_view gxd_genotype_summary_view voc_evidence_view
# all_allele_cellline_view voc_term_view all_allele_view all_allele_mutation_view prb_strain_view
# all_summary_view mrk_summary_view mgi_note_vocevidence_view acc_logicaldb_view mgi_note_strain_view
# prb_strain_acc_view prb_strain_summary_view prb_strain_marker_view
    tables = [
        'mgi_dbinfo',
        'gxd_genotype_view'
    ]

    namespaces = {
        'MP': 'http://purl.obolibrary.org/obo/MP_',
        'MA': 'http://purl.obolibrary.org/obo/MA_',
        'MGI' : 'http://www.informatics.jax.org/accession/MGI:'  #note that MGI genotypes will not resolve at this address
    }

    def __init__(self):
        Source.__init__(self, 'mgi')

        #assemble all the curie mappings from the imported models
        self.namespaces.update(Assoc.curie_map)
        self.namespaces.update(Genotype.curie_map)
        self.namespaces.update(G2PAssoc.curie_map)

        #update the dataset object with details about this resource
        self.dataset = Dataset('mgi', 'MGI', 'http://www.informatics.jax.org/')

        #check if config exists; if it doesn't, error out and let user know
        if (not (('dbauth' in config.get_config()) and ('mgi' in config.get_config()['dbauth']))):
            print("ERROR: not configured with PG user/password.")
        return

        #source-specific warnings.  will be cleared when resolved.
        #print("WARN: we are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self):
        '''
        For the MGI resource, we connect to the remote database, and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        '''

        #create the connection details for MGI
        cxn = config.get_config()['dbauth']['mgi']
        cxn.update({'host' : 'adhoc.informatics.jax.org', 'database' : 'mgd', 'port' : 5432 })

        self.dataset.setFileAccessUrl(('').join(('jdbc:postgresql://',cxn['host'],':',str(cxn['port']),'/',cxn['database'])))

        #process the tables
        self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        #self.fetch_from_pgdb(self.tables,cxn)

        datestamp=ver=None
        #get the resource version information from table mgi_dbinfo, already fetched above
        outfile=('/').join((self.rawdir,'mgi_dbinfo'))
        if os.path.exists(outfile):
            st = os.stat(outfile)
            with open(outfile, 'r') as f:
                f.readline() #read the header row; skip
                info = f.readline()
                cols = info.split('\t')
                ver = cols[0] #col 0 is public_version
                ver = ver.replace('MGI ','')  #MGI 5.20 --> 5.20
                #MGI has a datestamp for the data within the database; use it instead of the download date
                #datestamp in the table: 2014-12-23 00:14:20
                d = cols[7].strip()  #modification date
                datestamp = datetime.strptime(d, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
                f.close()

        self.dataset.setVersion(datestamp,ver)

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes: (none)
        :return: None
        '''

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")



        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)


        ##### Write it out #####

        filewriter = open(self.outfile, 'w')
        print(self.graph.serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        filewriter = open(self.datasetfile, 'w')
        print(self.dataset.getGraph().serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        print("Wrote", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):
        print("Processing Genotypes")
        #TODO
        return



    def _process_g2p(self, raw, out, g, limit=None):
        '''
        This module currently filters for only wild-type environments, which clearly excludes application
        of morpholinos.  Very stringent filter.  To be updated at a later time.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        '''
        print("Processing G2P")
        line_counter = 0
        # hardcode
        eco_id = "ECO:0000059"  #experimental_phenotypic_evidence

        #TODO


        return


    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


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
            for i, l in enumerate(f):
                pass
        return i + 1