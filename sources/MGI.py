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
        'gxd_genotype_view',
        'gxd_genotype_summary_view',
        'gxd_allelepair_view',
        'all_summary_view',
        'all_allele_mutation_view'
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
        #self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables,cxn)

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

        # What needs to be done? Above code is grabbing the tables from the MGI database
        # Do you need to assemble a table, or just grab the individual bits from various tables and match them up?
        # First test: grab the mgiid

        # Grab mgiid for genotype_id





        self._process_genotype_features(('/').join((self.rawdir,self.tables[1])), self.outfile, self.graph, limit)




        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Loaded", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):
        print("Processing Genotypes")
        #TODO
        # Genotype basics: Grab the genotype ID (mgiid) and create the genotype label
        # Will need to process the genotype label based on zygosity? No, MGI has allele pairs, just need to combine them,
        # and combine if there is more than one set of allele pairs for a genotype (more than one locus)
        # Grab data from the gdx_allelepair_view and process, allele1 and allele 2 to get the genotype label
        # Add in background? Assuming so, assuming we want the level of intrinsic/effective genotype.
        # Need to handle multiple loci (perhaps make an array, where entry 0 is the ID, append all additional entries, combine all additional entries with ';' when passing to genotype_label variable)

        # Update: Just creating a temporary genotype label for now. Will have to figure out how to get it later and in the right format.

        # NOTE: If this nested loop approach is the best route to go, see about adding breaks to speed up the processing


        #Is this the most efficient method? Requires looping through the allelepair data once for every genotype...
        #Requires > 4.8 billion allele loops so far...
        line_counter = 0
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1
                cols = line.split('\t')
                genotype_key = cols[0]  # genotype key for connecting with alleles, first column.
                genotype_id = cols[10]  # mgiid, column 11.
                #print("Capture", genotype_id, "and", genotype_key)

                # Map the genotype ID to the genotype label (temporary)
                table = (('/').join((self.rawdir,self.tables[2])))  # gxd_genotype_summary_view
                with open(table, 'r') as f2:
                    f2.readline()  # read the header row; skip
                    line_counter2 = 0
                    genotype_label = "Not specified"
                    for line in f2:
                        line_counter2 += 1
                        cols2 = line.split('\t')
                        #print(cols)
                        genotype_id_match = cols2[13]  # mgiid, column 14. Temporarily used for genotype label matching.
                        #print("Genotype ID", genotype_id, "GenotypeID Match", genotype_id_match)
                        genotype_label_match = cols2[16]



                        if len(genotype_label_match.strip()) != 0 and genotype_id == genotype_id_match.strip() :
                            #genotype_label = cols2[16] #short_description is column 17, which will work for a temporary genotype label.
                            genotype_label = genotype_label_match

                            #print("Genotype ID", genotype_label)
                            genotype_label = genotype_label.strip()
                            #print("Captured", genotype_id, "for genotype", genotype_label)


                    #print("Captured", genotype_id, "for genotype", genotype_label)
                    geno = Genotype(genotype_id, genotype_label, self.namespaces)
                    #self.graph.__iadd__(geno.getGraph())

                f2.close()


                #Add alleles (aka sequence alterations not variant locus, although they are the same for this resource)
                #Use the allelepair_view to bring in allele data for each genotype.
                #Can grab the two alleles for a locus, and the zygosity label.
                #However, this table does not have the official MGI ID used in the view in DISCO.
                #sequence_alteration_id/variant_locus_id (which are the same) are obtained from all_summary_view table,
                #which can be used for the allele id.

                # Could create an allele object first. Can do this from the all_summary_view
                # Has mgiid (column 14, for sequence_alteration_id/variant_locus_id),
                # a label (short_description, column 17),
                # and an allele subtype (column 15) which could be converted to an allele type. (NOT IDEAL)

                # Need to get from the _genotype_key -> _allele_key
                # Get mgiid from the all_summary_view as allele_id
                # Get short_description from all_summary_view as allele_name, stripping out the text before the <XYZ>.
                # Question: Do we keep the wild type alleles?
                # Get mutation from the all_allele_mutation_view as allele_type, map to GO/SO terms

                # Connect _genotype_key to _genotype_key in gxd_allele_pair_view
                # Capture allele1_key and allele2_key, use to match individually to all_summary_view (_object_key) and all_allele_mutation_view (_allele_key)
                # Connect gxd_allelepair_view to all_summary_view on _allele_key = _object_key for allele_id and allele_name
                # Connect gxd_allelepair_view to all_allele_mutation_view on _allele_key = _allele_key for allele_type
                # Map allele_type
                #TODO Zygosity label.




                table = (('/').join((self.rawdir,self.tables[3]))) #gxd_allelepair_view
                with open(table, 'r') as f2:
                    f2.readline()  # read the header row; skip
                    for line in f2:
                        allele1_key = ''
                        allele2_key = ''
                        allele1_id = ''
                        allele2_id = ''
                        allele1_name = ''
                        allele2_name = ''
                        allele1_type = ''
                        allele2_type = ''
                        cols2 = line.split('\t')
                        genotype_key_match = cols2[1]  # _genotype_key, column 2
                        if genotype_key == genotype_key_match:
                            allele1_key = cols2[2]  # _allele_key_1, column 3
                            allele2_key = cols2[3]  # _allele_key_2, column 4

                            #Note: Filter out the nulls for allele 2 at some point

                            #Next, match each allele_key to the all_summary_view table for allele_name and allele_id
                            table = (('/').join((self.rawdir,self.tables[4]))) #all_summary_view
                            with open(table, 'r') as f3:
                                f3.readline()  # read the header row; skip
                                for line in f3:
                                    cols3 = line.split('\t')
                                    allele_match = cols3[5]  # _object_key, column 6
                                    if allele1_key == allele_match:
                                        allele1_id = cols3[13]  # mgiid, column 14
                                        allele1_name = cols3[16]  # short_description, column 17
                                        allele1_name = allele1_name.strip()
                                    if allele2_key == allele_match:
                                        allele2_id = cols3[13]  # mgiid, column 14
                                        allele2_name = cols3[16]  # short_description, column 17
                                        allele2_name = allele2_name.strip()
                            f3.close()




                            #Next, match each allele_key to the all_allele_mutation_view for allele_type
                            table = (('/').join((self.rawdir,self.tables[5]))) #all_allele_mutation_view
                            with open(table, 'r') as f4:
                                f4.readline()  # read the header row; skip
                                for line in f4:
                                    cols4 = line.split('\t')
                                    allele_match = cols4[0]
                                    if allele1_key == allele_match:
                                        allele1_type = cols4[4]
                                    if allele2_key == allele_match:
                                        allele2_type = cols4[4]
                            f4.close()



                            #Map the allele types to GENO/SO terms and add to genotype
                            #How should we handle wild type alleles?
                            #print("Allele_1_Type", allele1_type, "Allele_2_Type", allele2_type)
                            allele1_type = self._map_allele_type_to_geno(allele1_type)
                            geno.addAllele(allele1_id, allele1_name, allele1_type)
                            geno.addAlleleToGenotype(genotype_id, allele1_id)




                            if allele2_id is not None and len(allele2_id.strip()) != 0:
                                #print('allele2_id is', allele2_id, 'allele2_name is', allele2_name, 'allele2_type is ', allele2_type)
                                if allele2_type.strip() != '':
                                    allele2_type = self._map_allele_type_to_geno(allele2_type)
                                #FIXME leaving allele type empty is throwing an error: TypeError: argument of type 'NoneType' is not iterable
                                #Setting empty types to wild type for the moment.
                                else:
                                    allele2_type = 'wild type'
                                    allele2_type = self._map_allele_type_to_geno(allele2_type)
                                #print('allele2_id is', allele2_id, 'allele2_name is', allele2_name, 'allele2_type is ', allele2_type)
                                geno.addAllele(allele2_id, allele2_name, allele2_type)
                                geno.addAlleleToGenotype(genotype_id, allele2_id)


                            #Add each allele to the genotype object



                            #if len(allele2_id.strip()) != 0:



                        #allele_id = cols2[13]  # mgiid for the allele, 14th column
                        #allele_label = cols2[16]  # short_description, 17th column
                        #allele_type = cols2[14]  # subtype, 15th column, not ideal but will use temporarily.
                        # Proper allele type needs to be extracted from all_allele_mutation_view,
                        # match using the _allele_key

                    #Create the add the alleles to the genotype


                f2.close()
                self.graph.__iadd__(geno.getGraph())




        f1.close()






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


    #TODO: Finish identifying SO/GENO terms for mappings for those found in MGI
    def _map_allele_type_to_geno(self, allele_type):
        type = None
        type_map = {
            'Deletion': 'SO:0000159',  # deletion
            'Disruption caused by insertion of vector': 'SO:XXXXXXX',
            'Duplication': 'SO:XXXXXXX',
            'Insertion': 'SO:XXXXXXX',
            'Insertion of gene trap vector': 'SO:XXXXXXX',
            'Intergenic deletion': 'SO:XXXXXXX',
            'Intragenic deletion': 'SO:XXXXXXX',
            'Inversion': 'SO:XXXXXXX',
            'Not Applicable': 'SO:XXXXXXX',
            'Not Specified': 'SO:XXXXXXX',
            'Nucleotide repeat expansion': 'SO:XXXXXXX',
            'Nucleotide substitutions': 'SO:XXXXXXX',
            'Other': 'SO:XXXXXXX',
            'Single point mutation': 'SO:1000008',  #point_mutation
            'Translocation': 'SO:0000199',  #translocation
            'Transposon insertion': 'SO:XXXXXXX',
            'Undefined': 'SO:XXXXXXX',
            'Viral insertion': 'SO:XXXXXXX',
            'wild type': 'SO:XXXXXXX'




            #'complex_substitution': 'SO:1000005',  # complex substitution
            #'deficiency': 'SO:1000029',  # incomplete chromosome
            #'deletion': 'SO:0000159',  # deletion
            #'indel': 'SO:1000032',  #indel
            #'insertion': 'SO:0000667',  #insertion
            #'point_mutation': 'SO:1000008',  #point_mutation
            #'sequence_variant': 'SO:0001060',  #sequence variant
            #'transgenic_insertion': 'SO:0001218',  #transgenic insertion
            #'transgenic_unspecified': 'SO:0000781',  #transgenic unspecified
            #'transloc': 'SO:0000199',  #translocation
            #            'unspecified' : None

        }
        if (allele_type.strip() in type_map):
            type = type_map.get(allele_type)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(allele_type)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Allele Type (", allele_type, ") not mapped")

        return type