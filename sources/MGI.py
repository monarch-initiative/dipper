import csv
import os
from datetime import datetime
from stat import *
import re
import psycopg2


from rdflib import Literal
from rdflib.namespace import RDFS, OWL, RDF
from rdflib import Namespace, URIRef, BNode

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
        'all_allele_view'
#        'all_allele_mutation_view'
    ]


    relationship = {
        'is_mutant_of' : 'GENO:0000440',
        'derives_from' : 'RO:0001000',
        'has_alternate_part' : 'GENO:0000382',
        'has_reference_part' : 'GENO:0000385',
        'in_taxon' : 'RO:0000216',
        'has_zygosity' : 'GENO:0000608',   #what exactly "has zygosity"?  is it the allele?  genotype?
        'is_sequence_variant_instance_of' : 'GENO:0000408',
    }


    def __init__(self):
        Source.__init__(self, 'mgi')
        self.namespaces.update(curie_map.get())
        #assemble all the curie mappings from the imported models
        #self.namespaces.update(Assoc.curie_map)
        #self.namespaces.update(Genotype.curie_map)
        #self.namespaces.update(G2PAssoc.curie_map)

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
        #TODO any scrubbing needed for this resource?
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes: (none)
        :return: None
        '''

        #Should the wild type alleles be removed?

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





        #self._process_genotype_features(('/').join((self.rawdir,self.tables[1])), self.outfile, self.graph, limit)

        self._process_genotypes_new(('/').join((self.rawdir,'gxd_genotype_view')),limit)
        self._process_gxd_genotype_summary_view(('/').join((self.rawdir,'gxd_genotype_summary_view')),limit)
        self._process_all_summary_view(('/').join((self.rawdir,'all_summary_view')),limit)
        self._process_all_allele_view(('/').join((self.rawdir,'all_allele_view')),limit)
        #self._process_gxd_allele_pair_view(('/').join((self.rawdir,'gxd_allelepair_view')),limit)


        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Loaded", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):
        print("Processing Genotypes")
        #TODO




        line_counter = 0
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1
                cols = line.split('\t')
                genotype_key = cols[0]  # genotype key for connecting with alleles, first column.
                genotype_id = cols[10]  # mgiid, column 11.
                genotype_label = 'temporary label'
                #print("Capture", genotype_id, "and", genotype_key)
                geno = Genotype(genotype_id, genotype_label, self.namespaces)

        return

    def _process_genotypes_new(self, raw, limit=None):
        #need to make triples:
        #1.  genotype is a class  (or instance?)
        #2.  genotype has equivalentClass mgi internal identifier?  -- maybe make the internal identifier an anonymous node?
        #3.  genotype subclass of intrinsic_genotype
        #4.  genotype has_genomic_background strain_key
        #5.  strainkey has_label strain
        #6.  strainkey in_taxon taxon_id  #not part of this table
        #7.  strainkey is a class (or an instance?)

        has_reference_part = 'GENO:0000385'
        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())
        line_counter = 0
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1
                (genotype_key,strain_key,isconditional,note,existsas_key,createdby_key,modifiedby_key,creation_date,
                 modification_date,strain,mgiid,dbname,createdbymodifiedby,existsas,empty) = line.split('\t')

                #we can make these proper methods later
                gt = URIRef(cu.get_uri(mgiid))
                igt = BNode('genotypekey'+genotype_key)
                self.graph.add((gt,RDF['type'],Assoc.OWLCLASS))
                self.graph.add((gt,OWL['equivalentClass'],igt))
                self.graph.add((gt,Assoc.SUBCLASS,URIRef(cu.get_uri('GENO:0000000'))))
                istrain = BNode('strainkey'+strain_key)
                self.graph.add((istrain,RDF['type'],Assoc.OWLCLASS))
                self.graph.add((istrain,RDFS['label'],Literal(strain)))
                self.graph.add((gt,URIRef(cu.get_uri(has_reference_part)),istrain))
                #temporary assignment to Mus musculus
                self.graph.add((istrain,URIRef(cu.get_uri(self.relationship['in_taxon'])),URIRef(cu.get_uri('NCBITaxon:10090'))))

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_gxd_genotype_summary_view(self,raw,limit=None):
        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,createdby_key,modifiedby_key,
                 creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')
                #note the short_description is the GVC

                #we can make these proper methods later
                gt = URIRef(cu.get_uri(mgiid))
                igt = BNode('genotypekey'+object_key)
                self.graph.add((gt,RDF['type'],Assoc.OWLCLASS))
                self.graph.add((gt,OWL['equivalentClass'],igt))
                self.graph.add((gt,Assoc.SUBCLASS,URIRef(cu.get_uri('GENO:0000000'))))
                self.graph.add((gt,RDFS['label'],Literal(description)))  #the 'description' is the full genotype label

                if (limit is not None and line_counter > limit):
                    break

        return
    #NOTE: might be best to process alleles initially from the all_allele_view, as this does not have any repeats of alleles!
    def _process_all_summary_view(self,raw,limit):
        #Things to process:
        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                 createdby_key,modifiedby_key,creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')
                #NOTE:May want to filter alleles based on the preferred field (preferred = 1) or will get duplicates
                ## (24288, to be exact... Reduced to 480 if filtered on preferred = 1)
                #NOTE:Decision to make: use the short_description here or use the labels from the allelepair_view?
                #Use this table. allelepair_view will have more repititions of alleles due to inclusion in genotypes.

                #we can make these proper methods later
                #If we want to filter on preferred, use this if statement:
                if preferred == '1':
                    allele = URIRef(cu.get_uri(mgiid))
                    iallele = BNode('allelekey'+object_key)
                    self.graph.add((iallele,OWL['equivalentClass'],allele))
                    #self.graph.add(allele,RDFS['label'],Literal(short_description))
                    #This above statement isn't needed, as it would be redundant with the label provided for the allele Bnode, correct?
                    self.graph.add((allele,RDF['type'],Assoc.OWLCLASS))
                    #self.graph.add((al,OWL['equivalentClass'],ial)) #redundant, yes?
                    self.graph.add((allele,Assoc.SUBCLASS,URIRef(cu.get_uri('GENO:0000008'))))  # GENO:0000008 = allele
                if (limit is not None and line_counter > limit):
                    break

        return


    def _process_all_allele_view(self,raw,limit):
        #Need to make triples:
        variant_of = 'GENO:0000408' #FIXME:is_sequence_variant_instance_of. Is this correct?
        #GENO:0000440=is_mutant_of
        reference_of = 'GENO:0000409'#FIXME:is_reference_locus_instance_of. Is this correct?

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allele_key,marker_key,strain_key,mode_key,allele_type_key,allele_status_key,transmission_key,
                 collection_key,symbol,name,nomensymbol,iswildtype,isextinct,ismixed,createdby_key,modifiedby_key,
                 approvedby_key,approval_date,creation_date,modification_date,markersymbol,term,statusnum,strain,createby,modifiedby,approvedby) = line.split('\t')
                #NOTE:

                #we can make these proper methods later

                #Process the  alleles

                #al = URIRef(cu.get_uri(mgiid))
                iallele = BNode('allelekey'+allele_key)
                imarker = BNode('markerkey'+marker_key)
                if iswildtype == '0':
                    self.graph.add((iallele,URIRef(cu.get_uri(variant_of)),imarker))
                elif iswildtype == '1':
                    self.graph.add((iallele,URIRef(cu.get_uri(reference_of)),imarker))

                if re.match(".*<.*>.*", symbol):
                    #print(symbol)
                    symbol = re.sub(".*<", "<", symbol)
                    #print(symbol)
                    self.graph.add((iallele,RDFS['label'],Literal(symbol)))
                elif re.match("\+", symbol):
                    #TODO: Check to see if this is the proper handling, as while symbol is just +, marker symbol has entries without any <+>.
                    symbol = '<+>'
                    print(symbol)
                    self.graph.add((iallele,RDFS['label'],Literal(symbol)))
                else:
                    #print(symbol)
                    self.graph.add((iallele,RDFS['label'],Literal(symbol)))

                self.graph.add((iallele,OWL['OIO:hasExactSynonym'],Literal(symbol))) #FIXME: syntax correct for the OIO statement?)
                    #FIXME Need to handle alleles not in the *<*> format, such as many gene traps, induced mutations, and transgenics
                    #Handling by regex matching at the moment, but need to confirm that the format is dependent on subtype...
                if (limit is not None and line_counter > limit):
                    break

        return



    def _process_gxd_allele_pair_view(self,raw,limit):
        #Things to process:
        #allele_pair_key, genotype_key, allele1_key, allele2_key, allele_state
        #should already have symbol from the all_allele_view. marker_key?
        #Process the locus by combining allele1/allele2, with processing of those labels, or use the allele_keys to
        # link to the graph, grabbing the label there?
        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allelepair_key,genotype_key,allele_key_1,allele_key_2,marker_key,mutantcellline_key_1,mutantcellline_key_2,
                 pairstate_key,compound_key,sequencenum,createdby_key,modifiedby_key,creation_date,modification_date,symbol,
                 chromosome,allele1,allele2,allelestate,compound) = line.split('\t')
                #NOTE: symbol = gene/marker, allele1 + allele2 = VSLC, allele1/allele2 = variant locus, allelestate = zygosity
                #FIXME Need to handle alleles not in the *<*> format, such as many gene traps, induced mutations, and transgenics

                #we can make these proper methods later
                gt = URIRef(cu.get_uri(mgiid))
                igt = BNode('genotypekey'+object_key)
                self.graph.add((gt,RDF['type'],Assoc.OWLCLASS))
                self.graph.add((gt,OWL['equivalentClass'],igt))
                self.graph.add((gt,Assoc.SUBCLASS,URIRef(cu.get_uri('GENO:0000000'))))
                self.graph.add((gt,RDFS['label'],Literal(description)))  #the 'description' is the full genotype label

                if (limit is not None and line_counter > limit):
                    break

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
            'Disruption caused by insertion of vector': 'SO:0000667',  # insertion - correct?
            'Duplication': 'SO:1000035',  # duplication
            'Insertion': 'SO:0000667',  # insertion
            'Insertion of gene trap vector': 'SO:0000667',  # insertion - correct?
            'Intergenic deletion': 'SO:0000159',  # deletion
            'Intragenic deletion': 'SO:0000159',  # deletion
            'Inversion': 'SO:1000036',  # inversion
            'Not Applicable': 'SO:0001060',  # sequence variant - correct?
            'Not Specified': 'SO:0001060',  # sequence variant - correct?
            'Nucleotide repeat expansion': 'SO:0000667',  # insertion - correct?
            'Nucleotide substitutions': 'SO:1000002',  # substitution - Correct? Or another term indicating more than one?
            'Other': 'SO:0001060',  # sequence variant - correct?
            'Single point mutation': 'SO:1000008',  # point_mutation
            'Translocation': 'SO:0000199',  # translocation
            'Transposon insertion': 'SO:0000101',  # transposable_element
            'Undefined': 'SO:0001060',  # sequence variant - correct?
            'Viral insertion': 'SO:0000667',  # insertion - correct?
            'wild type': 'SO:0000817'  # wild type
        }
        if (allele_type.strip() in type_map):
            type = type_map.get(allele_type)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(allele_type)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Allele Type (", allele_type, ") not mapped")

        return type