import csv
from utils import pysed
import os, datetime
from datetime import datetime
from stat import *
from rdflib import Graph, Literal
from rdflib.namespace import RDFS, OWL, RDF
import gzip,os.path

from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from rdflib import Namespace, URIRef
import re
from utils.CurieUtil import CurieUtil
import curie_map


class IMPC(Source):

    files = {
        'impc' : {'file' : 'IMPC_genotype_phenotype.csv.gz',
                  'url' : 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/IMPC_genotype_phenotype.csv.gz'},
        'euro' : {'file' : 'EuroPhenome_genotype_phenotype.csv.gz',
                  'url' : 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/EuroPhenome_genotype_phenotype.csv.gz'},
        'mgd' : {'file' : 'MGP_genotype_phenotype.csv.gz',
                 'url' : 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/MGP_genotype_phenotype.csv.gz'}
    }


    def __init__(self):
        Source.__init__(self, 'impc')

        #update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('impc', 'IMPC', 'http://www.mousephenotype.org')

        #source-specific warnings.  will be cleared when resolved.
        #print("WARN: we are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self, is_dl_forced):
        #this is fetching the standard files, not from the API/REST service
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                ('/').join((self.rawdir,file['file'])),
                                is_dl_forced)
            self.dataset.setFileAccessUrl(file['url'])
            st = os.stat(('/').join((self.rawdir,file['file'])))

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
        #pysed.replace("\\\\", '', ('/').join((self.rawdir,self.files['geno']['file'])))

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")
        # IMPC data is delivered in three separate csv files. Need to iterate over each process (genotype, G2P, etc.)
        # for each file, unless instead we merge all three files into one, then process once.


        self._process_genotype_features(('/').join((self.rawdir,self.files['impc']['file'])), self.outfile, self.graph, limit)
        self._process_g2p(('/').join((self.rawdir,self.files['impc']['file'])), self.outfile, self.graph, limit)


        #self._process_genotype_features(('/').join((self.rawdir,self.files['geno']['file'])), self.outfile, self.graph, limit)
        #self._process_g2p(('/').join((self.rawdir,self.files['pheno']['file'])), self.outfile, self.graph, limit)
        #self._process_pubinfo(('/').join((self.rawdir,self.files['pubs']['file'])), self.outfile, self.graph, limit)

        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Found", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):


        print("Processing Genotypes")

        line_counter = 0
        with gzip.open(raw, 'rt') as csvfile:
        #with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                line_counter += 1
                #print(row)

                (resource_name,phenotyping_center,colony_id,strain_name,strain_accession_id,marker_symbol,
                 marker_accession_id,allele_symbol,allele_accession_id,zygosity,sex,pipeline_name,pipeline_stable_id,
                 procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,mp_term_name,mp_term_id,p_value,
                 effect_size) = row

                # For IMPC, we put together genotype IDs from other parts (marker, allele, and strain IDs, with
                # adjustments based on zygosity).



                # Making genotype IDs from the various parts, can change later if desired.
                # Should these just be the genotype IDs, or should we use the self._makeInternalIdentifier()?
                # Or perhaps this, adjusted base on zygosity:
                # genotype_id = self.make_id((marker_accession_id+allele_accession_id+strain_accession_id))
                # OR, even more simplified:
                # genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id))
                # If we use this even more simplified method for the genotype ID,
                # may want to move the zygosity warning to the zygosity processing portion.
                genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id))

                # Make the variant locus name/label
                if re.match(".*<.*>.*", allele_symbol):
                    variant_locus_name = allele_symbol
                else:
                    variant_locus_name = allele_symbol+'<'+allele_symbol+'>'

                # Making VSLC labels from the various parts, can change later if desired.
                if zygosity == 'heterozygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<+>',variant_locus_name)
                    #print(vslc_name)
                elif zygosity == 'homozygote':
                    vslc_name = variant_locus_name+'/'+variant_locus_name
                    #print(vslc_name)
                elif zygosity == 'hemizygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<0>',variant_locus_name)
                    #print(vslc_name)
                elif zygosity == 'not_applicable':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<?>',variant_locus_name)
                    #print(vslc_name)
                # Do we need to handle an unknown zygosity label in another fashion?
                # How do we want to handle the "not_applicable" zygosity labels?
                # Is the question mark allele the correct way?
                else:
                    print("INFO: found unknown zygosity :",zygosity)
                    break

                # Create the GVC labels
                # For IMPC, the GVC == VSLC
                gvc_name = vslc_name


                # Making genotype labels from the various parts, can change later if desired.
                genotype_name = vslc_name+'['+strain_name+']'


                # Create the genotype
                geno = Genotype(genotype_id, genotype_name)
                #print('ID='+genotype_id+'     LABEL='+genotype_name)

                #FIXME
                #IMPC contains targeted mutations with either gene traps, knockouts, insertion/intragenic deletions.
                #Currently hard-coding to the insertion type. Does this need to be adjusted?
                allele_type = 'SO:0000667'


                #This is for handling any of the alleles that do not have an MGI ID or IMPC has not yet
                # updated their data for the allele with the MGI ID. The IDs look like NULL-<10-digit string>.
                if re.match("MGI:.*",allele_accession_id):
                    allele_id = allele_accession_id
                else:
                    allele_id = 'IMPC:'+allele_accession_id

                #Add allele to genotype
                geno.addAllele(allele_id, allele_symbol, allele_type)

                #Hard coding gene_type as gene.
                gene_type = 'SO:0000704'# gene

                if re.match('MGI:.*',marker_accession_id):
                    gene_id = marker_accession_id
                else:
                    gene_id = 'IMPC:'+marker_accession_id


                #Add gene to genotype
                #need gene_id, gene_label, gene_type. No gene description in data.
                #All marker IDs have the MGI: prefix, but should we include an error message should a
                # future data set not be in that format?
                geno.addGene(gene_id, marker_symbol,gene_type)

                # Add allele to gene
                geno.addAlleleOfGene(allele_id, gene_id)

                # genotype has_part allele
                geno.addAlleleToGenotype(genotype_id, allele_id)

                # need to make some attributes of this relationship for allele zygosity within the genotype
                #or do we just make the allele_complement here, based on the zygosity?
                #allele_in_gene_id=self.make_id(genotype_id+allele_id+zygosity)
                #allele has_disposition zygosity?
                #                g.add(())

                #TODO
                #The following parts are test code for creating the more complex parts of the genotype partonomy.
                #Will work on abstracting these code snippets to generalized classes for use by any resource.

                # Create the VSLC ID

                vslc_id = self.make_id((marker_accession_id+allele_accession_id+zygosity))
                #Alternatively could use an if statement based on zygosity like for the label, but this is simpler.




                # Link the VSLC to the allele

                # Add the zygosity to the VSLC

                # Create the GVC

                # Add the GVC to the intrinsic genotype

                # Link the GVC to the VSLC

                #TODO: Create the effective genotype label/id by adding the sex of the mouse.



                #add the specific genotype subgraph to the overall graph
                if (limit is not None and line_counter > limit):
                    break
                self.graph = geno.getGraph().__iadd__(self.graph)
        print("INFO: Done with genotypes")
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

        #NOTE: No evidence provided in the data file, no environment data

        #QUESTION: A matter of efficiency, but is it better to process the genotype ang G2P sections separately,
        # or in a case like IMPC with single files, we could potentially process the G2P portions within the Genotype call.
        # Current implementation results in processing the files multiple times, which isn't terrible given that
        # the files for IMPC are not large files.


        print("Processing G2P")


        line_counter = 0
        with gzip.open(raw, 'rt') as csvfile:
        #with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                line_counter += 1
                #print(row)

                (resource_name,phenotyping_center,colony_id,strain_name,strain_accession_id,marker_symbol,
                 marker_accession_id,allele_symbol,allele_accession_id,zygosity,sex,pipeline_name,pipeline_stable_id,
                 procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,mp_term_name,mp_term_id,p_value,
                 effect_size) = row


                phenotype_id = mp_term_id

                #Not currently using the phenotype name, but add it somewhere as the label?
                phenotype_name = mp_term_name

                # Making genotype IDs from the various parts, can change later if desired.
                # Should these just be the genotype IDs, or should we use the self._makeInternalIdentifier()?
                genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id))

                # Make the variant locus name/label
                if re.match(".*<.*>.*", allele_symbol):
                    variant_locus_name = allele_symbol
                else:
                    variant_locus_name = allele_symbol+'<'+allele_symbol+'>'

                # Making VSLC labels from the various parts, can change later if desired.
                if zygosity == 'heterozygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<+>',variant_locus_name)
                elif zygosity == 'homozygote':
                    vslc_name = variant_locus_name+'/'+variant_locus_name
                elif zygosity == 'hemizygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<0>',variant_locus_name)
                elif zygosity == 'not_applicable':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<?>',variant_locus_name)
                else:
                    print("INFO: found unknown zygosity :",zygosity)
                    break

                # Create the GVC labels
                # For IMPC, the GVC == VSLC
                gvc_name = vslc_name

                # Making genotype labels from the various parts, can change later if desired.
                genotype_name = gvc_name+'['+strain_name+']'

                geno = Genotype(genotype_id, genotype_name)
                self.graph.__iadd__(geno.getGraph())

                #FIXME
                #In the NIF/DISCO view, we don't really have a publication id. We have a publication URL.
                #pub_url example: http://www.mousephenotype.org/data/genes/MGI:108077
                #There are duplicates, and will be duplicates of the genotype/phenotype/pub_id combination.
                #but could add the testing parameter?
                pub_id = 'IMPC:'+marker_accession_id

                #FIXME
                #No evidence code provided. NIF/DISCO view was hard coded to null. However,
                # isn't all of the IMPC data based on experimental evidence? Or is that too general?
                # Could use 'EXP': 'ECO:0000006', although the code below used in ZFIN might be more appropriate.
                eco_id = "ECO:0000059"  #experimental_phenotypic_evidence This was used in ZFIN

                assoc_id = self.make_id((genotype_id+phenotype_id+pub_id))

                assoc = G2PAssoc(assoc_id, genotype_id, phenotype_id, pub_id, eco_id)
                self.graph = assoc.addAssociationNodeToGraph(self.graph)




                if (limit is not None and line_counter > limit):
                    break



        return



    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


