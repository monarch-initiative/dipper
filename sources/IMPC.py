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



        #self._process_genotype_features(('/').join((self.rawdir,self.files['geno']['file'])), self.outfile, self.graph, limit)
        #self._process_g2p(('/').join((self.rawdir,self.files['pheno']['file'])), self.outfile, self.graph, limit)
        #self._process_pubinfo(('/').join((self.rawdir,self.files['pubs']['file'])), self.outfile, self.graph, limit)

        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Found", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):

        #TODO make this more efficient
        #the problem with this implementation is that it creates many genotypes over and over, if the
        #same genotype has many features (on many rows) then the same genotype is recreated, then it must be
        #merged.  We should probably just create the genotype once, and then find the other
        #items that belong to that genotype.
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
                 marker_accesssion_id,allele_symbol,allele_accession_id,zygosity,sex,pipeline_name,pipeline_stable_id,
                 procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,mp_term_name,mp_term_id,p_value,
                 effect_size) = row

                # For IMPC, we put together genotype IDs from other parts (marker, allele, and strain IDs, with
                # adjustments based on zygosity).



                # Making genotype IDs from the various parts, can change later if desired.
                # Should these just be the genotype IDs, or should we use the self._makeInternalIdentifier()?
                if zygosity == 'heterozygote':
                    genotype_id = 'IMPC:'+marker_accesssion_id.replace(':','_')+'_'+allele_accession_id.replace(':','_')+'_W_'+strain_accession_id.replace(':','_')
                    #print(genotype_id)
                elif zygosity == 'homozygote':
                    genotype_id = 'IMPC:'+marker_accesssion_id.replace(':','_')+'_'+allele_accession_id.replace(':','_')+'_'+marker_accesssion_id.replace(':','_')+'_'+allele_accession_id.replace(':','_')+'_'+strain_accession_id.replace(':','_')
                    #print(genotype_id)
                elif zygosity == 'hemizygote':
                    genotype_id = 'IMPC:'+marker_accesssion_id.replace(':','_')+'_'+allele_accession_id.replace(':','_')+'_0_'+strain_accession_id.replace(':','_')
                    #print(genotype_id)
                    #Do something completely different
                elif zygosity == 'not_applicable':
                    genotype_id = 'IMPC:'+marker_accesssion_id.replace(':','_')+'_'+allele_accession_id.replace(':','_')+'_?_'+strain_accession_id.replace(':','_')
                    #print(genotype_id)

                # Do we need to handle an unknown zygosity label in another fashion?
                # How do we want to handle the "not_applicable" zygosity labels?
                # Is the question mark allele the correct way?
                else:
                    print("INFO: found unknown zygosity :",zygosity)



                # Make the variant locus name/label
                if re.match(".*<.*>.*", allele_symbol):
                    variant_locus_name = allele_symbol
                else:
                    variant_locus_name = allele_symbol+'<'+allele_symbol+'>'




                # Making genotype labels from the various parts, can change later if desired.
                if zygosity == 'heterozygote':
                    genotype_name = variant_locus_name+'/'+re.sub('<.*','<+>',variant_locus_name)+'['+strain_name+']'
                    #print(genotype_name)
                elif zygosity == 'homozygote':
                    genotype_name = variant_locus_name+'/'+variant_locus_name+'['+strain_name+']'
                    #print(genotype_name)
                elif zygosity == 'hemizygote':
                    genotype_name = variant_locus_name+'/'+re.sub('<.*','<0>',variant_locus_name)+'['+strain_name+']'
                    #print(genotype_name)
                    #Do something completely different
                elif zygosity == 'not_applicable':
                    genotype_name = variant_locus_name+'/'+re.sub('<.*','<?>',variant_locus_name)+'['+strain_name+']'
                    print(genotype_name)

                # Do we need to handle an unknown zygosity label in another fashion?
                # How do we want to handle the "not_applicable" zygosity labels?
                # Is the question mark allele the correct way?
                else:
                    print("INFO: found unknown zygosity :",zygosity)




                '''
                genotype_id = 'ZFIN:' + genotype_id.strip()
                geno = Genotype(genotype_id, genotype_name)

                # reassign the allele_type to a proper GENO or SO class
                allele_type = self._map_allele_type_to_geno(allele_type)

                allele_id = 'ZFIN:' + allele_id.strip()
                geno.addAllele(allele_id, allele_name, allele_type)

                if (gene_id is not None and gene_id.strip() != ''):
                    gene_id = 'ZFIN:' + gene_id.strip()
                    geno.addGene(gene_id, gene_symbol)

                    # if it's a transgenic construct, then we'll have to get the other bits
                    if (construct_id is not None and construct_id.strip() != ''):
                        construct_id = 'ZFIN:' + construct_id.strip()
                        geno.addAlleleDerivesFromConstruct(allele_id, construct_id)

                    # allele to gene
                    geno.addAlleleOfGene(allele_id, gene_id)


                # genotype has_part allele
                geno.addAlleleToGenotype(genotype_id, allele_id)
                # need to make some attributes of this relationship for allele zygosity within the genotype
                #or do we just make the allele_complement here, based on the zygosity?
                #allele_in_gene_id=self.make_id(genotype_id+allele_id+zygosity)
                #allele has_disposition zygosity?
                #                g.add(())

                if (limit is not None and line_counter > limit):
                    break

                #add the specific genotype subgraph to the overall graph
                self.graph = geno.getGraph().__iadd__(self.graph)

                '''

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
        print("Processing G2P")




        return



    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


