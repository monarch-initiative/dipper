import csv
import os
from datetime import datetime
from stat import *
import re

from utils import pysed
from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
from conf import curie_map


class ZFIN(Source):

    files = {
        'geno' : {'file' : 'genotype_features.txt', 'url' : 'http://zfin.org/downloads/genotype_features.txt'},
        'pheno' : {'file' : 'phenotype.txt', 'url' : 'http://zfin.org/downloads/phenotype.txt'},
        'pubs' : {'file' : 'zfinpubs.txt', 'url' : 'http://zfin.org/downloads/zfinpubs.txt'},
        'zpmap' : {'file' : 'zp-mapping.txt', 'url' : 'https://phenotype-ontologies.googlecode.com/svn/trunk/src/ontology/zp/zp-mapping.txt'},
    }


    def __init__(self):
        Source.__init__(self, 'zfin')


        #update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('zfin', 'ZFIN', 'http://www.zfin.org')

        #source-specific warnings.  will be cleared when resolved.
        print("WARN: we are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self, is_dl_forced):

        #fetch all the files
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                ('/').join((self.rawdir,file['file'])),
                                is_dl_forced)
            self.dataset.setFileAccessUrl(file['url'])
            # zfin versions are set by the date of download.
            st = os.stat(('/').join((self.rawdir,file['file'])))
        self.scrub()

        #this will set the version based on the last-ingested file.
        #TODO should be a date-stamp for each file?  how to track that prov?
        self.dataset.setVersion(datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d"))

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:
        * remove oddities where there are "\" instead of empty strings
        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        pysed.replace("\\\\", '', ('/').join((self.rawdir,self.files['geno']['file'])))

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")

        self._load_zp_mappings()

        self._process_genotype_features(limit)
        self._process_g2p(('/').join((self.rawdir,self.files['pheno']['file'])), self.outfile, self.graph, limit)
        self._process_pubinfo(('/').join((self.rawdir,self.files['pubs']['file'])), self.outfile, self.graph, limit)

        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Found", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, limit=None):
        """
        We don't actually know the allele pairs, so we'll store some info in a hashmap for post-processing
        :param limit:
        :return:
        """
        raw = ('/').join((self.rawdir,self.files['geno']['file']))
        out = self.outfile

        geno_hash = {}
        print("Processing Genotypes")
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (genotype_id, genotype_name, genotype_unique_name, allele_id, allele_name, allele_ab,
                 allele_type, allele_disp_type, gene_symbol, gene_id, zygosity,
                 construct_name, construct_id, other) = row

                genotype_id = 'ZFIN:' + genotype_id.strip()
                geno = Genotype(self.graph)
                gt = geno.addGenotype(genotype_id, genotype_name)
                if genotype_id not in geno_hash:
                    geno_hash[genotype_id] = {};
                genoparts = geno_hash[genotype_id]

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
                    if gene_id not in genoparts:
                        genoparts[gene_id] = [allele_id]
                    else:
                        genoparts[gene_id].append(allele_id)

                    if (zygosity == 'homozygous'):
                        genoparts[gene_id].append(allele_id)  #add the allele again
                    elif (zygosity == 'unknown'):
                        genoparts[gene_id].append('?')  #we'll just use this as a convention for unknown
                    #elif (zygosity == 'complex'):  #what are these?
                    #    genoparts[gene_id].append('complex')
                    geno_hash[genotype_id] = genoparts
                else:
                    #if the gene is not known, still need to add the allele to the genotype hash
                    #these will be added as sequence alterations.
                    genoparts[allele_id] = [allele_id]
                    if zygosity == 'homozygous':
                        genoparts[allele_id].append(allele_id)
                    elif zygosity == 'unknown':
                        genoparts[allele_id].append('?')
                    #elif zygosity == 'complex':  #not sure what to do with these?
                    #    genoparts[allele_id].append('complex')

                    geno_hash[allele_id] = genoparts

                if (limit is not None and line_counter > limit):
                    break

                #end loop through file
            #now loop through the geno_hash, and build the vslcs

            for gt in geno_hash:
                for gene_id in geno_hash.get(gt):
                    variant_locus_parts = geno_hash.get(gt).get(gene_id)
                    if gene_id in variant_locus_parts:
                        #reset the gene_id to none
                        gene_id = None

                    allele1_id = variant_locus_parts[0]
                    allele2_id = None
                    zygosity_id = None
                    #making the assumption that there are not more than 2 variant_locus_parts
                    if len(variant_locus_parts) > 1:
                        allele2_id = variant_locus_parts[1]
                    if allele2_id is not None:
                        if allele2_id == '?':
                            zygosity_id = geno.zygosity['indeterminate']
                            allele2_id = None
                        elif allele2_id == 'complex':
                            pass #not sure what to assign here
                        elif allele1_id != allele2_id:
                            zygosity_id = geno.zygosity['heterozygous']
                        elif allele1_id == allele2_id:
                            zygosity_id = geno.zygosity['homozygous']
                    else:
                        zygosity_id = geno.zygosity['indeterminate']

                    #create the vslc
                    if gene_id is None:
                        g = ''
                    else:
                        g = gene_id
                    if (allele2_id is None):
                        a2 = ''
                    else:
                        a2 = allele2_id

                    vslc_id = self.make_id(('-').join((g,allele1_id,a2)))
                    geno.addPartsToVSLC(vslc_id,allele1_id,allele2_id,zygosity_id)
                    geno.addVSLCtoParent(vslc_id,gt)


            print("INFO: Done with genotypes")
        return

    def _map_allele_type_to_geno(self, allele_type):
        type = 'SO:0001059'  #default: sequence_alteration
        type_map = {
            'complex_substitution': 'SO:1000005',  # complex substitution
            'deficiency': 'SO:1000029',  # incomplete chromosome
            'deletion': 'SO:0000159',  # deletion
            'indel': 'SO:1000032',  #indel
            'insertion': 'SO:0000667',  #insertion
            'point_mutation': 'SO:1000008',  #point_mutation
            'sequence_variant': 'SO:0001060',  #sequence variant
            'transgenic_insertion': 'SO:0001218',  #transgenic insertion
            'transgenic_unspecified': 'SO:0000781',  #transgenic unspecified
            'transloc': 'SO:0000199',  #translocation
            'unspecified' : 'SO:0001059' #sequence alteration
        }
        if (allele_type.strip() in type_map):
            type = type_map.get(allele_type)
        else:
            # TODO add logging
            print("ERROR: Allele Type (", allele_type, ") not mapped")

        return type

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

        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (genotype_id, genotype_name,
                 start_stage_id, start_stage_name,
                 end_stage_id, end_stage_name,
                 subterm1_id, subterm1_name,
                 postcomp1_rel_id, postcomp1_rel_name,
                 superterm1_id, superterm1_name,
                 quality_id, quality_name, modifier,
                 subterm2_id, subterm2_name,
                 postcomp2_rel_id, postcomp2_rel_name,
                 superterm2_id, superterm2_name,
                 pub_id, env_id, empty) = row


                #deal with environments
                #FIXME i am only dealing with 'wild-type' environments for now
                if (not re.match('ZDB-EXP-041102-1', env_id)):
                    print("INFO: Skipping non-wildtype environment", env_id, "for", genotype_id)
                    continue

                genotype_id = 'ZFIN:' + genotype_id.strip()

                geno = Genotype(self.graph)
                geno.addGenotype(genotype_id,genotype_name)
                #because we are using only w.t. environments, the genotype is just intrinsic.

                phenotype_id = self._map_sextuple_to_phenotype(superterm1_id, subterm1_id, quality_id,
                                                               superterm2_id, subterm2_id, modifier)

                if (phenotype_id is None):
                    continue

                #add abnormal phenotypes
                if (not re.match('^normal', modifier)):
                    assoc_id = self.make_id((genotype_id+env_id+phenotype_id+pub_id))
                    pub_id = 'ZFIN:' + pub_id.strip()
                    assoc = G2PAssoc(assoc_id, genotype_id, phenotype_id, pub_id, eco_id)
                    self.graph = assoc.addAssociationNodeToGraph(self.graph)
                else:
                    #add normal phenotypes
                    print("WARN: found normal phenotype; skipping for now")

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_pubinfo(self, raw, out, g, limit=None):
        '''
        This will pull the zfin internal publication information, and map them to their equivalent
        pmid, and make labels.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        '''
        line_counter = 0
        cu = CurieUtil(curie_map.get())
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r', encoding="latin-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pub_id, pubmed_id, authors, title, journal, year, vol, pages, empty) = row

                pub_id = 'ZFIN:'+pub_id.strip()
                pub_label = ('; ').join((authors, title, journal, year, vol, pages))
                gu.addIndividualToGraph(self.graph,pub_id,pub_label)

                pubmed_id = 'PMID:'+pubmed_id.strip()
                if (pubmed_id != '' and pubmed_id is not None):
                    gu.addIndividualToGraph(self.graph,pubmed_id,None)
                    gu.addSameIndividual(self.graph,pub_id,pubmed_id)

                if (limit is not None and line_counter > limit):
                    break

        return


    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status

    def _map_sextuple_to_phenotype(self, superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier):
        '''
        This will take the 6-part EQ-style annotation used by ZFIN and return the ZP id.
        Currently relies on an external mapping file, but the method may be swapped out in the future
        :param superterm1_id:
        :param subterm1_id:
        :param quality_id:
        :param superterm2_id:
        :param subterm2_id:
        :param modifier:
        :return: ZP id
        '''
        zp_id = None
        #FIXME hardcode
        mod_id=modifier
        #zfin uses free-text modifiers, but we need to convert them to proper PATO classes for the mapping
        modifiers = {
            'abnormal' : 'PATO:0000460',
            'normal' : 'PATO:0000461'
        }
        if (modifier in modifiers.keys()):
            mod_id = modifiers.get(modifier)

        key = self._make_zpkey(superterm1_id,subterm1_id,quality_id,superterm2_id,subterm2_id,mod_id)
        mapping = self.zp_map.get(key)

        if (mapping is None):
            print("WARN: Couldn't map ZP id to",("_").join((superterm1_id,subterm1_id,quality_id,superterm2_id,subterm2_id,mod_id)))
        else:
            zp_id = mapping['zp_id']

        return zp_id


    def _load_zp_mappings(self):
        '''
        Given a file that defines the mapping between ZFIN-specific EQ definitions and the automatically
        derived ZP ids, create a mapping here.
        This may be deprecated in the future
        :return:
        '''
        self.zp_map = {}
        print("Loading ZP-to-EQ mappings")
        line_counter = 0
        file=('/').join((self.rawdir,self.files['zpmap']['file']))
        with open(file, 'r', encoding="utf-8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (zp_id, zp_label, superterm1_id, subterm1_id,
                 quality_id, modifier, superterm2_id, subterm2_id) = row
                key = self._make_zpkey(superterm1_id,subterm1_id,quality_id,superterm2_id,subterm2_id,modifier)
                self.zp_map[key] = {
                    'zp_id' : zp_id,
                    'label' : zp_label,
                    'superterm1_id' : superterm1_id,
                    'subterm1_id' : subterm1_id,
                    'quality_id' : quality_id,
                    'modifier' : modifier,
                    'superterm2_id' : superterm2_id,
                    'subterm2_id' : subterm2_id,
                }
        print("Loaded",self.zp_map.__len__(),"zp terms")

        return

    def _make_zpkey(self,superterm1_id,subterm1_id,quality_id,superterm2_id,subterm2_id,modifier):
        key = self.make_id(('_').join((superterm1_id,subterm1_id,quality_id,superterm2_id,subterm2_id,modifier)))
        return key
