import os, csv
from stat import *
import re
from datetime import datetime
import gzip
import os.path
import unicodedata
from dipper.utils import pysed

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Assoc import Assoc
from dipper.models.Genotype import Genotype
from dipper.models.G2PAssoc import G2PAssoc
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper.models.GenomicFeature import Feature,makeChromID


class ClinVar(Source):
    """
    ClinVar is a host of clinically relevant variants, both directly-submitted and curated from the literature.
    We process the variant_summary file here, which is a digested version of their full xml.  We add all
    variants (and coordinates/build) from their system.
    """

    files = {
        'variant_summary' : {
            'file' : 'variant_summary.txt.gz',
            'url' : 'http://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz'
        },
        'variant_citations' : {
            'file' : 'variant_citations.txt',
            'url' : 'http://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt'
        }
        #TODO work through xml

    }

    testmode = False

    def __init__(self, tax_ids=None, gene_ids=None):
        Source.__init__(self, 'clinvar')

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids
        self.filter = 'taxids'
        self.load_bindings()

        self.dataset = Dataset('ClinVar', 'National Center for Biotechnology Information', 'http://www.ncbi.nlm.nih.gov/clinvar/')
        # data-source specific warnings (will be removed when issues are cleared)


        if self.testmode:
            self.gene_ids = [17151, 100008564, 17005, 11834, 14169]
            self.filter = 'geneids'

        self.properties = Feature.properties

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
        """
        The var_citations file has a bad row in it with > 6 cols.  I will comment these out.

        :return:
        """
        #awk  -F"\t" '{if (NF <= 6) print $1, $2, $3, $4, $5, $6 ; OFS = "\t"}' variant_citations.txt
        f = ('/').join((self.rawdir,self.files['variant_citations']['file']))
        print('INFO: removing the line that has too many cols (^15091)')
        pysed.replace("^15091", '#15091', f)



        return

    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows")

        self.scrub()

        print("Parsing files...")

        self._get_variants(limit)
        self._get_var_citations(limit)

        self.load_core_bindings()
        self.load_bindings()

        print("Done parsing files.")

        return

    def _get_variants(self,limit):
        '''
        Currently loops through the variant_summary file.


        :param limit:
        :return:
        '''
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)

        #not unzipping the file
        print("INFO: Processing Variant records")
        line_counter=0
        myfile=('/').join((self.rawdir,self.files['variant_summary']['file']))
        print("FILE:",myfile)
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                if (re.match('^#',line)):
                    continue
                # AlleleID               integer value as stored in the AlleleID field in ClinVar  (//Measure/@ID in the XML)
                # Type                   character, the type of variation
                # Name                   character, the preferred name for the variation
                # GeneID                 integer, GeneID in NCBI's Gene database
                # GeneSymbol             character, comma-separated list of GeneIDs overlapping the variation
                # ClinicalSignificance   character, comma-separated list of values of clinical significance reported for this variation
                #                          for the mapping between the terms listed here and the integers in the .VCF files, see
                #                          http://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
                # RS# (dbSNP)            integer, rs# in dbSNP
                # nsv (dbVar)            character, the NSV identifier for the region in dbVar
                # RCVaccession           character, list of RCV accessions that report this variant
                # TestedInGTR            character, Y/N for Yes/No if there is a test registered as specific to this variation in the NIH Genetic Testing Registry (GTR)
                # PhenotypeIDs           character, list of db names and identifiers for phenotype(s) reported for this variant
                # Origin                 character, list of all allelic origins for this variation
                # Assembly               character, name of the assembly on which locations are based
                # Chromosome             character, chromosomal location
                # Start                  integer, starting location, in pter->qter orientation
                # Stop                   integer, end location, in pter->qter orientation
                # Cytogenetic            character, ISCN band
                # ReviewStatus           character, highest review status for reporting this measure. For the key to the terms,
                #                            and their relationship to the star graphics ClinVar displays on its web pages,
                #                            see http://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/#interpretation
                # HGVS(c.)               character, RefSeq cDNA-based HGVS expression
                # HGVS(p.)               character, RefSeq protein-based HGVS expression
                # NumberSubmitters       integer, number of submissions with this variant
                # LastEvaluated          datetime, the latest time any submitter reported clinical significance
                # Guidelines             character, ACMG only right now, for the reporting of incidental variation in a Gene
                #                                (NOTE: if ACMG, not a specific to the allele but to the Gene)
                # OtherIDs               character, list of other identifiers or sources of information about this variant
                # VariantID              integer, the value used to build the URL for the current default report,
                #                            e.g. http://www.ncbi.nlm.nih.gov/clinvar/variation/1756/
                #


                (allele_num,allele_type,allele_name,gene_num,gene_symbol,clinical_significance,
                 dbsnp_num,dbvar_num,rcv_num,tested_in_gtr,phenotype_ids,origin,
                 assembly,chr,start,stop,cytogenetic_loc,
                 review_status,hgvs_c,hgvs_p,number_of_submitters,last_eval,guidelines,other_ids,variant_num) = line.split('\t')

                ##### set filter=None in init if you don't want to have a filter
                #if self.filter is not None:
                #    if ((self.filter == 'taxids' and (int(tax_num) not in self.tax_ids))
                #            or (self.filter == 'geneids' and (int(gene_num) not in self.gene_ids))):
                #        continue
                ##### end filter

                #print(line)

                line_counter += 1

                seqalt_id = (':').join(('ClinVarVariant',variant_num))
                gene_id = (':').join(('NCBIGene',gene_num))
                tax_id = 'NCBITaxon:9606'   #all are human, so hardcode
                build_id = (':').join(('NCBIGenome',assembly))

                #make the build
                geno.addReferenceGenome(build_id,assembly,tax_id)

                allele_type_id = self._map_type_of_allele(allele_type)

                if (str(chr) == ''):
                    pass
                    #print(line)
                else:
                    geno.addChromosome(str(chr),tax_id,build_id,assembly)
                    chrinbuild_id = makeChromID(str(chr),build_id)


                #todo clinical significance needs to be mapped to a list of terms
                #first, make the variant:
                f = Feature(seqalt_id,allele_name,allele_type_id)

                if start != '-' and start.strip() != '':
                    f.addFeatureStartLocation(start,chrinbuild_id)
                if stop != '-' and stop.strip() != '':
                    f.addFeatureEndLocation(stop,chrinbuild_id)

                f.addFeatureToGraph(self.graph)

                #add the hgvs as synonyms
                if hgvs_c != '-' and hgvs_c.strip() != '':
                    gu.addSynonym(self.graph,seqalt_id,hgvs_c)
                if hgvs_p != '-' and hgvs_p.strip() != '':
                    gu.addSynonym(self.graph,seqalt_id,hgvs_p)

                #add the dbsnp and dbvar ids as equivalent
                if dbsnp_num != '-':
                    dbsnp_id = 'dbSNP:rs'+dbsnp_num
                    gu.addIndividualToGraph(self.graph,dbsnp_id,None)
                    gu.addEquivalentClass(self.graph,seqalt_id,dbsnp_id)
                if dbvar_num != '-':
                    dbvar_id = 'dbVAR:'+dbvar_num
                    gu.addIndividualToGraph(self.graph,dbvar_id,None)
                    gu.addEquivalentClass(self.graph,seqalt_id,dbvar_id)
                #TODO - not sure if this is right
                #the rcv is like the combo of the phenotype with the variant
                #if rcv_num != '-':
                #    rcv_id = 'ClinVar:'+rcv_num
                #    gu.addIndividualToGraph(self.graph,rcv_id,None)
                #    gu.addEquivalentClass(self.graph,seqalt_id,rcv_id)

                #add the gene
                gu.addClassToGraph(self.graph,gene_id,gene_symbol)

                gu.addTriple(self.graph,seqalt_id,geno.object_properties['is_sequence_variant_instance_of'],gene_id)
                #make the variant locus
                #vl_id = ':'+gene_id+'-'+variant_num
                #geno.addSequenceAlterationToVariantLocus(seqalt_id,vl_id)
                #geno.addAlleleOfGene(seqalt_id,gene_id,geno.properties['has_alternate_part'])

                f.loadAllProperties(self.graph)   #FIXME inefficient

                #parse the list of "phenotypes" which are diseases.  add them as an association
                #;GeneReviews:NBK1440,MedGen:C0392514,OMIM:235200,SNOMED CT:35400008;MedGen:C3280096,OMIM:614193;MedGen:CN034317,OMIM:612635;MedGen:CN169374
                #the list is both semicolon delimited and comma delimited, but i don't know why!
                if phenotype_ids != '-':
                    #trim any leading/trailing semicolons/commas
                    phenotype_ids = re.sub('^[;,]','',phenotype_ids)
                    phenotype_ids = re.sub('[;,]$','',phenotype_ids)
                    pheno_list =  re.split('[,;]',phenotype_ids)
                    #print('list:',pheno_list)
                    for p in pheno_list:
                        if (re.match('Orphanet:ORPHA',p)):
                            p = re.sub('Orphanet:ORPHA','ORPHANET:',p)
                        elif (re.match('SNOMED CT',p)):
                            p = re.sub('SNOMED CT','SNOMED',p)
                        assoc_id = self.make_id(seqalt_id+p.strip())
                        #print('assoc:',(',').join((assoc_id,seqalt_id,p.strip())))
                        assoc = G2PAssoc(assoc_id,seqalt_id,p.strip(),None,None)
                        assoc.addAssociationToGraph(self.graph)

                if other_ids != '-':
                    id_list = other_ids.split(',')
                    #TODO make xrefs

                if (limit is not None and line_counter > limit):
                    break

        return

    def _get_var_citations(self,limit):

        # Generated weekly, the first of the week
        # A tab-delimited report of citations associated with data in ClinVar, connected to the AlleleID, the VariationID, and either rs# from dbSNP or nsv in dbVar.
        #
        # AlleleID          integer value as stored in the AlleleID field in ClinVar  (//Measure/@ID in the XML)
        # VariationID       The identifier ClinVar uses to anchor its default display. (in the XML,  //MeasureSet/@ID)
        # rs			    rs identifier from dbSNP
        # nsv				nsv identifier from dbVar
        # citation_source	The source of the citation, either PubMed, PubMedCentral, or the NCBI Bookshelf
        # citation_id		The identifier used by that source


        gu = GraphUtils(curie_map.get())
        print("INFO: Processing Citations for variants")
        line_counter=0
        myfile=('/').join((self.rawdir,self.files['variant_citations']['file']))
        print("FILE:",myfile)
        with open(myfile, 'r',encoding="utf8") as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')

            for line in filereader:
                #skip comments
                line=line
                if (re.match('^#',line[0])):
                    continue
                (allele_num,variant_num,rs_num,nsv_num,citation_source,citation_id) = line

                line_counter += 1

                #the citation for a variant is made to the allele+variant+rs/nsv id
                #we don't know what the citation is for exactly, other than the variant.  so use mentions

                var_id = 'ClinVarVariant:'+variant_num

                #citation source: PubMed | PubMedCentral | citation_source
                #citation id:
                #format the citation id:
                ref_id = None
                if citation_source == 'PubMed':
                    ref_id = 'PMID:'+str(citation_id)
                elif citation_source == 'PubMedCentral':
                    ref_id = 'PMCID:'+str(citation_id)
                if ref_id is not None:
                    gu.addTriple(self.graph,ref_id,self.properties['is_about'],var_id)

                if (limit is not None and line_counter > limit):
                    break



        return

    def _map_type_of_allele(self,type):
        so_id = 'SO:0001060'
        type_to_so_map = {
            'NT expansion' : 'SO:1000039',  #direct tandem duplication
            'copy number gain' : 'SO:0001742',
            'copy number loss' : 'SO:0001743',
            'deletion' : 'SO:0000159',
            'duplication' : 'SO:1000035',
            'fusion' : 'SO:0000806',
            'indel' : 'SO:1000032',
            'insertion' : 'SO:0000667',
            'inversion' : 'SO:1000036',
            'protein only' : 'SO:0001580',   #coding sequence variant.  check me
            'short repeat' : 'SO:0000657',   #repeat region - not sure if this is what's intended.
            'single nucleotide variant' : 'SO:0001483',
            'structural variant' : 'SO:0001537',
            'undetermined variant' : 'SO:0001060'    #sequence variant

        }

        if (type in type_to_so_map):
            so_id = type_to_so_map.get(type)
        else:
            print("WARN: unmapped code",type,". Defaulting to 'SO:sequence_variant'.")

        return so_id


        return so_id

    def _cleanup_id(self,i):
        cleanid = i
        #MIM:123456 --> #OMIM:123456
        cleanid = re.sub('^MIM','OMIM',cleanid)

        #HGNC:HGNC --> HGNC
        cleanid = re.sub('^HGNC:HGNC','HGNC',cleanid)

        #Ensembl --> ENSEMBL
        cleanid = re.sub('^Ensembl','ENSEMBL',cleanid)

        #MGI:MGI --> MGI
        cleanid = re.sub('^MGI:MGI','MGI',cleanid)

        return cleanid


    def remove_control_characters(self,s):
        return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")