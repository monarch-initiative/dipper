import os
from stat import *
import re
from datetime import datetime
import gzip
import os.path
import unicodedata
import logging

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Assoc import Assoc
from dipper.models.Genotype import Genotype

from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper import config
from dipper.models.GenomicFeature import Feature,makeChromID

logger = logging.getLogger(__name__)


class NCBIGene(Source):
    """
    This is the processing module for the National Center for Biotechnology Information.  It includes parsers
    for the gene_info (gene names, symbols, ids, equivalent ids), gene history (alt ids), and
    gene2pubmed publication references about a gene.

    This creates Genes as classes, when they are properly typed as such.  For those entries where it is an
    'unknown significance', it is added simply as an instance of a sequence feature.  It will add equivalentClasses
    for a subset of external identifiers, including:  ENSEMBL, HGMD, MGI, ZFIN, and gene product links for HPRD.
    They are additionally located to their Chromosomal band (until we process actual genomic coords in
    a separate file).

    We process the genes from the filtered taxa, starting with those configured by default (human, mouse, fish).
    This can be overridden in the calling script to include additional taxa, if desired.
    The gene ids in the conf.json will be used to subset the data when testing.

    All entries in the gene_history file are added as deprecated classes, and linked to the current gene id, with
    "replaced_by" relationships.

    Since we do not know much about the specific link in the gene2pubmed; we simply create a "mentions" relationship.

    """

    files = {
        'gene_info' : {
            'file' : 'gene_info.gz',
            'url' : 'http://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'
        },
        'gene_history' : {
            'file': 'gene_history.gz',
            'url' : 'http://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz'
        },
        'gene2pubmed' : {
            'file': 'gene2pubmed.gz',
            'url' : 'http://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz'
        },
    }


    # relationships = {
    #     'gene_product_of' : 'RO:0002204',
    #     'has_gene_product' : 'RO:0002205',
    #     'is_about' : 'IAO:0000136'
    # }


    def __init__(self, tax_ids=None, gene_ids=None):
        Source.__init__(self, 'ncbigene')

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids
        self.filter = 'taxids'
        self.load_bindings()

        self.dataset = Dataset('ncbigene', 'National Center for Biotechnology Information', 'http://ncbi.nih.nlm.gov/gene',None,'http://www.ncbi.nlm.nih.gov/About/disclaimer.html','https://creativecommons.org/publicdomain/mark/1.0/')
        # data-source specific warnings (will be removed when issues are cleared)

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]

        #if self.testmode:
        #    self.gene_ids = [17151, 100008564, 17005, 11834, 14169]
        #    self.filter = 'geneids'
        self.gene_ids = []
        if (not (('test_ids' in config.get_config()) and ('gene' in config.get_config()['test_ids']))):
            logger.warn("not configured with gene test ids.")
        else:
            self.gene_ids = config.get_config()['test_ids']['gene']


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



    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._get_gene_info(limit)
        self._get_gene_history(limit)
        self._get_gene2pubmed(limit)

        self.load_core_bindings()
        self.load_bindings()

        logger.info("Done parsing files.")

        return

    def _get_gene_info(self,limit):
        '''
        Currently loops through the gene_info file and creates the genes as classes, typed with SO.  It will add their
        label, any alternate labels as synonyms, alternate ids as equivlaent classes.  HPRDs get added as
        protein products.  The chromosome and chr band get added as blank node regions, and the gene is faldo:located
        on the chr band.
        :param limit:
        :return:
        '''
        gu = GraphUtils(curie_map.get())

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)

        #not unzipping the file
        logger.info("Processing Gene records")
        line_counter=0
        myfile=('/').join((self.rawdir,self.files['gene_info']['file']))
        logger.info("FILE: %s",myfile)
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                if (re.match('^#',line)):
                    continue
                (tax_num,gene_num,symbol,locustag,
                 synonyms,xrefs,chr,map_loc,desc,
                 gtype,authority_symbol,name,
                 nomenclature_status,other_designations,modification_date) = line.split('\t')

                ##### set filter=None in init if you don't want to have a filter
                #if self.filter is not None:
                #    if ((self.filter == 'taxids' and (int(tax_num) not in self.tax_ids))
                #            or (self.filter == 'geneids' and (int(gene_num) not in self.gene_ids))):
                #        continue
                ##### end filter


                if ((self.testMode and (int(gene_num) not in self.gene_ids)) or
                    (not self.testMode and (int(tax_num) not in self.tax_ids))):
                    continue

                line_counter += 1

                gene_id = (':').join(('NCBIGene',gene_num))
                tax_id = (':').join(('NCBITaxon',tax_num))
                gene_type_id = self._map_type_of_gene(gtype)

                if (symbol == 'NEWENTRY'):
                    label = None
                else:
                    label = symbol

                #TODO might have to figure out if things aren't genes, and make them individuals
                gu.addClassToGraph(g,gene_id,label,gene_type_id,desc)

                #we have to do special things here for genes, because they're classes not individuals
                #f = Feature(gene_id,label,gene_type_id,desc)

                if (name != '-'):
                    gu.addSynonym(g,gene_id,name)
                if (synonyms.strip() != '-'):
                    for s in synonyms.split('|'):
                        gu.addSynonym(g,gene_id,s.strip(),Assoc.properties['hasRelatedSynonym'])
                if (other_designations.strip() != '-'):
                    for s in other_designations.split('|'):
                        gu.addSynonym(g,gene_id,s.strip(),Assoc.properties['hasRelatedSynonym'])


                #deal with the xrefs
                #MIM:614444|HGNC:HGNC:16851|Ensembl:ENSG00000136828|HPRD:11479|Vega:OTTHUMG00000020696
                if (xrefs.strip() != '-'):
                    for r in xrefs.strip().split('|'):
                        fixedr = self._cleanup_id(r)
                        if ((fixedr is not None) and (fixedr.strip() != '')):
                            if (re.match('HPRD',fixedr)):
                                #proteins are not == genes.
                                gu.addTriple(g,gene_id,self.properties['has_gene_product'],fixedr)
                            else:
                                if (fixedr.split(':')[0] not in ['Vega','IMGT/GENE-DB']):  #skip these for now
                                    gu.addEquivalentClass(g,gene_id,fixedr)

                if (str(chr) != '-' and str(chr) != ''):
                    if (re.search('\|',str(chr))):
                        #this means that there's uncertainty in the mapping.  skip it
                        #TODO we'll need to figure out how to deal with >1 loc mapping
                        logger.info('%s is non-uniquely mapped to %s.  Skipping for now.', gene_id, str(chr))
                        continue
                    #if (not re.match('(\d+|(MT)|[XY]|(Un)$',str(chr).strip())):
                    #    print('odd chr=',str(chr))
                    geno.addGenome(tax_id,tax_num)   #tax label can get added elsewhere
                    geno.addChromosome(str(chr),tax_id,tax_num)  #temporarily use the taxnum for the disambiguating label
                    mychrom = makeChromID(str(chr),tax_id)
                    if (tax_num == '9606' and map_loc != '-'):
                        #this matches the regular kind of chrs, so make that kind of band
                        #not sure why this matches? chrX|Y or 10090chr12|Un"
                        #TODO we probably need a different regex per organism
                        if re.match('[0-9A-Z]+[pq](\d+)?(\.\d+)?$',map_loc):
                            #the maploc_id already has the numeric chromosome in it, strip it first
                            bid = re.sub('^'+str(chr),'',map_loc)
                            maploc_id = makeChromID(str(chr)+bid,tax_num)  #the generic location (no coordinates)
                            #print(map_loc,'-->',bid,'-->',maploc_id)
                            band = Feature(maploc_id,None,None)  #Assume it's type will be added elsewhere
                            band.addFeatureToGraph(g)
                            #add the band as the containing feature
                            gu.addTriple(g,gene_id,Feature.object_properties['is_subsequence_of'],maploc_id)
                        else:
                            gu.addTriple(g,gene_id,Feature.object_properties['is_subsequence_of'],mychrom)
                            #TODO handle these cases
                            #examples are: 15q11-q22, Xp21.2-p11.23, 15q22-qter, 10q11.1-q24,
                            ## 12p13.3-p13.2|12p13-p12, 1p13.3|1p21.3-p13.1,  12cen-q21, 22q13.3|22q13.3
                            logger.info('not regular band pattern for %s: %s', gene_id, map_loc)
                    else:
                        #add the gene as a subsequence of the chromosome
                        gu.addTriple(g,gene_id,Feature.object_properties['is_subsequence_of'],mychrom)
                else:
                    #since the gene is unlocated, add the taxon as a property of it
                    #TODO should we add the gene to the "genome" instead of letting it dangle?
                    geno.addTaxon(tax_id,gene_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

            gu.loadProperties(g,Feature.object_properties,gu.OBJPROP)
            gu.loadProperties(g,Feature.data_properties,gu.DATAPROP)
            gu.loadProperties(g,Genotype.object_properties,gu.OBJPROP)
            gu.loadAllProperties(g)


        return


    def _get_gene_history(self,limit):
        '''
        Loops through the gene_history file and adds the old gene ids as deprecated classes, where the new
        gene id is the replacement for it.  The old gene symbol is added as a synonym to the gene.
        :param limit:
        :return:
        '''
        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        logger.info("Processing Gene records")
        line_counter=0
        myfile=('/').join((self.rawdir,self.files['gene_history']['file']))
        logger.info("FILE: %s",myfile)
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                if (re.match('^#',line)):
                    continue
                (tax_num,gene_num,discontinued_num,discontinued_symbol,discontinued_date) = line.split('\t')

                ##### set filter=None in init if you don't want to have a filter
                #if self.filter is not None:
                #    if ((self.filter == 'taxids' and (int(tax_num) not in self.tax_ids))
                #            or (self.filter == 'geneids' and (int(gene_num) not in self.gene_ids))):
                #        continue
                ##### end filter


                if (gene_num == '-' or discontinued_num =='-'):
                    continue

                if ((self.testMode and (int(gene_num) not in self.gene_ids)) or
                    (not self.testMode and (int(tax_num) not in self.tax_ids))):
                    continue

                line_counter += 1
                gene_id = (':').join(('NCBIGene',gene_num))
                discontinued_gene_id = (':').join(('NCBIGene',discontinued_num))
                tax_id = (':').join(('NCBITaxon',tax_num))

                #add the two genes
                gu.addClassToGraph(g,gene_id,None)
                gu.addClassToGraph(g,discontinued_gene_id,discontinued_symbol)

                #add the new gene id to replace the old gene id
                gu.addDeprecatedClass(g,discontinued_gene_id,[gene_id])

                #also add the old symbol as a synonym of the new gene
                gu.addSynonym(g,gene_id,discontinued_symbol)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _get_gene2pubmed(self,limit):
        '''
        Loops through the gene2pubmed file and adds a simple triple to say that a given publication
        is_about a gene.  Publications are added as NamedIndividuals.
        :param limit:
        :return:
        '''

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        is_about = gu.getNode(Assoc.properties['is_about'])

        logger.info("Processing Gene records")
        line_counter=0
        myfile=('/').join((self.rawdir,self.files['gene2pubmed']['file']))
        logger.info("FILE: %s",myfile)
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                if (re.match('^#',line)):
                    continue
                (tax_num,gene_num,pubmed_num) = line.split('\t')

                ##### set filter=None in init if you don't want to have a filter
                #if self.filter is not None:
                #    if ((self.filter == 'taxids' and (int(tax_num) not in self.tax_ids))
                #       or (self.filter == 'geneids' and (int(gene_num) not in self.gene_ids))):
                #        continue
                ##### end filter

                if ((self.testMode and (int(gene_num) not in self.gene_ids)) or
                    (not self.testMode and (int(tax_num) not in self.tax_ids))):
                    continue


                if (gene_num == '-' or pubmed_num =='-'):
                    continue
                line_counter += 1
                gene_id = (':').join(('NCBIGene',gene_num))
                pubmed_id = (':').join(('PMID',pubmed_num))

                #add the gene, in case it hasn't before
                gu.addClassToGraph(g,gene_id,None)
                #add the publication as a NamedIndividual
                gu.addIndividualToGraph(g,pubmed_id,None,None)  #add type publication
                self.graph.add((gu.getNode(pubmed_id),is_about,gu.getNode(gene_id)))

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _map_type_of_gene(self,type):
        so_id = 'SO:0000704'
        type_to_so_map = {
            'ncRNA': 'SO:0001263',
            'other': 'SO:0000704',
            'protein-coding': 'SO:0001217',
            'pseudo': 'SO:0000336',
            'rRNA': 'SO:0001637',
            'snRNA': 'SO:0001268',
            'snoRNA': 'SO:0001267',
            'tRNA': 'SO:0001272',
            'unknown': 'SO:0000704',
            'scRNA' : 'SO:0000013',
            'miscRNA' : 'SO:0000233', #mature transcript - there is no good mapping
            'chromosome' : 'SO:0000340',
            'chromosome_arm' : 'SO:0000105',
            'chromosome_band' : 'SO:0000341',
            'chromosome_part' : 'SO:0000830'
        }

        if (type in type_to_so_map):
            so_id = type_to_so_map.get(type)
        else:
            logger.warn("unmapped code %s. Defaulting to 'SO:0000704'.", type)

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


    #TODO move to util
    def remove_control_characters(self,s):
        return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")



    def getTestSuite(self):
        import unittest
        from tests.test_ncbi import NCBITestCase
        #TODO test genes
        #from tests.test_genotypes import GenotypeTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(NCBITestCase)

        return test_suite
