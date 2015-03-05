import os
from stat import *
import re
from datetime import datetime
import gzip
import os.path
import logging

from dipper.sources.Source import Source
from dipper.models.GenomicFeature import Feature,makeChromID
from dipper.models.Dataset import Dataset
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


class UCSCBands(Source):
    '''
    This will take the UCSC defintions of cytogenic bands and create the nested structures to enable
    overlap and containment queries.
    For example,
    13q21.31 ==>  13q21.31,  13q21.3,  13q21,  13q2,  13q, 13
    We leverage the Faldo model here for region definitions, and map each of the chromosomal parts to SO.
    At the moment, this only computes the bands for Human.
    We differentiate the species by adding the species id to the identifier prior to the chromosome number.
    The identifers therefore are generically created like:
    <species number>chr<num><band>
    We will create triples for a given band like:
    :9606chr1p36.33 rdf[type] SO:chromosome_band, faldo:Region
    :9606chr1p36 subsequence_of :9606chr1p36.3
    :9606chr1p36 faldo:location [ faldo:begin 0 faldo:end 2300000]

    where any band in the file is an instance of a chr_band (or a more specific type), is a subsequence
    of it's containing region, and is located in the specified coordinates.

    we determine the containing regions of the band by parsing the band-string; since each alphanumeric
    is a significant "place", we can split it with the shorter strings being parents of the longer string

    TODO: will can then sort the locations of the annotated bands, and propagate them to the
    intermediate/parent regions

    TODO: any species by commandline argument

    TODO: abstract this out into a model
    '''

    files = {
        '9606' : {
            'file' : '9606cytoBand.txt.gz',
            'url' : 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz',
        },
    }


    # relationships = {
    #     'gene_product_of' : 'RO:0002204',
    #     'has_gene_product' : 'RO:0002205',
    #     'is_about' : 'IAO:00000136',
    #     'has_subsequence' : 'RO:0002524',
    #     'is_subsequence_of' : 'RO:0002525',
    # }


    def __init__(self, tax_ids=None):
        super().__init__('ucscbands')

        self.tax_ids = tax_ids
        self.load_bindings()
        self.gu = GraphUtils(curie_map.get())
        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606]

        self._check_tax_ids()

        self.dataset = Dataset('ucscbands', 'UCSC Cytogenic Bands', 'http://hgdownload.cse.ucsc.edu')


        #data-source specific warnings (will be removed when issues are cleared)

        return

    def fetch(self, is_dl_forced):

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

        for taxon in self.tax_ids:
            self._get_chrbands(limit, str(taxon))
        #TODO as a post-processing step, we must propagate the coordinates to the upper-level features

        self.load_core_bindings()
        self.load_bindings()


        logger.info("Done parsing files.")

        return

    def _get_chrbands(self,limit,taxon):
        '''
        :param limit:
        :return:
        '''
        line_counter=0
        myfile=('/').join((self.rawdir,self.files[taxon]['file']))
        logger.info("Processing Chr bands from FILE: %s", myfile)

        mybands = {}


        #add the build - ideally it would be taken from the file itself, but it's not in it.
        build_id = 'UCSC:hg19'
        f = Feature(build_id,build_id,Feature.types['reference_genome'])
        f.addTaxonToFeature(self.graph,'NCBITaxon:9606')
        f.addFeatureToGraph(self.graph)
        #f.loadAllProperties(self.graph)

        with gzip.open(myfile, 'rb') as f:
            for line in f:
                #skip comments
                line=line.decode().strip()
                if (re.match('^#',line)):
                    continue

                #chr13	4500000	10000000	p12	stalk
                (chrom,start,stop,band,rtype) = line.split('\t')
                line_counter += 1

                cid = makeChromID(chrom,taxon)
                tax_id = (':').join(('NCBITaxon',taxon))
                region_type_id = self._map_type_of_region(rtype)
                self.gu.addMember(self.graph,build_id,cid)

                #add the chr
                cfeature = Feature(cid,chrom,self._map_type_of_region('chromosome'))
                cfeature.addFeatureStartLocation(0,build_id)
                cfeature.addFeatureToGraph(self.graph) #fixme - should probably only add this at the end?
                cfeature.addTaxonToFeature(self.graph,tax_id)
                #add the chr to the hashmap of coordinates
                if chrom not in mybands.keys():
                    mybands[chrom] = {'min' : 0, 'max' : 0, 'chr' : build_id }

                #add the region and it's location
                maploc_id = cid+band
                bfeature = Feature(maploc_id,chrom+band,region_type_id,rtype)
                bfeature.addFeatureStartLocation(start,cid)
                bfeature.addFeatureEndLocation(stop,cid)
                bfeature.addFeatureToGraph(self.graph)

                #add the staining intensity of the band
                if re.match('g(neg|pos|var)',rtype):
                    self.gu.addTriple(self.graph,maploc_id,Feature.properties['has_staining_intensity'],Feature.types.get(rtype))


                #get the parent bands, and make them unique
                parents = list(self._make_parent_bands(band,set()))
                #alphabetical sort will put them in smallest to biggest
                parents.sort(reverse=True)
                #print('parents of',chrom,band,':',parents)

                #add the parents to the graph, in hierarchical order
                #TODO this is somewhat inefficient due to re-adding upper-level nodes when iterating over the file
                for i in range(len(parents)):
                    pid = cid+parents[i]
                    if (re.match('[pq]$',parents[i])):
                        rti = self._map_type_of_region('chromosome_arm')
                    if (re.match('p$',parents[i])):
                        rti = Feature.types['short_chromosome_arm']
                    elif (re.match('q$',parents[i])):
                        rti = Feature.types['long_chromosome_arm']
                    elif (re.match('[pq]\d$',parents[i])):
                        rti = Feature.types['chromosome_region']
                    elif (re.match('[pq]\d\d',parents[i])):
                        rti = Feature.types['chromosome_band']
                    elif (re.match('[pq]\d\d\.\d+',parents[i])):
                        rti = Feature.types['chromosome_subband']
                    else:
                        rti = self._map_type_of_region('chromosome_part')

                    pfeature = Feature(pid,chrom+parents[i],rti)
                    pfeature.addFeatureToGraph(self.graph)
                    #add the relationships to the parent
                    if (i < len(parents)-1):
                        ppid = cid+parents[i+1]
                        pfeature.addSubsequenceOfFeature(self.graph,ppid)
                    else:
                        #add the last one (p or q usually) as attached to the chromosome
                        pfeature.addSubsequenceOfFeature(self.graph,cid)

                #connect the band here to the first one in the parent list
                bfeature.addSubsequenceOfFeature(self.graph,cid+parents[0])

                #Here, we add the parents to a hashmap of chr bands to propagate the chromosomal coords
                for p in parents:
                    k = chrom+p
                    sta=int(start)
                    sto=int(stop)
                    if k not in mybands.keys():
                        b = {'min' : min(sta,sto), 'max' : max(sta,sto), 'chr' : chrom}
                        mybands[k] = b
                    else:
                        b = mybands.get(k)
                        b['min'] = min(sta,sto,b['min'])
                        b['max'] = max(sta,sto,b['max'])
                        mybands[k] = b
                        #also, set the max for the chrom
                        c = mybands.get(chrom)
                        c['max'] = max(sta,sto,c['max'])
                        mybands[chrom] = c

                if (limit is not None and line_counter > limit):
                    break

        #add the band locations to the graph.
        for b in mybands.keys():
            bid = makeChromID(b,taxon)
            chr = makeChromID(mybands.get(b)['chr'],taxon)
            bfeature = Feature(bid,None,None)  #for now, hopefully it won't overwrite
            bfeature.addFeatureStartLocation(mybands.get(b)['min'],chr)
            bfeature.addFeatureEndLocation(mybands.get(b)['max'],chr)
            bfeature.addFeatureToGraph(self.graph)

        #TODO figure out the staining intensities for the encompassing bands

        return

    def _make_parent_bands(self,band,child_bands):
        '''
        #this will determine the grouping bands that it belongs to, recursively
        #13q21.31 ==>  13, 13q, 13q2, 13q21, 13q21.3, 13q21.31

        :param chrom:
        :param band:
        :param child_bands:
        :return:
        '''
        m=re.match('([pq]\d+(?:\.\d+)?)',band)
        if (len(band) > 0):
            if (m):
                p=str(band[0:len(band)-1])
                p = re.sub('\.$','',p)
                if p is not None:
                    child_bands.add(p)
                    self._make_parent_bands(p,child_bands)
        else:
            child_bands = set()
        return child_bands


    def _map_type_of_region(self,type):
        '''
        Note that "stalk" refers to the short arm of acrocentric chromosomes chr13,14,15,21,22 for human.
        :param type:
        :return:
        '''
        so_id = 'SO:0000830'
        types = Feature.types
        type_to_so_map = {
            'acen' : types['centromere'],
            'gvar' : types['chromosome_part']  , #chromosomal structural element
            'stalk' : types['short_chromosome_arm'], #FIXME using chromosome part for now
            'gneg' : types['chromosome_band'],
            'gpos100' : types['chromosome_band'],
            'gpos25' : types['chromosome_band'],
            'gpos50' : types['chromosome_band'],
            'gpos75' : types['chromosome_band'],
            'chromosome' : types['chromosome'],
            'chromosome_arm' : types['chromosome_arm'],
            'chromosome_band' : types['chromosome_band'],
            'chromosome_part' : types['chromosome_part']
        }

        if (type in type_to_so_map):
            so_id = type_to_so_map.get(type)
        else:
            logger.warn("Unmapped code %s. "
                        "Defaulting to chr_part 'SO:0000830'.", type)

        return so_id

    def _check_tax_ids(self):
        for taxon in self.tax_ids:
            if str(taxon) not in self.files:
                raise Exception("Taxon " + str(taxon) + " not supported"
                                " by source UCSCBands")


