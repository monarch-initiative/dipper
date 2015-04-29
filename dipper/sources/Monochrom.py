import re
import gzip
import logging

from dipper.sources.Source import Source
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.models.Dataset import Dataset
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Genotype import Genotype
from dipper import curie_map

logger = logging.getLogger(__name__)


class Monochrom(Source):
    """
    This class will leverage the GENO ontology and modeling patterns to build an ontology of chromosomes
    for any species. These classes represent major structural pieces of Chromosomes which are often universally
    referenced, using physical properties/observations that remain constant across different genome builds
    (such as banding patterns and arms). The idea is to create a scaffold upon which we can hang
    build-specific chromosomal coordinates, and reason across them.

    In general, this will take the cytogenic bands files from UCSC, and create missing grouping classes, in order
    to build the partonomy from a very specific chromosomal band up through the chromosome itself and enable
    overlap and containment queries.  We use RO:subsequence_of as our relationship between nested chromosomal parts.
    For example,
    13q21.31 ==>  13q21.31,  13q21.3,  13q21,  13q2,  13q, 13

    At the moment, this only computes the bands for Human, Mouse, Zebrafish, and Rat,
    but will be expanding in the future as needed.

    Because this is a universal framework to represent the chromosomal structure of any species, we must
    mint identifiers for each chromosome and part.  We differentiate species by first creating a
    species-specific genome, then for each species-specific chromosome we include the NCBI taxon number together with
    the chromosome number, like:  ```<species number>chr<num><band>```.  For 13q21.31, this would be
    9606chr13q21.31.
    We then create triples for a given band like:
    <pre>
    CHR:9606chr1p36.33 rdf[type] SO:chromosome_band
    CHR:9606chr1p36 subsequence_of :9606chr1p36.3
    </pre>
    where any band in the file is an instance of a chr_band (or a more specific type), is a subsequence
    of it's containing region.

    We determine the containing regions of the band by parsing the band-string; since each alphanumeric
    is a significant "place", we can split it with the shorter strings being parents of the longer string

    Since this is small, and we have not limited other items in our test set to a small region, we
    simply use the whole graph (genome) for testing purposes, and copy the main graph to the test graph.

    Remember, this Dipper class is building an ONTOLOGY, therefore we must also include domain and range
    constraints, and other owl-isms.  Some of these functions can be generalized and

    TODO: any species by commandline argument

    We are currently mapping these to the **CHR idspace**, but this is NOT YET APPROVED and is subject to change.
    """

    files = {
        '9606': {
            'file': '9606cytoBand.txt.gz',
            'url': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz',
            'build_num': 'hg19',
            'genome_label': 'Human'
        },
        # Note that there are no bands, arms or staining components for zfish at the moment
        '7955': {
            'file': '7955cytoBand.txt.gz',
            'url': 'http://hgdownload.cse.ucsc.edu/goldenPath/danRer10/database/cytoBandIdeo.txt.gz',
            'build_num': 'danRer10',
            'genome_label': 'Zebrafish'
        },
        '10090': {
            'file': '10090cytoBand.txt.gz',
            'url': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBandIdeo.txt.gz',
            'build_num': 'mm10',
            'genome_label': 'Mouse'
        },
        # Note that there are no bands, arms or staining components for rat at the moment
        '10116': {
            'file': '10116cytoBand.txt.gz',
            'url': 'http://hgdownload.cse.ucsc.edu/goldenPath/rn6/database/cytoBandIdeo.txt.gz',
            'build_num': 'rn6',
            'genome_label': 'Rat'
        },
    }

    def __init__(self, tax_ids=None):
        super().__init__('monochrom')

        self.tax_ids = tax_ids
        self.load_bindings()
        self.gu = GraphUtils(curie_map.get())

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 7955, 10090, 10116]

        self._check_tax_ids()

        # TODO add license
        self.dataset = Dataset('monochrom', 'Monarch Chromosome Ontology', 'http://monarchinitiative.org',
                               None, 'http://creativecommons.org/licenses/by/4.0/')

        # data-source specific warnings (will be removed when issues are cleared)

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):

        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        for taxon in self.tax_ids:
            self._get_chrbands(limit, str(taxon))

        self.load_core_bindings()
        self.load_bindings()

        # using the full graph as the test here
        self.testgraph = self.graph
        logger.info("Found %d nodes", len(self.graph))
        logger.info("Done parsing files.")

        return

    def _get_chrbands(self, limit, taxon):
        """
        For the given taxon, it will fetch the chr band file.  We will not deal with the coordinate information
        as part of this parser.  Here, we only are concerned with building the partonomy.
        :param limit:
        :return:
        """
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files[taxon]['file']))
        logger.info("Processing Chr bands from FILE: %s", myfile)
        geno = Genotype(self.graph)

        # build the organism's genome from the taxon
        genome_label = self.files[taxon]['genome_label']
        taxon_id = 'NCBITaxon:'+taxon
        genome_id = geno.makeGenomeID(taxon_id)
        geno.addGenome(taxon_id, genome_label)

        with gzip.open(myfile, 'rb') as f:
            for line in f:
                # skip comments
                line = line.decode().strip()
                if re.match('^#', line):
                    continue

                # chr13	4500000	10000000	p12	stalk
                (chrom, start, stop, band, rtype) = line.split('\t')
                line_counter += 1

                cclassid = makeChromID(chrom, taxon)  # the chrom class, taxon as the reference

                # add the generic and build-specific chromosome
                geno.addChromosome(chrom, taxon_id, genome_label)
                self.gu.addOWLPropertyClassRestriction(self.graph, cclassid,
                                                       self.gu.object_properties['member_of'], genome_id)

                # add the band(region) as a class
                maplocclass_id = cclassid+band
                maplocclass_label = makeChromLabel(chrom+band, genome_label)
                region_type_id = self._map_type_of_region(rtype)
                self.gu.addClassToGraph(self.graph, maplocclass_id, maplocclass_label, region_type_id)

                # add the staining intensity of the band
                if re.match('g(neg|pos|var)', rtype):
                    stain_type = Feature.types.get(rtype)
                    if stain_type is not None:
                        self.gu.addTriple(self.graph, maplocclass_id, Feature.properties['has_staining_intensity'],
                                          Feature.types.get(rtype))
                    else:
                        logger.warn('staining type not found: %s', rtype)

                # get the parent bands, and make them unique
                parents = list(self._make_parent_bands(band, set()))
                # alphabetical sort will put them in smallest to biggest
                parents.sort(reverse=True)

                # print("PARENTS of",maplocclass_id,"=",parents)
                # add the parents to the graph, in hierarchical order
                # TODO this is somewhat inefficient due to re-adding upper-level nodes when iterating over the file
                for i in range(len(parents)):
                    pclassid = cclassid+parents[i]  # class chr parts
                    pclass_label = makeChromLabel(chrom+parents[i], genome_label)

                    if re.match('p$', parents[i]):
                        rti = Feature.types['short_chromosome_arm']
                    elif re.match('q$', parents[i]):
                        rti = Feature.types['long_chromosome_arm']
                    elif re.match('[pq][A-H\d]$', parents[i]):
                        rti = Feature.types['chromosome_region']
                    elif re.match('[pq][A-H\d]\d', parents[i]):
                        rti = Feature.types['chromosome_band']
                    elif re.match('[pq][A-H\d]\d\.\d+', parents[i]):
                        rti = Feature.types['chromosome_subband']
                    else:
                        rti = self._map_type_of_region('chromosome_part')

                    self.gu.addClassToGraph(self.graph, pclassid, pclass_label, rti)

                    # for canonical chromosomes, then the subbands are subsequences of the full band
                    # add the subsequence stuff as restrictions
                    if i < len(parents) - 1:
                        pid = cclassid+parents[i+1]   # the instance
                        self.gu.addOWLPropertyClassRestriction(self.graph, pclassid,
                                                               Feature.object_properties['is_subsequence_of'], pid)
                        self.gu.addOWLPropertyClassRestriction(self.graph, pid,
                                                               Feature.object_properties['has_subsequence'], pclassid)

                    else:
                        # add the last one (p or q usually) as attached to the chromosome
                        self.gu.addOWLPropertyClassRestriction(self.graph, pclassid,
                                                               Feature.object_properties['is_subsequence_of'], cclassid)
                        self.gu.addOWLPropertyClassRestriction(self.graph, cclassid,
                                                               Feature.object_properties['has_subsequence'], pclassid)

                # connect the band here to the first one in the parent list
                if len(parents) > 0:
                    self.gu.addOWLPropertyClassRestriction(self.graph, maplocclass_id,
                                                           Feature.object_properties['is_subsequence_of'],
                                                           cclassid+parents[0])
                    self.gu.addOWLPropertyClassRestriction(self.graph, cclassid+parents[0],
                                                           Feature.object_properties['has_subsequence'], maplocclass_id)

                if limit is not None and line_counter > limit:
                    break

        self.gu.loadAllProperties(self.graph)

        # TODO figure out the staining intensities for the encompassing bands

        return

    def _make_parent_bands(self, band, child_bands):
        """
        this will determine the grouping bands that it belongs to, recursively
        13q21.31 ==>  13, 13q, 13q2, 13q21, 13q21.3, 13q21.31

        :param band:
        :param child_bands:
        :return:
        """
        m = re.match('([pq][A-H\d]+(?:\.\d+)?)', band)
        if len(band) > 0:
            if m:
                p = str(band[0:len(band)-1])
                p = re.sub('\.$', '', p)
                if p is not None:
                    child_bands.add(p)
                    self._make_parent_bands(p, child_bands)
        else:
            child_bands = set()
        return child_bands

    def _map_type_of_region(self, regiontype):
        """
        Note that "stalk" refers to the short arm of acrocentric chromosomes chr13,14,15,21,22 for human.
        :param regiontype:
        :return:
        """
        so_id = 'SO:0000830'
        types = Feature.types
        type_to_so_map = {
            'acen': types['centromere'],
            'gvar': types['chromosome_band'],
            'stalk': types['chromosome_band'],  # FIXME using chromosome part for now
            'gneg': types['chromosome_band'],
            'gpos100': types['chromosome_band'],
            'gpos25': types['chromosome_band'],
            'gpos33': types['chromosome_band'],
            'gpos50': types['chromosome_band'],
            'gpos66': types['chromosome_band'],
            'gpos75': types['chromosome_band'],
            'chromosome': types['chromosome'],
            'chromosome_arm': types['chromosome_arm'],
            'chromosome_band': types['chromosome_band'],
            'chromosome_part': types['chromosome_part']
        }

        if regiontype in type_to_so_map:
            so_id = type_to_so_map.get(regiontype)
        else:
            logger.warn("Unmapped code %s. "
                        "Defaulting to chr_part 'SO:0000830'.", regiontype)

        return so_id

    def _check_tax_ids(self):
        for taxon in self.tax_ids:
            if str(taxon) not in self.files:
                raise Exception("Taxon " + str(taxon) + " not supported"
                                " by source Monochrom")

    def getTestSuite(self):
        # import unittest
        # from tests.test_ucscbands import UCSCBandsTestCase
        test_suite = None
        # test_suite = unittest.TestLoader().loadTestsFromTestCase(UCSCBandsTestCase)

        return test_suite
