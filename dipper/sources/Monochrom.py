import re
import gzip
import logging

from dipper.sources.Source import Source
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.models.Dataset import Dataset
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model
from dipper.models.Family import Family


logger = logging.getLogger(__name__)
MCDL = 'http://hgdownload.cse.ucsc.edu/goldenPath'


class Monochrom(Source):
    """
    This class will leverage the GENO ontology and modeling patterns to build
    an ontology of chromosomes for any species. These classes represent major
    structural pieces of Chromosomes which are often universally referenced,
    using physical properties/observations that remain constant over different
    genome builds (such as banding patterns and arms). The idea is to create a
    scaffold upon which we can hang build-specific chromosomal coordinates,
    and reason across them.

    In general, this will take the cytogenic bands files from UCSC, and create
    missing grouping classes, in order to build the partonomy from a very
    specific chromosomal band up through the chromosome itself and enable
    overlap and containment queries.  We use RO:subsequence_of as our
    relationship between nested chromosomal parts. For example,
    13q21.31 ==>  13q21.31,  13q21.3,  13q21,  13q2,  13q, 13

    At the moment, this only computes the bands for
    Human, Mouse, Zebrafish, and Rat
    but will be expanding in the future as needed.

    Because this is a universal framework to represent the chromosomal
    structure of any species, we must mint identifiers for each chromosome
    and part. We differentiate species by first creating a species-specific
    genome, then for each species-specific chromosome we include the NCBI taxon
    number together with the chromosome number, like:
    ```<species number>chr<num><band>```.  For 13q21.31, this would be
    9606chr13q21.31.
    We then create triples for a given band like:
    <pre>
    CHR:9606chr1p36.33 rdf[type] SO:chromosome_band
    CHR:9606chr1p36 subsequence_of :9606chr1p36.3
    </pre>
    where any band in the file is an instance of a chr_band
    (or a more specific type), is a subsequence of it's containing region.

    We determine the containing regions of the band by parsing the band-string;
    since each alphanumeric is a significant "place", we can split it with the
    shorter strings being parents of the longer string

    Since this is small, and we have not limited other items in our test set to
    a small region, we simply use the whole graph (genome)
    for testing purposes, and copy the main graph to the test graph.

    Since this Dipper class is building an ONTOLOGY,
    rather than instance-level data, we must also include domain and range
    constraints, and other owl-isms.

    TODO: any species by commandline argument

    We are currently mapping these to the **CHR idspace**,
    but this is NOT YET APPROVED and is subject to change.
    """

    files = {
        '9606': {
            'file': '9606cytoBand.txt.gz',
            'url': MCDL + '/hg19/database/cytoBand.txt.gz',
            'build_num': 'hg19',
            'genome_label': 'Human'
        },
        '10090': {
            'file': '10090cytoBand.txt.gz',
            'url': MCDL + '/mm10/database/cytoBandIdeo.txt.gz',
            'build_num': 'mm10',
            'genome_label': 'Mouse'
        },
        # Note that there are no bands, arms or staining components
        # for the following genomes at the moment
        '7955': {
            'file': '7955cytoBand.txt.gz',
            'url': MCDL + '/danRer10/database/cytoBandIdeo.txt.gz',
            'build_num': 'danRer10',
            'genome_label': 'Zebrafish'
        },
        '10116': {
            'file': '10116cytoBand.txt.gz',
            'url': MCDL + '/rn6/database/cytoBandIdeo.txt.gz',
            'build_num': 'rn6',
            'genome_label': 'Rat'
        },
        '9913': {
            'file': 'bosTau7cytoBand.txt.gz',
            'url': MCDL + '/bosTau7/database/cytoBandIdeo.txt.gz',
            'build_num': 'bosTau7',
            'genome_label': 'cow'
        },
        '9031': {
            'file': 'galGal4cytoBand.txt.gz',
            'url': MCDL + '/galGal4/database/cytoBandIdeo.txt.gz',
            'build_num': 'galGal4',
            'genome_label': 'chicken'
        },
        '9823': {
            'file': 'susScr3cytoBand.txt.gz',
            'url': MCDL + '/susScr3/database/cytoBandIdeo.txt.gz',
            'build_num': 'susScr3',
            'genome_label': 'pig'
        },
        '9940': {
            'file': 'oviAri3cytoBand.txt.gz',
            'url': MCDL + '/oviAri3/database/cytoBandIdeo.txt.gz',
            'build_num': 'oviAri3',
            'genome_label': 'sheep'
        },
        '9796': {
            'file': 'equCab2cytoBand.txt.gz',
            'url': MCDL + '/equCab2/database/cytoBandIdeo.txt.gz',
            'build_num': 'equCab2',
            'genome_label': 'horse'
        },
    }

    region_type_map = {
        'acen': Feature.types['centromere'],
        'gvar': Feature.types['chromosome_band'],
        'stalk': Feature.types['chromosome_band'],
        'gneg': Feature.types['chromosome_band'],
        'gpos100': Feature.types['chromosome_band'],
        'gpos25': Feature.types['chromosome_band'],
        'gpos33': Feature.types['chromosome_band'],
        'gpos50': Feature.types['chromosome_band'],
        'gpos66': Feature.types['chromosome_band'],
        'gpos75': Feature.types['chromosome_band'],
        'chromosome': Feature.types['chromosome'],
        'chromosome_arm': Feature.types['chromosome_arm'],
        'chromosome_band': Feature.types['chromosome_band'],
        'chromosome_part': Feature.types['chromosome_part']
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'monochrom')

        self.tax_ids = tax_ids
        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [
                9606, 10090, 7955, 10116, 9913, 9031, 9823, 9940, 9796]

        self._check_tax_ids()

        # TODO add license
        self.dataset = Dataset(
            'monochrom', 'Monarch Chromosome Ontology',
            'http://monarchinitiative.org', None,
            'http://creativecommons.org/licenses/by/4.0/')

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

        # using the full graph as the test here
        self.testgraph = self.graph
        logger.info("Done parsing files.")

        return

    def _get_chrbands(self, limit, taxon):
        """
        For the given taxon, it will fetch the chr band file.
        We will not deal with the coordinate information with this parser.
        Here, we only are concerned with building the partonomy.
        :param limit:
        :return:

        """
        model = Model(self.graph)
        family = Family(self.graph)
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files[taxon]['file']))
        logger.info("Processing Chr bands from FILE: %s", myfile)
        geno = Genotype(self.graph)

        # build the organism's genome from the taxon
        genome_label = self.files[taxon]['genome_label']
        taxon_id = 'NCBITaxon:'+taxon

        # add the taxon as a class.  adding the class label elsewhere
        model.addClassToGraph(taxon_id, None)
        model.addSynonym(taxon_id, genome_label)

        genome_id = geno.makeGenomeID(taxon_id)
        geno.addGenome(taxon_id, genome_label)
        model.addOWLPropertyClassRestriction(
            genome_id, Genotype.object_properties['in_taxon'],
            taxon_id)

        with gzip.open(myfile, 'rb') as f:
            for line in f:
                # skip comments
                line = line.decode().strip()
                if re.match(r'^#', line):
                    continue

                # chr13	4500000	10000000	p12	stalk
                (chrom, start, stop, band, rtype) = line.split('\t')
                line_counter += 1

                # NOTE
                # some less-finished genomes have placed and unplaced scaffolds
                # * Placed scaffolds:
                #    Scaffold has an oriented location within a chromosome.
                # * Unlocalized scaffolds:
                #     scaffold 's chromosome  is known,
                #     scaffold's position, orientation or both is not known.
                # *Unplaced scaffolds:
                #   it is not known which chromosome the scaffold belongs to.

                # find out if the thing is a full on chromosome, or a scaffold:
                # ex: unlocalized scaffold: chr10_KL568008v1_random
                # ex: unplaced scaffold: chrUn_AABR07022428v1
                placed_scaffold_pattern = r'chr(\d+|X|Y|Z|W|MT|M)'

                # TODO unused
                # unlocalized_scaffold_pattern = \
                #    placed_scaffold_pattern + r'_(\w+)_random'
                # unplaced_scaffold_pattern = r'chrUn_(\w+)'

                m = re.match(placed_scaffold_pattern+r'$', chrom)
                if m is not None and len(m.groups()) == 1:
                    # the chromosome is the first match of the pattern
                    # ch = m.group(1)  # TODO unused
                    pass
                else:
                    # let's skip over anything that isn't a placed_scaffold
                    # at the class level
                    logger.info("Skipping non-placed chromosome %s", chrom)
                    continue
                # the chrom class, taxon as the reference
                cclassid = makeChromID(chrom, taxon, 'CHR')

                # add the chromosome as a class
                geno.addChromosomeClass(chrom, taxon_id, genome_label)
                model.addOWLPropertyClassRestriction(
                    cclassid, family.object_properties['member_of'], genome_id)

                # add the band(region) as a class
                maplocclass_id = cclassid+band
                maplocclass_label = makeChromLabel(chrom+band, genome_label)
                if band is not None and band.strip() != '':
                    region_type_id = self.map_type_of_region(rtype)
                    model.addClassToGraph(
                        maplocclass_id, maplocclass_label,
                        region_type_id)
                else:
                    region_type_id = Feature.types['chromosome']
                # add the staining intensity of the band
                if re.match(r'g(neg|pos|var)', rtype):
                    if region_type_id in [
                            Feature.types['chromosome_band'],
                            Feature.types['chromosome_subband']]:
                        stain_type = Feature.types.get(rtype)
                        if stain_type is not None:
                            model.addOWLPropertyClassRestriction(
                                maplocclass_id,
                                Feature.properties['has_staining_intensity'],
                                Feature.types.get(rtype))
                    else:
                        # usually happens if it's a chromosome because
                        # they don't actually have banding info
                        logger.info("feature type %s != chr band",
                                    region_type_id)
                else:
                    logger.warning('staining type not found: %s', rtype)

                # get the parent bands, and make them unique
                parents = list(self.make_parent_bands(band, set()))
                # alphabetical sort will put them in smallest to biggest
                parents.sort(reverse=True)

                # print("PARENTS of",maplocclass_id,"=",parents)
                # add the parents to the graph, in hierarchical order
                # TODO this is somewhat inefficient due to
                # re-adding upper-level nodes when iterating over the file
                # TODO PYLINT Consider using enumerate
                # instead of iterating with range and len
                for i in range(len(parents)):
                    pclassid = cclassid+parents[i]  # class chr parts
                    pclass_label = \
                        makeChromLabel(chrom+parents[i], genome_label)

                    rti = getChrPartTypeByNotation(parents[i])

                    model.addClassToGraph(pclassid, pclass_label, rti)

                    # for canonical chromosomes,
                    # then the subbands are subsequences of the full band
                    # add the subsequence stuff as restrictions
                    if i < len(parents) - 1:
                        pid = cclassid+parents[i+1]   # the instance
                        model.addOWLPropertyClassRestriction(
                            pclassid,
                            Feature.object_properties['is_subsequence_of'],
                            pid)
                        model.addOWLPropertyClassRestriction(
                            pid,
                            Feature.object_properties['has_subsequence'],
                            pclassid)

                    else:
                        # add the last one (p or q usually)
                        # as attached to the chromosome
                        model.addOWLPropertyClassRestriction(
                            pclassid,
                            Feature.object_properties['is_subsequence_of'],
                            cclassid)
                        model.addOWLPropertyClassRestriction(
                            cclassid,
                            Feature.object_properties['has_subsequence'],
                            pclassid)

                # connect the band here to the first one in the parent list
                if len(parents) > 0:
                    model.addOWLPropertyClassRestriction(
                        maplocclass_id,
                        Feature.object_properties['is_subsequence_of'],
                        cclassid+parents[0])
                    model.addOWLPropertyClassRestriction(
                        cclassid+parents[0],
                        Feature.object_properties['has_subsequence'],
                        maplocclass_id)

                if limit is not None and line_counter > limit:
                    break

        # TODO figure out the staining intensities for the encompassing bands

        return

    def make_parent_bands(self, band, child_bands):
        """
        this will determine the grouping bands that it belongs to, recursively
        13q21.31 ==>  13, 13q, 13q2, 13q21, 13q21.3, 13q21.31

        :param band:
        :param child_bands:
        :return:

        """
        m = re.match(r'([pq][A-H\d]+(?:\.\d+)?)', band)
        if len(band) > 0:
            if m:
                p = str(band[0:len(band)-1])
                p = re.sub(r'\.$', '', p)
                if p is not None:
                    child_bands.add(p)
                    self.make_parent_bands(p, child_bands)
        else:
            child_bands = set()
        return child_bands

    def map_type_of_region(self, regiontype):
        """
        Note that "stalk" refers to the short arm of acrocentric chromosomes
        chr13,14,15,21,22 for human.
        :param regiontype:
        :return:

        """
        so_id = Feature.types['chromosome_part']

        if regiontype in self.region_type_map.keys():
            so_id = self.region_type_map.get(regiontype)
        else:
            logger.warning(
                "Unmapped code %s. Defaulting to chr_part 'SO:0000830'.",
                regiontype)

        return so_id

    def _check_tax_ids(self):
        for taxon in self.tax_ids:
            if str(taxon) not in self.files:
                raise Exception("Taxon " + str(taxon) +
                                " not supported by source Monochrom")

    def getTestSuite(self):
        # import unittest
        # from tests.test_ucscbands import UCSCBandsTestCase
        test_suite = None
        # test_suite = \
        #   unittest.TestLoader().loadTestsFromTestCase(UCSCBandsTestCase)

        return test_suite


def getChrPartTypeByNotation(notation):
    """
    This method will figure out the kind of feature that a given band
    is based on pattern matching to standard karyotype notation.
    (e.g. 13q22.2 ==> chromosome sub-band)

    This has been validated against human, mouse, fish, and rat nomenclature.
    :param notation: the band (without the chromosome prefix)
    :return:

    """

    # Note that for mouse,
    # they don't usually include the "q" in their notation,
    # though UCSC does. We may need to adjust for that here

    if re.match(r'p$', notation):
        rti = Feature.types['short_chromosome_arm']
    elif re.match(r'q$', notation):
        rti = Feature.types['long_chromosome_arm']
    elif re.match(r'[pq][A-H\d]$', notation):
        rti = Feature.types['chromosome_region']
    elif re.match(r'[pq][A-H\d]\d', notation):
        rti = Feature.types['chromosome_band']
    elif re.match(r'[pq][A-H\d]\d\.\d+', notation):
        rti = Feature.types['chromosome_subband']
    else:
        rti = Feature.types['chromosome_part']

    return rti
