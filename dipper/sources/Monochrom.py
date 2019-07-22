import re
import gzip
import logging
from dipper.sources.Source import Source
from dipper.models.GenomicFeature import makeChromID, makeChromLabel
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model


LOG = logging.getLogger(__name__)
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
    and part.
    (note: in truth we create blank nodes and then pretend they are something else. TEC)

    We differentiate species by first creating a species-specific
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

    Since this is small, and we have limited other items in our test set to
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

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'monochrom',
            ingest_title='Monarch Chromosome Ontology',
            ingest_url='https://monarchinitiative.org',   # TODO can we be more specific
            license_url='http://creativecommons.org/licenses/by/4.0/'  # make our lic
            # data_rights=None,
            # file_handle=None
        )

        self.tax_ids = tax_ids
        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955, 10116, 9913, 9031, 9823, 9940, 9796]

        self.tax_ids = [str(x) for x in self.tax_ids]  # needed, if they are passed in

        self._check_tax_ids()

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)

    def parse(self, limit=None):

        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        for taxon in self.tax_ids:
            self._get_chrbands(limit, str(taxon))

        # using the full graph as the test here
        self.testgraph = self.graph
        LOG.info("Done parsing files.")

    def _get_chrbands(self, limit, taxon, genome_id=None):
        """
        For the given taxon, it will fetch the chr band file.
        We will not deal with the coordinate information with this parser.
        Here, we only are concerned with building the partonomy.
        :param limit:
        :param: taxon:
        :param: genome
        :return:

        """
        model = Model(self.graph)
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files[taxon]['file']))
        LOG.info("Processing Chr bands from FILE: %s", myfile)
        geno = Genotype(self.graph)

        # build the organism's genome from the taxon
        genome_label = self.files[taxon]['genome_label']
        taxon_id = 'NCBITaxon:' + taxon

        # add the taxon as a class.  adding the class label elsewhere
        model.addClassToGraph(taxon_id, None)
        model.addSynonym(taxon_id, genome_label)

        if genome_id is None:
            genome_id = geno.makeGenomeID(taxon_id)  # makes a blank node allways
        geno.addGenome(taxon_id, genome_label)
        model.addOWLPropertyClassRestriction(
            genome_id, self.globaltt['in taxon'], taxon_id)

        placed_scaffold_pattern = r'chr(\d+|X|Y|Z|W|MT|M)'
        # currently unused patterns
        # unlocalized_scaffold_pattern = placed_scaffold_pattern + r'_(\w+)_random'
        # unplaced_scaffold_pattern = r'chrUn_(\w+)'

        col = ['chrom', 'start', 'stop', 'band', 'rtype']
        with gzip.open(myfile, 'rb') as reader:
            for line in reader:
                line_counter += 1
                # skip comments
                line = line.decode().strip()
                if line[0] == '#':
                    continue
                # chr13	4500000	10000000	p12	stalk
                row = line.split('\t')
                chrom = row[col.index('chrom')]
                band = row[col.index('band')]
                rtype = row[col.index('rtype')]
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

                mch = re.match(placed_scaffold_pattern+r'$', chrom)
                if mch is not None and len(mch.groups()) == 1:
                    # the chromosome is the first match of the pattern
                    # chrom = m.group(1)  # TODO unused
                    pass
                else:
                    # let's skip over anything that isn't a placed_scaffold
                    LOG.info("Skipping non-placed chromosome %s", chrom)
                    continue
                # the chrom class, taxon as the reference
                cclassid = makeChromID(chrom, taxon, 'CHR')

                # add the chromosome as a class
                geno.addChromosomeClass(chrom, taxon_id, genome_label)
                model.addOWLPropertyClassRestriction(
                    cclassid, self.globaltt['member of'], genome_id)

                # add the band(region) as a class
                maplocclass_id = cclassid+band
                maplocclass_label = makeChromLabel(chrom+band, genome_label)
                if band is not None and band.strip() != '':
                    region_type_id = self.map_type_of_region(rtype)
                    model.addClassToGraph(
                        maplocclass_id, maplocclass_label,
                        region_type_id)
                else:
                    region_type_id = self.globaltt['chromosome']
                # add the staining intensity of the band
                if re.match(r'g(neg|pos|var)', rtype):
                    if region_type_id in [
                            self.globaltt['chromosome_band'],
                            self.globaltt['chromosome_subband']]:
                        stain_type = self.resolve(rtype)
                        if stain_type is not None:
                            model.addOWLPropertyClassRestriction(
                                maplocclass_id,
                                self.globaltt['has_sequence_attribute'],
                                self.resolve(rtype))
                    else:
                        # usually happens if it's a chromosome because
                        # they don't actually have banding info
                        LOG.info("feature type %s != chr band", region_type_id)
                else:
                    LOG.warning('staining type not found: %s', rtype)

                # get the parent bands, and make them unique
                parents = list(self.make_parent_bands(band, set()))
                # alphabetical sort will put them in smallest to biggest
                parents.sort(reverse=True)

                # print("PARENTS of", maplocclass_id, "=", parents)
                # add the parents to the graph, in hierarchical order
                # TODO this is somewhat inefficient due to
                # re-adding upper-level nodes when iterating over the file
                for prnt in parents:
                    parent = prnt.strip()
                    if parent is None or parent == "":
                        continue
                    pclassid = cclassid + parent  # class chr parts
                    pclass_label = makeChromLabel(chrom + parent, genome_label)
                    rti = getChrPartTypeByNotation(parent, self.graph)
                    model.addClassToGraph(pclassid, pclass_label, rti)

                    # for canonical chromosomes,
                    # then the subbands are subsequences of the full band
                    # add the subsequence stuff as restrictions

                    if prnt != parents[-1]:
                        grandparent = 1 + parents.index(prnt)
                        pid = cclassid + parents[grandparent]   # the instance
                        model.addOWLPropertyClassRestriction(
                            pclassid, self.globaltt['is subsequence of'], pid)
                        model.addOWLPropertyClassRestriction(
                            pid, self.globaltt['has subsequence'], pclassid)
                    else:
                        # add the last one (p or q usually)
                        # as attached to the chromosome
                        model.addOWLPropertyClassRestriction(
                            pclassid, self.globaltt['is subsequence of'], cclassid)
                        model.addOWLPropertyClassRestriction(
                            cclassid, self.globaltt['has subsequence'], pclassid)

                # connect the band here to the first one in the parent list
                if len(parents) > 0:
                    model.addOWLPropertyClassRestriction(
                        maplocclass_id, self.globaltt['is subsequence of'],
                        cclassid + parents[0])
                    model.addOWLPropertyClassRestriction(
                        cclassid + parents[0], self.globaltt['has subsequence'],
                        maplocclass_id)

                if limit is not None and line_counter > limit:
                    break

        # TODO figure out the staining intensities for the encompassing bands

    def make_parent_bands(self, band, child_bands):
        """
        this will determine the grouping bands that it belongs to, recursively
        13q21.31 ==>  13, 13q, 13q2, 13q21, 13q21.3, 13q21.31

        :param band:
        :param child_bands:
        :return:

        """
        mch = re.match(r'([pq][A-H\d]+(?:\.\d+)?)', band)
        if len(band) > 0:
            if mch:
                prnt = str(band[0:len(band)-1])
                prnt = re.sub(r'\.$', '', prnt)
                if prnt is not None:
                    child_bands.add(prnt)
                    self.make_parent_bands(prnt, child_bands)
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

        if regiontype in self.localtt:
            so_id = self.resolve(regiontype)
        else:
            so_id = self.globaltt['chromosome_part']
            LOG.warning(
                "Unmapped code %s. Defaulting to chr_part '" +
                self.globaltt['chromosome_part'] + "'.", regiontype)

        return so_id

    def _check_tax_ids(self):
        for taxon in self.tax_ids:
            if str(taxon) not in self.files:
                raise Exception(
                    "Taxon " + str(taxon) + " not supported by source Monochrom")

    def getTestSuite(self):
        # import unittest
        # from tests.test_ucscbands import UCSCBandsTestCase
        test_suite = None
        # test_suite = unittest.TestLoader().loadTestsFromTestCase(UCSCBandsTestCase)

        return test_suite


def getChrPartTypeByNotation(notation, graph=None):
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
        rti = graph.globaltt['short_chromosome_arm']
    elif re.match(r'q$', notation):
        rti = graph.globaltt['long_chromosome_arm']
    elif re.match(r'[pq][A-H\d]$', notation):
        rti = graph.globaltt['chromosome_region']
    elif re.match(r'[pq][A-H\d]\d', notation):
        rti = graph.globaltt['chromosome_band']
    elif re.match(r'[pq][A-H\d]\d\.\d+', notation):
        rti = graph.globaltt['chromosome_subband']
    else:
        rti = graph.globaltt['chromosome_part']

    return rti
