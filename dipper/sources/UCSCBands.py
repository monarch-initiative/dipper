import re
import gzip
import logging

from dipper.sources.Source import Source
from dipper.sources.Monochrom import Monochrom, getChrPartTypeByNotation
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.models.Dataset import Dataset
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model


logger = logging.getLogger(__name__)


class UCSCBands(Source):
    '''
    This will take the UCSC defintions of cytogenic bands and create the nested
    structures to enable overlap and containment queries. We use
    ```Monochrom.py``` to create the OWL-classes of the chromosomal parts.
    Here, we simply worry about the instance-level values for particular genome
    builds.

    Given a chr band definition, the nested containment structures look like:
    13q21.31 ==>  13q21.31,  13q21.3,  13q21,  13q2,  13q, 13

    We determine the containing regions of the band by parsing the band-string;
    since each alphanumeric is a significant "place", we can split it
    with the shorter strings being parents of the longer string.
    #
    Here we create build-specific chroms, which are instances of the classes
    produced from ```Monochrom.py```.
    You can instantiate any number of builds for a genome.

    We leverage the Faldo model here for region definitions,
    and map each of the chromosomal parts to SO.

    We differentiate the build by adding the build id to the identifier prior
    to the chromosome number.
    These then are instances of the species-specific chromosomal class.

    The build-specific chromosomes are created like:
    <pre>
    <build number>chr<num><band>
    with triples for a given band like:
    _:hg19chr1p36.33
        rdfs:type SO:chromosome_band,
        faldo:Region, CHR:9606chr1p36.33,
        subsequence_of _:hg19chr1p36.3,
        faldo:location
        [ a faldo:BothStrandPosition
                faldo:begin 0,
                faldo:end 2300000,
                faldo:reference 'hg19'
        ] .
    </pre>
    where any band in the file is an instance of a chr_band
    (or a more specific type), is a subsequence of it's containing region,
    and is located in the specified coordinates.

    We do not have a separate graph for testing.

    TODO: any species by commandline argument

    '''

    HGGP = 'http://hgdownload.cse.ucsc.edu/goldenPath'
    files = {
        # TODO accommodate multiple builds per species
        '9606': {
            'file': 'hg19cytoBand.txt.gz',
            'url': HGGP + '/hg19/database/cytoBand.txt.gz',
            'build_num': 'hg19',
            'genome_label': 'Human'
        },
        '10090': {
            'file': 'mm10cytoBand.txt.gz',
            'url': HGGP + '/mm10/database/cytoBandIdeo.txt.gz',
            'build_num': 'mm10',
            'genome_label': 'Mouse'
        },
        # Note that there are no bands,
        # arms or staining components for the species below at the moment
        '7955': {
            'file': 'danRer10cytoBand.txt.gz',
            'url': HGGP + '/danRer10/database/cytoBandIdeo.txt.gz',
            'build_num': 'danRer10',
            'genome_label': 'Zebrafish'
        },
        '9913': {
            'file': 'bosTau7cytoBand.txt.gz',
            'url': HGGP + '/bosTau7/database/cytoBandIdeo.txt.gz',
            'build_num': 'bosTau7',
            'genome_label': 'cow'
        },
        '9031': {
            'file': 'galGal4cytoBand.txt.gz',
            'url': HGGP + '/galGal4/database/cytoBandIdeo.txt.gz',
            'build_num': 'galGal4',
            'genome_label': 'chicken'
        },
        '9823': {
            'file': 'susScr3cytoBand.txt.gz',
            'url': HGGP + '/susScr3/database/cytoBandIdeo.txt.gz',
            'build_num': 'susScr3',
            'genome_label': 'pig'
        },
        '9940': {
            'file': 'oviAri3cytoBand.txt.gz',
            'url': HGGP + '/oviAri3/database/cytoBandIdeo.txt.gz',
            'build_num': 'oviAri3',
            'genome_label': 'sheep'
        },
        '9796': {
            'file': 'equCab2cytoBand.txt.gz',
            'url': HGGP + '/equCab2/database/cytoBandIdeo.txt.gz',
            'build_num': 'equCab2',
            'genome_label': 'horse'
        },
        # TODO rainbow trout, 8022, when available
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'ucscbands')

        self.tax_ids = tax_ids

        # Defaults
        if self.tax_ids is None:
            # self.tax_ids = [9606, 10090, 7955]
            self.tax_ids = [9606, 10090, 7955, 9913, 9031, 9823, 9940, 9796]

        # TODO add other species as defaults

        self._check_tax_ids()

        self.dataset = Dataset('ucscbands', 'UCSC Cytogenic Bands',
                               'http://hgdownload.cse.ucsc.edu', None,
                               'http://genome.ucsc.edu/license/')

        # data-source specific warnings
        # (will be removed when issues are cleared)

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

        self._create_genome_builds()

        # using the full graph as the test here
        self.testgraph = self.graph
        logger.info("Done parsing files.")

        return

    def _get_chrbands(self, limit, taxon):
        """
        :param limit:
        :return:

        """
        model = Model(self.graph)
        # TODO PYLINT figure out what limit was for and why it is unused
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files[taxon]['file']))
        logger.info("Processing Chr bands from FILE: %s", myfile)
        geno = Genotype(self.graph)
        monochrom = Monochrom(self.graph_type, self.are_bnodes_skized)

        # used to hold band definitions for a chr
        # in order to compute extent of encompasing bands

        mybands = {}
        # build the organism's genome from the taxon
        genome_label = self.files[taxon]['genome_label']
        taxon_id = 'NCBITaxon:'+taxon

        # add the taxon as a class.  adding the class label elsewhere
        model.addClassToGraph(taxon_id, None)
        model.addSynonym(taxon_id, genome_label)

        geno.addGenome(taxon_id, genome_label)

        # add the build and the taxon it's in
        build_num = self.files[taxon]['build_num']
        build_id = 'UCSC:'+build_num
        geno.addReferenceGenome(build_id, build_num, taxon_id)

        # process the bands
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                # skip comments
                line = line.decode().strip()
                if re.match('^#', line):
                    continue

                # chr13	4500000	10000000	p12	stalk
                (scaffold, start, stop, band_num, rtype) = line.split('\t')
                line_counter += 1

                # NOTE some less-finished genomes have
                # placed and unplaced scaffolds
                # * Placed scaffolds:
                #       the scaffolds have been placed within a chromosome.
                # * Unlocalized scaffolds:
                #   although the chromosome within which the scaffold occurs
                #   is known, the scaffold's position or orientation
                #   is not known.
                # * Unplaced scaffolds:
                #   it is not known which chromosome the scaffold belongs to
                #
                # find out if the thing is a full on chromosome, or a scaffold:
                # ex: unlocalized scaffold: chr10_KL568008v1_random
                # ex: unplaced scaffold: chrUn_AABR07022428v1
                placed_scaffold_pattern = r'(chr(?:\d+|X|Y|Z|W|M))'
                unlocalized_scaffold_pattern = \
                    placed_scaffold_pattern+r'_(\w+)_random'
                unplaced_scaffold_pattern = r'chr(Un(?:_\w+)?)'

                m = re.match(placed_scaffold_pattern+r'$', scaffold)
                if m is not None and len(m.groups()) == 1:
                    # the chromosome is the first match of the pattern
                    chrom_num = m.group(1)
                else:
                    # skip over anything that isn't a placed_scaffold
                    # at the class level
                    logger.info("Found non-placed chromosome %s", scaffold)
                    chrom_num = None

                m_chr_unloc = re.match(unlocalized_scaffold_pattern, scaffold)
                m_chr_unplaced = re.match(unplaced_scaffold_pattern, scaffold)

                scaffold_num = None
                if m:
                    pass
                elif m_chr_unloc is not None and\
                        len(m_chr_unloc.groups()) == 2:
                    chrom_num = m_chr_unloc.group(1)
                    scaffold_num = chrom_num+'_'+m_chr_unloc.group(2)
                elif m_chr_unplaced is not None and\
                        len(m_chr_unplaced.groups()) == 1:
                    scaffold_num = m_chr_unplaced.group(1)
                else:
                    logger.error(
                        "There's a chr pattern that we aren't matching: %s",
                        scaffold)

                if chrom_num is not None:
                    # the chrom class (generic) id
                    chrom_class_id = makeChromID(chrom_num, taxon, 'CHR')

                    # first, add the chromosome class (in the taxon)
                    geno.addChromosomeClass(
                        chrom_num, taxon_id, self.files[taxon]['genome_label'])

                    # then, add the chromosome instance (from the given build)
                    geno.addChromosomeInstance(chrom_num, build_id, build_num,
                                               chrom_class_id)

                    # add the chr to the hashmap of coordinates for this build
                    # the chromosome coordinate space is itself
                    if chrom_num not in mybands.keys():
                        mybands[chrom_num] = {
                            'min': 0,
                            'max': int(stop),
                            'chr': chrom_num,
                            'ref': build_id,
                            'parent': None,
                            'stain': None,
                            'type': Feature.types['chromosome']}

                if scaffold_num is not None:
                    # this will put the coordinates of the scaffold
                    # in the scaffold-space and make sure that the scaffold
                    # is part of the correct parent.
                    # if chrom_num is None,
                    # then it will attach it to the genome,
                    # just like a reg chrom
                    mybands[scaffold_num] = {
                        'min': start,
                        'max': stop,
                        'chr': scaffold_num,
                        'ref': build_id,
                        'parent': chrom_num,
                        'stain': None,
                        'type': Feature.types['assembly_component'],
                        'synonym': scaffold}

                if band_num is not None and band_num.strip() != '':
                    # add the specific band
                    mybands[chrom_num+band_num] = {'min': start,
                                                   'max': stop,
                                                   'chr': chrom_num,
                                                   'ref': build_id,
                                                   'parent': None,
                                                   'stain': None,
                                                   'type': None}

                    # add the staining intensity of the band
                    if re.match(r'g(neg|pos|var)', rtype):
                        mybands[chrom_num+band_num]['stain'] = \
                            Feature.types.get(rtype)

                    # get the parent bands, and make them unique
                    parents = list(
                        monochrom.make_parent_bands(band_num, set()))
                    # alphabetical sort will put them in smallest to biggest,
                    # so we reverse
                    parents.sort(reverse=True)
                    # print('parents of',chrom,band,':',parents)

                    if len(parents) > 0:
                        mybands[chrom_num+band_num]['parent'] = \
                            chrom_num+parents[0]
                else:
                    # TODO PYLINT why is 'parent'
                    # a list() a couple of lines up and a set() here?
                    parents = set()

                # loop through the parents and add them to the hash
                # add the parents to the graph, in hierarchical order
                # TODO PYLINT Consider using enumerate
                # instead of iterating with range and len
                for i in range(len(parents)):
                    rti = getChrPartTypeByNotation(parents[i])

                    pnum = chrom_num+parents[i]
                    sta = int(start)
                    sto = int(stop)
                    if pnum not in mybands.keys():
                        # add the parental band to the hash
                        b = {'min': min(sta, sto),
                             'max': max(sta, sto),
                             'chr': chrom_num,
                             'ref': build_id,
                             'parent': None,
                             'stain': None,
                             'type': rti}
                        mybands[pnum] = b
                    else:
                        # band already in the hash means it's a grouping band
                        # need to update the min/max coords
                        b = mybands.get(pnum)
                        b['min'] = min(sta, sto, b['min'])
                        b['max'] = max(sta, sto, b['max'])
                        mybands[pnum] = b

                        # also, set the max for the chrom
                        c = mybands.get(chrom_num)
                        c['max'] = max(sta, sto, c['max'])
                        mybands[chrom_num] = c

                    # add the parent relationships to each
                    if i < len(parents) - 1:
                        mybands[pnum]['parent'] = chrom_num+parents[i+1]
                    else:
                        # add the last one (p or q usually)
                        # as attached to the chromosome
                        mybands[pnum]['parent'] = chrom_num

        f.close()  # end looping through file

        # loop through the hash and add the bands to the graph
        for b in mybands.keys():
            myband = mybands.get(b)
            band_class_id = makeChromID(b, taxon, 'CHR')
            band_class_label = makeChromLabel(b, genome_label)
            band_build_id = makeChromID(b, build_num, 'MONARCH')
            band_build_label = makeChromLabel(b, build_num)
            # the build-specific chrom
            chrom_in_build_id = makeChromID(
                myband['chr'], build_num, 'MONARCH')
            # if it's != part, then add the class
            if myband['type'] != Feature.types['assembly_component']:
                model.addClassToGraph(band_class_id,
                                      band_class_label, myband['type'])
                bfeature = Feature(self.graph, band_build_id, band_build_label,
                                   band_class_id)
            else:
                bfeature = Feature(self.graph, band_build_id, band_build_label,
                                   myband['type'])
                if 'synonym' in myband:
                    model.addSynonym(band_build_id, myband['synonym'])

            if myband['parent'] is None:
                if myband['type'] == Feature.types['assembly_component']:
                    # since we likely don't know the chr,
                    # add it as a part of the build
                    geno.addParts(band_build_id, build_id)
            elif myband['type'] == Feature.types['assembly_component']:
                # geno.addParts(band_build_id, chrom_in_build_id)
                parent_chrom_in_build = makeChromID(myband['parent'],
                                                    build_num, 'MONARCH')
                bfeature.addSubsequenceOfFeature(parent_chrom_in_build)

            # add the band as a feature
            # (which also instantiates the owl:Individual)
            bfeature.addFeatureStartLocation(myband['min'], chrom_in_build_id)
            bfeature.addFeatureEndLocation(myband['max'], chrom_in_build_id)
            if 'stain' in myband and myband['stain'] is not None:
                # TODO 'has_staining_intensity' being dropped by MB
                bfeature.addFeatureProperty(
                    Feature.properties['has_staining_intensity'],
                    myband['stain'])

            # type the band as a faldo:Region directly (add_region=False)
            # bfeature.setNoBNodes(self.nobnodes)
            # to come when we merge in ZFIN.py
            bfeature.addFeatureToGraph(False)

        return

    def _create_genome_builds(self):
        """
        Various resources will map variations to either UCSC (hg*)
        or to NCBI assemblies. Here we create the equivalences between them.
        Data taken from:
        https://genome.ucsc.edu/FAQ/FAQreleases.html#release1

        :return:

        """

        # TODO add more species
        ucsc_assembly_id_map = {
            "9606": {
                "UCSC:hg38": "NCBIGenome:GRCh38",
                "UCSC:hg19": "NCBIGenome:GRCh37",
                "UCSC:hg18": "NCBIGenome:36.1",
                "UCSC:hg17": "NCBIGenome:35",
                "UCSC:hg16": "NCBIGenome:34",
                "UCSC:hg15": "NCBIGenome:33",
                },
            "7955": {
                "UCSC:danRer10": "NCBIGenome:GRCz10",
                "UCSC:danRer7":	"NCBIGenome:Zv9",
                "UCSC:danRer6": "NCBIGenome:Zv8",
                },
            "10090": {
                "UCSC:mm10": "NCBIGenome:GRCm38",
                "UCSC:mm9":	"NCBIGenome:37"
            },
            "9031": {
                "UCSC:galGal4": "NCBIAssembly:317958",
                },
            "9913": {
                "UCSC:bosTau7": "NCBIAssembly:GCF_000003205.5",
                },
            "9823": {
                "UCSC:susScr3": "NCBIAssembly:304498",
                },
            "9940": {
                "UCSC:oviAri3": "NCBIAssembly:GCF_000298735.1",
                },
            "9796": {
                "UCSC:equCab2": "NCBIAssembly:GCF_000002305.2",
                }
        }
        g = self.graph
        geno = Genotype(g)
        model = Model(g)
        logger.info("Adding equivalent assembly identifiers")
        for sp in ucsc_assembly_id_map:
            tax_num = sp
            tax_id = 'NCBITaxon:'+tax_num
            mappings = ucsc_assembly_id_map[sp]
            for i in mappings:
                ucsc_id = i
                ucsc_label = re.split(':', i)[1]
                mapped_id = mappings[i]
                mapped_label = re.split(':', mapped_id)[1]
                mapped_label = 'NCBI build '+str(mapped_label)
                geno.addReferenceGenome(ucsc_id, ucsc_label, tax_id)
                geno.addReferenceGenome(mapped_id, mapped_label, tax_id)
                model.addSameIndividual(ucsc_id, mapped_id)

        return

    def _check_tax_ids(self):
        for taxon in self.tax_ids:
            if str(taxon) not in self.files:
                raise Exception("Taxon " + str(taxon) + " not supported"
                                " by source UCSCBands")

    def getTestSuite(self):
        import unittest
        from tests.test_ucscbands import UCSCBandsTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(
            UCSCBandsTestCase)

        return test_suite
