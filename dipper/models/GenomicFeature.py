import logging
import re
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.graph.Graph import Graph


logger = logging.getLogger(__name__)


class Feature():
    """
    Dealing with genomic features here.  By default they are all faldo:Regions.
    We use SO for typing genomic features. At the moment,
    RO:has_subsequence is the default relationship
    between the regions, but this should be tested/verified.

    TODO:
    the graph additions are in the addXToFeature functions,
    but should be separated.
    TODO:
    this will need to be extended to properly deal with
    fuzzy positions in faldo.

    """

    object_properties = {
        'location': 'faldo:location',
        'begin': 'faldo:begin',
        'end': 'faldo:end',
        'reference': 'faldo:reference',
        'gene_product_of': 'RO:0002204',
        'has_gene_product': 'RO:0002205',
        'is_about': 'IAO:0000136',
        'has_subsequence': 'RO:0002524',
        'is_subsequence_of': 'RO:0002525',
        'has_staining_intensity': 'GENO:0000207',
        'upstream_of_sequence_of': 'RO:0002528',
        'downstream_of_sequence_of': 'RO:0002529'
    }

    data_properties = {
        'position': 'faldo:position',
    }

    annotation_properties = {}

    properties = object_properties.copy()
    properties.update(data_properties)
    properties.update(annotation_properties)

    types = {
        'region': 'faldo:Region',
        'Position': 'faldo:Position',
        # big P for Position type.  little p for position property
        'FuzzyPosition': 'faldo:FuzzyPosition',
        'chromosome': 'SO:0000340',
        'chromosome_arm': 'SO:0000105',
        'chromosome_band': 'SO:0000341',
        'chromosome_part': 'SO:0000830',
        'long_chromosome_arm': 'GENO:0000629',
        'short_chromosome_arm': 'GENO:0000628',
        'chromosome_region': 'GENO:0000614',
        'chromosome_subband': 'GENO:0000616',
        'centromere': 'SO:0000577',
        'plus_strand': 'faldo:PlusStrandPosition',
        'minus_strand': 'faldo:MinusStrandPosition',
        'both_strand': 'faldo:BothStrandPosition',
        'score': 'SO:0001685',
        # FIXME - score is not a good solution, too generic
        'reference_genome': 'SO:0001505',
        'genome': 'SO:0001026',
        'assembly_component': 'SO:0000143',
        'SNP': 'SO:0000694',
        'haplotype': 'GENO:0000871',

        # the following are sequence attributes:
        'band_intensity':  'GENO:0000618',
        'gneg': 'GENO:0000620',
        'gpos': 'GENO:0000619',
        'gpos100': 'GENO:0000622',
        'gpos75': 'GENO:0000623',
        'gpos50': 'GENO:0000624',
        'gpos25': 'GENO:0000625',
        'gvar': 'GENO:0000621',
        'gpos33': 'GENO:0000633',
        'gpos66': 'GENO:0000632'
    }

    def __init__(self, graph, feature_id=None, label=None,
                 feature_type=None, description=None):
        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)
        self.id = feature_id
        self.label = label
        self.type = feature_type
        self.description = description
        self.start = None
        self.stop = None
        return

    def addFeatureStartLocation(
            self, coordinate, reference_id, strand=None,
            position_types=None):
        """
        Adds coordinate details for the start of this feature.
        :param coordinate:
        :param reference_id:
        :param strand:
        :param position_types:

        :return:

        """

        # make an object for the start, which has:
        # {coordinate : integer, reference : reference_id, types = []}
        self.start = self._getLocation(coordinate, reference_id, strand,
                                       position_types)

        return

    def addFeatureEndLocation(
            self, coordinate, reference_id, strand=None,
            position_types=None):
        """
        Adds the coordinate details for the end of this feature
        :param coordinate:
        :param reference_id:
        :param strand:

        :return:

        """

        self.stop = self._getLocation(coordinate, reference_id, strand,
                                      position_types)

        return

    def _getLocation(self, coordinate, reference_id, strand, position_types):
        """
        Make an object for the location, which has:
        {coordinate : integer, reference : reference_id, types = []}
        where the strand is indicated in the type array
        :param coordinate:
        :param reference_id:
        :param strand:
        :param position_types:

        :return:

        """

        loc = dict()
        loc['coordinate'] = coordinate
        loc['reference'] = reference_id
        loc['type'] = []
        strand_id = self._getStrandType(strand)
        if strand_id is not None:
            loc['type'].append(strand_id)
        if position_types is not None:
            loc['type'] += position_types
        if position_types == []:
            loc['type'].append(self.types['Position'])

        return loc

    def _getStrandType(self, strand):
        """
        :param strand:
        :return:
        """

        # TODO make this a dictionary/enum:  PLUS, MINUS, BOTH, UNKNOWN
        strand_id = None
        if strand == '+':
            strand_id = self.types['plus_strand']
        elif strand == '-':
            strand_id = self.types['minus_strand']
        elif strand == '.':
            strand_id = self.types['both_strand']
        elif strand is None:  # assume this is Unknown
            pass
        else:
            logger.warning("strand type could not be mapped: %s", str(strand))

        return strand_id

    def addFeatureToGraph(
            self, add_region=True, region_id=None,
            feature_as_class=False):
        """
        We make the assumption here that all features are instances.
        The features are located on a region,
        which begins and ends with faldo:Position
        The feature locations leverage the Faldo model,
        which has a general structure like:
        Triples:
        feature_id a feature_type (individual)
        faldo:location region_id
        region_id a faldo:region
        faldo:begin start_position
        faldo:end end_position
        start_position a
        (any of: faldo:(((Both|Plus|Minus)Strand)|Exact)Position)
        faldo:position Integer(numeric position)
        faldo:reference reference_id
        end_position a
        (any of: faldo:(((Both|Plus|Minus)Strand)|Exact)Position)
        faldo:position Integer(numeric position)
        faldo:reference reference_id

        :param graph:

        :return:

        """

        if feature_as_class:
            self.model.addClassToGraph(self.id, self.label, self.type,
                                       self.description)
        else:
            self.model.addIndividualToGraph(self.id, self.label, self.type,
                                            self.description)

        if self.start is None and self.stop is None:
            add_region = False

        if add_region:
            # create a region that has the begin/end positions
            regionchr = re.sub(r'\w+\:_?', '', self.start['reference'])
            if region_id is None:
                # in case the values are undefined
                # if we know only one of the coordinates,
                # then we'll add an "unknown" other.
                st = sp = 'UN'
                strand = None
                if self.start is not None and \
                        self.start['coordinate'] is not None:
                    st = str(self.start['coordinate'])
                    strand = self._getStrandStringFromPositionTypes(
                        self.start['type'])
                if self.stop is not None and\
                        self.stop['coordinate'] is not None:
                    sp = str(self.stop['coordinate'])
                    if strand is not None:
                        strand = self._getStrandStringFromPositionTypes(
                            self.stop['type'])
                # assume that the strand is the same for both start and stop.
                # this will need to be fixed in the future
                region_items = [regionchr, st, sp]
                if strand is not None:
                    region_items += [strand]
                region_id = '-'.join(region_items)
                rid = region_id
                rid = re.sub(r'\w+\:', '', rid, 1)  # replace the id prefix
                rid = '_:'+rid+"-Region"
                region_id = rid

            self.graph.addTriple(self.id, self.properties['location'],
                                 region_id)
            self.model.addIndividualToGraph(region_id, None, 'faldo:Region')
        else:
            region_id = self.id
            self.model.addType(region_id, 'faldo:Region')

        # add the start/end positions to the region
        beginp = endp = None
        if self.start is not None:
            beginp = self._makePositionId(self.start['reference'],
                                          self.start['coordinate'],
                                          self.start['type'])
            self.addPositionToGraph(self.start['reference'],
                                    self.start['coordinate'],
                                    self.start['type'])

        if self.stop is not None:
            endp = self._makePositionId(self.stop['reference'],
                                        self.stop['coordinate'],
                                        self.stop['type'])
            self.addPositionToGraph(self.stop['reference'],
                                    self.stop['coordinate'],
                                    self.stop['type'])

        self.addRegionPositionToGraph(region_id, beginp, endp)

        # {coordinate : integer, reference : reference_id, types = []}

        return

    def _getStrandStringFromPositionTypes(self, tylist):
        strand = None
        if self.types['plus_strand'] in tylist:
            strand = 'plus'
        elif self.types['minus_strand'] in tylist:
            strand = 'minus'
        elif self.types['both_strand'] in tylist:
            strand = 'both'
        else:
            strand = None  # it is stranded, but we don't know what it is

        return strand

    def _makePositionId(self, reference, coordinate, types=None):
        """
        Note that positions should have a reference (we will enforce).
        Only exact positions need a coordinate.
        :param reference:
        :param coordinate:
        :param types:
        :return:
        """

        if reference is None:
            logger.error("Trying to make position with no reference.")
            return None

        curie = '_:'
        reference = re.sub(r'\w+\:', '', reference, 1)
        if re.match(r'^_', reference):
            # this is in the case if the reference is a bnode
            reference = re.sub(r'^_', '', reference)
        curie += reference
        if coordinate is not None:
            # just in case it isn't a string already
            curie = '-'.join((curie, str(coordinate)))
        if types is not None:
            tstring = self._getStrandStringFromPositionTypes(types)
            if tstring is not None:
                curie = '-'.join((curie, tstring))

        return curie

    def addRegionPositionToGraph(
            self, region_id, begin_position_id,
            end_position_id):

        if begin_position_id is None:
            pass
            # logger.warn(
            #   "No begin position specified for region %s", region_id)
        else:
            self.graph.addTriple(region_id, self.properties['begin'],
                                 begin_position_id)

        if end_position_id is None:
            pass
            # logger.warn("No end position specified for region %s", region_id)
        else:
            self.graph.addTriple(region_id, self.properties['end'],
                                 end_position_id)

        return

    def addPositionToGraph(
            self, reference_id, position,
            position_types=None, strand=None):
        """
        Add the positional information to the graph, following the faldo model.
        We assume that if the strand is None,
        we give it a generic "Position" only.
        Triples:
        my_position a (any of: faldo:(((Both|Plus|Minus)Strand)|Exact)Position)
        faldo:position Integer(numeric position)
        faldo:reference reference_id

        :param graph:
        :param reference_id:
        :param position:
        :param position_types:
        :param strand:

        :return:  Identifier of the position created

        """
        pos_id = self._makePositionId(reference_id, position, position_types)
        if position is not None:
            self.graph.addTriple(pos_id, self.properties['position'],
                                 position, object_is_literal=True,
                                 literal_type="xsd:integer")
        self.graph.addTriple(
            pos_id, self.properties['reference'], reference_id)
        if position_types is not None:
            for pos_type in position_types:
                self.model.addType(pos_id, pos_type)
        s = None
        if strand is not None:
            s = strand
            if not re.match(r'faldo', strand):
                # not already mapped to faldo, so expect we need to map it
                s = self._getStrandType(strand)
        # else:
        #    s = self.types['both_strand']
        if s is None and (position_types is None or position_types == []):
            s = self.types['Position']

        if s is not None:
            self.model.addType(pos_id, s)

        return pos_id

    def addSubsequenceOfFeature(self, parentid):
        """
        This will add reciprocal triples like:
        feature is_subsequence_of parent
        parent has_subsequence feature
        :param graph:
        :param parentid:

        :return:

        """
        self.graph.addTriple(
            self.id, self.properties['is_subsequence_of'], parentid)
        self.graph.addTriple(
            parentid, self.properties['has_subsequence'], self.id)

        return

    def addTaxonToFeature(self, taxonid):
        """
        Given the taxon id, this will add the following triple:
        feature in_taxon taxonid
        :param graph:
        :param taxonid:
        :return:
        """
        # TEC: should taxon be set in __init__()?
        self.taxon = taxonid
        self.graph.addTriple(
            self.id, Assoc.properties['in_taxon'], self.taxon)

        return

    def addFeatureProperty(self, property_type, property):
        self.graph.addTriple(self.id, property_type, property)
        return


def makeChromID(chrom, reference=None, prefix=None):
    """
    This will take a chromosome number and a NCBI taxon number,
    and create a unique identifier for the chromosome.  These identifiers
    are made in the @base space like:
    Homo sapiens (9606) chr1 ==> :9606chr1
    Mus musculus (10090) chrX ==> :10090chrX

    :param chrom: the chromosome (preferably without any chr prefix)
    :param reference: the numeric portion of the taxon id

    :return:

    """
    if reference is None:
        logger.warning(
            'No reference for this chr. You may have conflicting ids')

    # replace any chr-like prefixes with blank to standardize
    c = re.sub(r'ch(r?)[omse]*', '', str(chrom))

    # remove the build/taxon prefixes to look cleaner
    r = reference
    if re.match(r'.+:.+', reference):
        r = re.match(r'.*:(.*)', reference).group(1)
    id = ''.join((r, 'chr', c))
    if prefix is None:
        id = ''.join(('_', id))
    else:
        id = ':'.join((prefix, id))

    return id


def makeChromLabel(chrom, reference=None):
    c = re.sub(r'ch(r?)[omse\.]*', '', str(chrom))
    c = 'chr'+c
    if reference is None:
        label = c
    else:
        label = c+' ('+reference+')'

    return label
