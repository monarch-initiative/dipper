import logging
import re
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

    def __init__(
            self, graph, feature_id=None, label=None,
            feature_type=None, description=None):

        if isinstance(graph, Graph):
            self.graph = graph
        else:
            raise ValueError("{} is not a graph".graph)
        self.model = Model(self.graph)
        self.globaltt = self.graph.globaltt
        self.globaltcid = self.graph.globaltcid
        self.curie_map = self.graph.curie_map
        self.fid = feature_id
        self.label = label
        self.ftype = feature_type
        self.description = description
        self.start = None
        self.stop = None
        return

    def addFeatureStartLocation(
            self, coordinate, reference_id, strand=None, position_types=None):
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
        self.start = self._getLocation(coordinate, reference_id, strand, position_types)

        return

    def addFeatureEndLocation(
            self, coordinate, reference_id, strand=None, position_types=None):
        """
        Adds the coordinate details for the end of this feature
        :param coordinate:
        :param reference_id:
        :param strand:

        :return:

        """

        self.stop = self._getLocation(coordinate, reference_id, strand, position_types)

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

        loc = {}
        loc['coordinate'] = coordinate
        loc['reference'] = reference_id
        loc['type'] = []
        strand_id = self._getStrandType(strand)
        if strand_id is not None:
            loc['type'].append(strand_id)
        if position_types is not None:
            loc['type'] += position_types
        if position_types == []:
            loc['type'].append(self.globaltt['Position'])

        return loc

    def _getStrandType(self, strand):
        """
        :param strand:
        :return:
        """

        # TODO make this a dictionary/enum:  PLUS, MINUS, BOTH, UNKNOWN
        strand_id = None
        if strand == '+':
            strand_id = self.globaltt['plus_strand']
        elif strand == '-':
            strand_id = self.globaltt['minus_strand']
        elif strand == '.':
            strand_id = self.globaltt['both_strand']
        elif strand is None:  # assume this is Unknown
            pass
        else:
            logger.warning("strand type could not be mapped: %s", str(strand))

        return strand_id

    def addFeatureToGraph(
            self, add_region=True, region_id=None, feature_as_class=False):
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
            self.model.addClassToGraph(
                self.fid, self.label, self.ftype, self.description)
        else:
            self.model.addIndividualToGraph(
                self.fid, self.label, self.ftype, self.description)

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
                if self.start is not None and self.start['coordinate'] is not None:
                    st = str(self.start['coordinate'])
                    strand = self._getStrandStringFromPositionTypes(self.start['type'])
                if self.stop is not None and self.stop['coordinate'] is not None:
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

            self.graph.addTriple(self.fid, self.globaltt['location'], region_id)
            self.model.addIndividualToGraph(region_id, None, self.globaltt['Region'])
        else:
            region_id = self.fid
            self.model.addType(region_id, self.globaltt['region'])

        # add the start/end positions to the region
        beginp = endp = None
        if self.start is not None:
            beginp = self._makePositionId(
                self.start['reference'], self.start['coordinate'], self.start['type'])
            self.addPositionToGraph(
                self.start['reference'], self.start['coordinate'], self.start['type'])

        if self.stop is not None:
            endp = self._makePositionId(
                self.stop['reference'], self.stop['coordinate'], self.stop['type'])
            self.addPositionToGraph(
                self.stop['reference'], self.stop['coordinate'], self.stop['type'])

        self.addRegionPositionToGraph(region_id, beginp, endp)

        # {coordinate : integer, reference : reference_id, types = []}

        return

    def _getStrandStringFromPositionTypes(self, tylist):
        strand = None
        if self.globaltt['plus_strand'] in tylist:
            strand = 'plus'
        elif self.globaltt['minus_strand'] in tylist:
            strand = 'minus'
        elif self.globaltt['both_strand'] in tylist:
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

    def addRegionPositionToGraph(self, region_id, begin_position_id, end_position_id):

        if begin_position_id is None:
            pass
            # logger.warn("No begin position specified for region %s", region_id)
        else:
            self.graph.addTriple(region_id, self.globaltt['begin'], begin_position_id)

        if end_position_id is None:
            pass
            # logger.warn("No end position specified for region %s", region_id)
        else:
            self.graph.addTriple(region_id, self.globaltt['end'], end_position_id)

        return

    def addPositionToGraph(
            self, reference_id, position, position_types=None, strand=None):
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
            self.graph.addTriple(
                pos_id, self.globaltt['position'], position, object_is_literal=True,
                literal_type="xsd:integer")
        self.graph.addTriple(pos_id, self.globaltt['reference'], reference_id)
        if position_types is not None:
            for pos_type in position_types:
                self.model.addType(pos_id, pos_type)
        strnd = None
        if strand is not None:
            strnd = strand
            if not re.match(r'faldo', strand):
                # not already mapped to faldo, so expect we need to map it
                strnd = self._getStrandType(strand)
        # else:
        #    strnd = self.globaltt['both_strand']
        if strnd is None and (position_types is None or position_types == []):
            strnd = self.globaltt['Position']

        if strnd is not None:
            self.model.addType(pos_id, strnd)

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
        self.graph.addTriple(self.fid, self.globaltt['is subsequence of'], parentid)
        self.graph.addTriple(parentid, self.globaltt['has subsequence'], self.fid)

        return

    def addTaxonToFeature(self, taxonid):
        """
        Given the taxon id, this will add the following triple:
        feature in_taxon taxonid
        :param graph:
        :param taxonid:
        :return:
        """
        self.taxon = taxonid
        self.graph.addTriple(self.fid, self.globaltt['in taxon'], self.taxon)

        return

    def addFeatureProperty(self, property_type, property):
        self.graph.addTriple(self.fid, property_type, property)
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
    # blank nodes
    if reference is None:
        logger.warning('No reference for this chr. You may have conflicting ids')

    # replace any chr-like prefixes with blank to standardize
    chrid = re.sub(r'ch(r?)[omse]*', '', str(chrom))

    # remove the build/taxon prefixes to look cleaner
    ref = reference
    if re.match(r'.+:.+', reference):
        ref = re.match(r'.*:(.*)', reference).group(1)
    chrid = ''.join((ref, 'chr', chrid))
    if prefix is None:
        chrid = ''.join(('_', chrid))
    else:
        chrid = ':'.join((prefix, chrid))

    return chrid


def makeChromLabel(chrom, reference=None):
    chrm = re.sub(r'ch(r?)[omse\.]*', '', str(chrom))
    chrm = 'chr' + chrm
    if reference is None:
        label = chrm
    else:
        label = chrm + ' (' + reference + ')'

    return label
