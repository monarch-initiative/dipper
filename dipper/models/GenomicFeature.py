__author__ = 'nlw'

import logging
import re

from rdflib import Literal
from rdflib.namespace import RDF, XSD

from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.assoc.Association import Assoc


logger = logging.getLogger(__name__)


class Feature():
    """
    Dealing with genomic features here.  By default they are all faldo:Regions.
    We use SO for typing genomic features.
    At the moment, RO:has_subsequence is the default relationship between the regions, but this
    should be tested/verified.

    TODO: the graph additions are in the addXToFeature functions, but should be separated.
    TODO: this will need to be extended to properly deal with fuzzy positions in faldo.
    """

    object_properties = {
        'location': 'faldo:location',
        'begin': 'faldo:begin',
        'end': 'faldo:end',
        'reference': 'faldo:reference',
        'gene_product_of': 'RO:0002204',
        'has_gene_product': 'RO:0002205',
        'is_about': 'IAO:00000136',
        'has_subsequence': 'RO:0002524',
        'is_subsequence_of': 'RO:0002525',
        'has_staining_intensity': 'GENO:0000626',
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
        'Position': 'faldo:Position',  # big P for Position type.  little p for position property
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
        'score': 'SO:0001685',  # FIXME - this is not a good solution, too generic
        'reference_genome': 'SO:0001505',
        'genome': 'SO:0001026',
        'assembly_component': 'SO:0000143',
        'SNP': 'SO:0000694',

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

    def __init__(self, id, label, type, description=None):
        self.id = id
        self.label = label
        self.type = type
        self.description = description
        self.gu = GraphUtils(curie_map.get())
        self.start = None
        self.stop = None
        self.nobnodes = True  # TODO remove this before official release
        return

    def addFeatureStartLocation(self, coordinate, reference_id, strand=None, position_types=None):
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

    def addFeatureEndLocation(self, coordinate, reference_id, strand=None, position_types=None):
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
            logger.warn("strand type could not be mapped: %s", str(strand))

        return strand_id

    def addFeatureToGraph(self, graph, add_region=True, region_id=None, feature_as_class=False):
        """
        We make the assumption here that all features are instances.
        The features are located on a region, which begins and ends with faldo:Position
        The feature locations leverage the Faldo model, which has a general structure like:
        Triples:
        feature_id a feature_type (individual)
            faldo:location region_id
        region_id a faldo:region
            faldo:begin start_position
            faldo:end end_position
        start_position a (any of: faldo:(((Both|Plus|Minus)Strand)|Exact)Position)
            faldo:position Integer(numeric position)
            faldo:reference reference_id
        end_position a (any of: faldo:(((Both|Plus|Minus)Strand)|Exact)Position)
            faldo:position Integer(numeric position)
            faldo:reference reference_id

        :param graph:
        :return:
        """
        if feature_as_class:
            self.gu.addClassToGraph(graph, self.id, self.label, self.type, self.description)
        else:
            self.gu.addIndividualToGraph(graph, self.id, self.label, self.type, self.description)

        if self.start is None and self.stop is None:
            add_region = False

        if add_region:
            # create a region that has the begin/end positions
            regionchr = re.sub('\w+\:_?', '', self.start['reference'])
            if region_id is None:
                # in case the values are undefined
                # if we know only one of the coordinates, then we'll add an "unknown" other.
                st = sp = 'UN'
                strand = None
                if self.start is not None and self.start['coordinate'] is not None:
                    st = str(self.start['coordinate'])
                    strand = self._getStrandStringFromPositionTypes(self.start['type'])
                if self.stop is not None and self.stop['coordinate'] is not None:
                    sp = str(self.stop['coordinate'])
                    if strand is not None:
                        strand = self._getStrandStringFromPositionTypes(self.stop['type'])
                # assume that the strand is the same for both start and stop.  this will need to be fixed in the future
                region_items = [regionchr, st, sp]
                if strand is not None:
                    region_items += [strand]
                region_id = '-'.join(region_items)
                rid = region_id
                rid = re.sub('\w+\:', '', rid, 1)  # replace the id prefix
                rid = '_'+rid+"-Region"
                region_id = rid
                if self.nobnodes:
                    region_id = ':'+region_id
            self.gu.addTriple(graph, self.id, self.properties['location'], region_id)
            self.gu.addIndividualToGraph(graph, region_id, None, 'faldo:Region')
        else:
            region_id = self.id
            self.gu.addType(graph, region_id, 'faldo:Region')

        # add the start/end positions to the region
        beginp = endp = None
        if self.start is not None:
            beginp = self._makePositionId(self.start['reference'], self.start['coordinate'], self.start['type'])
            self.addPositionToGraph(graph, self.start['reference'], self.start['coordinate'], self.start['type'])

        if self.stop is not None:
            endp = self._makePositionId(self.stop['reference'], self.stop['coordinate'], self.stop['type'])
            self.addPositionToGraph(graph, self.stop['reference'], self.stop['coordinate'], self.stop['type'])

        self.addRegionPositionToGraph(graph, region_id, beginp, endp)

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
        Note that positions should have a reference (we will enforce).  Only exact positions need a coordinate.
        :param reference:
        :param coordinate:
        :param types:
        :return:
        """
        if reference is None:
            logger.error("Trying to make position with no reference.")
            return None

        i = '_'
        if self.nobnodes:
            i = ':'+i
        reference = re.sub('\w+\:', '', reference, 1)
        if re.match('^_', reference):
            reference = re.sub('^_','',reference)  # this is in the case if the reference is a bnode
        i += reference
        if coordinate is not None:
            i = '-'.join((i, str(coordinate)))      # just in case it isn't a string already
        if types is not None:
            tstring = self._getStrandStringFromPositionTypes(types)
            if tstring is not None:
                i = '-'.join((i, tstring))

        return i

    def addRegionPositionToGraph(self, graph, region_id, begin_position_id, end_position_id):

        if begin_position_id is None:
            pass
            # logger.warn("No begin position specified for region %s", region_id)
        else:
            self.gu.addTriple(graph, region_id, self.properties['begin'], begin_position_id)

        if end_position_id is None:
            pass
            # logger.warn("No end position specified for region %s", region_id)
        else:
            self.gu.addTriple(graph, region_id, self.properties['end'], end_position_id)

        return

    def addPositionToGraph(self, graph, reference_id, position, position_types=None, strand=None):
        """
        Add the positional information to the graph, following the faldo model.
        We assume that if the strand is None, we give it a generic "Position" only.
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

        iid = self._makePositionId(reference_id, position, position_types)
        n = self.gu.getNode(iid)
        pos = self.gu.getNode(self.properties['position'])
        ref = self.gu.getNode(self.properties['reference'])
        if position is not None:
            graph.add((n, pos, Literal(position, datatype=XSD['integer'])))
        graph.add((n, ref, self.gu.getNode(reference_id)))
        if position_types is not None:
            for t in position_types:
                graph.add((n, RDF['type'], self.gu.getNode(t)))
        s = None
        if strand is not None:
            s = strand
            if not re.match('faldo', strand):
                # not already mapped to faldo, so expect we need to map it
                s = self._getStrandType(strand)
        # else:
        #    s = self.types['both_strand']
        if s is None and (position_types is None or position_types == []):
            s = self.types['Position']

        if s is not None:
            graph.add((n, RDF['type'], self.gu.getNode(s)))

        return iid

    def addSubsequenceOfFeature(self, graph, parentid):
        """
        This will add reciprocal triples like:
        feature is_subsequence_of parent
        parent has_subsequence feature
        :param graph:
        :param parentid:
        :return:
        """
        self.gu.addTriple(graph, self.id, self.properties['is_subsequence_of'], parentid)
        self.gu.addTriple(graph, parentid, self.properties['has_subsequence'], self.id)

        return

    def addTaxonToFeature(self, graph, taxonid):
        """
        Given the taxon id, this will add the following triple:
        feature in_taxon taxonid
        :param graph:
        :param taxonid:
        :return:
        """
        self.taxon = taxonid
        self.gu.addTriple(graph, self.id, Assoc.properties['in_taxon'], self.taxon)

        return

    def loadAllProperties(self, graph):

        prop_dict = {
            Assoc(None).ANNOTPROP: self.annotation_properties,
            Assoc(None).OBJECTPROP: self.object_properties,
            Assoc(None).DATAPROP: self.data_properties
        }

        for p in prop_dict:
            self.gu.loadProperties(graph, prop_dict.get(p), p)

        return

    def addFeatureProperty(self, graph, property_type, property):
        self.gu.addTriple(graph, self.id, property_type, property)
        return

    def setNoBNodes(self, nobnodes):
        self.nobnodes = nobnodes
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
        logger.warn('no reference for this chrom.  you may have conflicting ids')
        taxon = ''

    # replace any chr-like prefixes with blank to standardize
    c = re.sub('ch(r?)[omse]*', '', str(chrom))

    # remove the build/taxon prefixes to look cleaner
    r = reference
    if re.match('.+:.+', reference):
        r = re.match('.*:(.*)', reference).group(1)
    id = ''.join((r, 'chr', c))
    if prefix is None:
        id = ''.join(('_', id))
    else:
        id = ':'.join((prefix, id))

    return id


def makeChromLabel(chrom, reference=None):
    label = ''
    c = re.sub('ch(r?)[omse\.]*', '', str(chrom))
    c = 'chr'+c
    if reference is None:
        label = c
    else:
        label = c+' ('+reference+')'

    return label
