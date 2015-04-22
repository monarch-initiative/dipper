__author__ = 'nlw'

import logging
import re

from rdflib import Literal
from rdflib.namespace import RDF, XSD

from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Assoc import Assoc

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
        'has_staining_intensity': 'GENO:0000626'
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

        # the following are sequence attributes:
        'band_intensity':  'GENO:0000618',
        'gneg': 'GENO:0000620',
        'gpos': 'GENO:0000619',
        'gpos100': 'GENO:0000622',
        'gpos75': 'GENO:0000623',
        'gpos50': 'GENO:0000624',
        'gpos25': 'GENO:0000625',
        'gvar': 'GENO:0000621'
    }

    def __init__(self, id, label, type, description=None):
        self.id = id
        self.label = label
        self.type = type
        self.description = description
        self.gu = GraphUtils(curie_map.get())
        self.start = None
        self.stop = None
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
        if position_types is None:
            loc['type'].append(self.types['Position'])
        else:
            loc['type'] += position_types

        return loc

    def _getStrandType(self, strand):
        """

        :param strand:
        :return:
        """
        # TODO make this a dictionary/enum:  PLUS, MINUS, BOTH
        strand_id = None
        if strand == '+':
            strand_id = self.types['plus_strand']
        elif strand == '-':
            strand_id = self.types['minus_strand']
        elif strand == '.' or strand is None:
            strand_id = self.types['both_strand']
        else:
            logger.warn("strand type could not be mapped: %s", strand)

        return strand_id

    def addFeatureToGraph(self, graph):
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
        self.gu.addIndividualToGraph(graph, self.id, self.label, self.type, self.description)

        # create a region that has the begin/end positions
        region_id = ':_'+self.id+'Region'  # FIXME make this anonymous
        self.gu.addTriple(graph, self.id, self.properties['location'], region_id)

        self.gu.addIndividualToGraph(graph, region_id, None, 'faldo:Region')
        # add the start/end positions to the region
        if self.start is not None:
            self.gu.addTriple(graph, region_id, self.properties['begin'],
                              self._makePositionId(self.start['reference'], self.start['coordinate']))
        if self.stop is not None:
            self.gu.addTriple(graph, region_id, self.properties['end'],
                              self._makePositionId(self.start['reference'], self.stop['coordinate']))

        # {coordinate : integer, reference : reference_id, types = []}

        if self.start is not None:
            self.addPositionToGraph(graph, self.start['reference'], self.start['coordinate'], self.start['type'])
        if self.stop is not None:
            self.addPositionToGraph(graph, self.stop['reference'], self.stop['coordinate'], self.stop['type'])

        return

    def _makePositionId(self, reference, coordinate, types=None):
        i = ':_'
        if reference is not None:
            i += reference + '-'
        i += str(coordinate)      # just in case it isn't a string already
        if types is not None:
            t = types.sort
            i += '-'.join(t)
        return i

    def addPositionToGraph(self,graph,reference_id,position,position_types=None,strand=None):
        """
        Add the positional information to the graph, following the faldo model.
        We assume that if the strand is None, it is meaning "Both".
        Triples:
        my_position a (any of: faldo:(((Both|Plus|Minus)Strand)|Exact)Position)
            faldo:position Integer(numeric position)
            faldo:reference reference_id

        :param graph:
        :param reference_id:
        :param position:
        :param position_types:
        :param strand:
        :return:
        """

        iid = self._makePositionId(reference_id, position)
        n = self.gu.getNode(iid)
        pos = self.gu.getNode(self.properties['position'])
        ref = self.gu.getNode(self.properties['reference'])
        graph.add((n, pos, Literal(position, datatype=XSD['integer'])))
        graph.add((n, ref, self.gu.getNode(reference_id)))
        if position_types is not None:
            for t in position_types:
                graph.add((n, RDF['type'], self.gu.getNode(t)))
        if strand is not None:
            s = strand
            if not re.match('faldo', strand):
                # not already mapped to faldo, so expect we need to map it
                s = self._getStrandType(strand)
        else:
            s = self.types['both_strand']
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
            Assoc().ANNOTPROP: self.annotation_properties,
            Assoc().OBJECTPROP: self.object_properties,
            Assoc().DATAPROP: self.data_properties
        }

        for p in prop_dict:
            self.gu.loadProperties(graph, prop_dict.get(p), p)

        return

    def addFeatureProperty(self, graph, property_type, property):
        self.gu.addTriple(graph, self.id, property_type, property)
        return


def makeChromID(chrom, reference=None):
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
    c = re.sub('ch(r?)[omse]*', '', chrom)

    # remove the build/taxon prefixes to look cleaner
    r = reference
    if re.match('.+:.+', reference):
        r = re.match('.*:(.*)', reference).group(1)
    id = ''.join((':', r, 'chr', c))
    return id


def makeChromLabel(chrom, reference=None):
    label = ''

    if reference is None:
        label = chrom
    else:
        label = chrom+' ('+reference+')'

    return label
