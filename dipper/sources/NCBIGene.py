import re
import gzip
import logging
import csv
import io
import requests

from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
from dipper.models.assoc.Association import Assoc
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Genotype import Genotype
from dipper import config
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.models.Reference import Reference
from dipper.utils.DipperUtil import DipperUtil


logger = logging.getLogger(__name__)


class NCBIGene(Source):
    """
    This is the processing module for the
    National Center for Biotechnology Information.  It includes parsers for
    the gene_info (gene names, symbols, ids, equivalent ids), gene history
    (alt ids), and gene2pubmed publication references about a gene.

    This creates Genes as classes, when they are properly typed as such.
    For those entries where it is an 'unknown significance', it is added simply
    as an instance of a sequence feature.  It will add equivalentClasses for
    a subset of external identifiers, including:
    ENSEMBL, HGMD, MGI, ZFIN, and gene product links for HPRD.
    They are additionally located to their Chromosomal band
    (until we process actual genomic coords in a separate file).

    We process the genes from the filtered taxa, starting with those configured
    by default (human, mouse, fish).
    This can be overridden in the calling script to include additional taxa,
    if desired.
    The gene ids in the conf.json will be used to subset the data when testing.

    All entries in the gene_history file are added as deprecated classes,
    and linked to the current gene id, with "replaced_by" relationships.

    Since we do not know much about the specific link in the gene2pubmed;
    we simply create a "mentions" relationship.

    """

    SCIGRAPH_BASE = 'https://scigraph-ontology-dev.monarchinitiative.org/scigraph/graph/'

    files = {
        'gene_info': {
            'file': 'gene_info.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'
        },
        'gene_history': {
            'file': 'gene_history.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz'
        },
        'gene2pubmed': {
            'file': 'gene2pubmed.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz'
        },
        'gene_group': {
            'file': 'gene_group.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_group.gz'}
    }

    resources = {
        'clique_leader': '../../resources/clique_leader.yaml'
    }

    def __init__(self, graph_type, are_bnodes_skolemized,
                 tax_ids=None, gene_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'ncbigene')

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids
        self.filter = 'taxids'

        self.dataset = Dataset(
            'ncbigene', 'National Center for Biotechnology Information',
            'http://ncbi.nih.nlm.gov/gene', None,
            'http://www.ncbi.nlm.nih.gov/About/disclaimer.html',
            'https://creativecommons.org/publicdomain/mark/1.0/')

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]
            logger.info("No taxa set.  Defaulting to %s", str(tax_ids))
        else:
            logger.info("Filtering on the following taxa: %s", str(tax_ids))

        self.gene_ids = []
        if 'test_ids' not in \
                config.get_config() or \
                'gene' not in config.get_config()['test_ids']:
            logger.warning("not configured with gene test ids.")
        else:
            self.gene_ids = config.get_config()['test_ids']['gene']

        self.properties = Feature.properties

        self.class_or_indiv = {}

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

        self._get_gene_info(limit)
        self._get_gene_history(limit)
        self._get_gene2pubmed(limit)

        logger.info("Done parsing files.")

        return

    def _get_gene_info(self, limit):
        """
        Currently loops through the gene_info file and
        creates the genes as classes, typed with SO.  It will add their label,
        any alternate labels as synonyms, alternate ids as equivlaent classes.
        HPRDs get added as protein products.
        The chromosome and chr band get added as blank node regions,
        and the gene is faldo:located
        on the chr band.
        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)
        model = Model(g)

        # not unzipping the file
        logger.info("Processing 'Gene Info' records")
        line_counter = 0
        gene_info = '/'.join((self.rawdir, self.files['gene_info']['file']))
        logger.info("FILE: %s", gene_info)
        # Add taxa and genome classes for those in our filter
        for tax_num in self.tax_ids:
            tax_id = ':'.join(('NCBITaxon', str(tax_num)))
            # tax label can get added elsewhere
            geno.addGenome(tax_id, str(tax_num))
            # label added elsewhere
            model.addClassToGraph(tax_id, None)
        with gzip.open(gene_info, 'rb') as f:
            row = f.readline().decode().strip().split('\t')
            logger.info("Header has %i columns", len(row))
            for line in f:
                # skip comments
                line = line.decode().strip()
                if re.match(r'^#', line):
                    continue
                (tax_num, gene_num, symbol, locustag, synonyms, xrefs, chrom,
                 map_loc, desc, gtype, authority_symbol, name,
                 nomenclature_status, other_designations,
                 modification_date, feature_type) = line.split('\t')

                # ##set filter=None in init if you don't want to have a filter
                # if self.filter is not None:
                #     if ((self.filter == 'taxids' and \
                #          (int(tax_num) not in self.tax_ids))
                #           or (self.filter == 'geneids' and \
                #               (int(gene_num) not in self.gene_ids))):
                #         continue
                # #### end filter

                if self.testMode and int(gene_num) not in self.gene_ids:
                    continue

                if not self.testMode and int(tax_num) not in self.tax_ids:
                    continue

                line_counter += 1

                gene_id = ':'.join(('NCBIGene', gene_num))
                tax_id = ':'.join(('NCBITaxon', tax_num))
                gene_type_id = self.map_type_of_gene(gtype.strip())

                if symbol == 'NEWENTRY':
                    label = None
                else:
                    label = symbol
                # sequence feature, not a gene
                if gene_type_id == 'SO:0000110':
                    self.class_or_indiv[gene_id] = 'I'
                else:
                    self.class_or_indiv[gene_id] = 'C'

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    continue

                if self.class_or_indiv[gene_id] == 'C':
                    model.addClassToGraph(gene_id, label, gene_type_id, desc)
                    # NCBI will be the default leader,
                    # so we will not add the leader designation here.
                else:
                    model.addIndividualToGraph(
                        gene_id, label, gene_type_id, desc)
                    # in this case, they aren't genes.
                    # so we want someone else to be the leader.

                if name != '-':
                    model.addSynonym(gene_id, name)
                if synonyms.strip() != '-':
                    for s in synonyms.split('|'):
                        model.addSynonym(
                            gene_id, s.strip(),
                            Assoc.annotation_properties['hasRelatedSynonym'])
                if other_designations.strip() != '-':
                    for s in other_designations.split('|'):
                        model.addSynonym(
                            gene_id, s.strip(),
                            Assoc.annotation_properties['hasRelatedSynonym'])
                if xrefs.strip() != '-':
                    self._add_gene_equivalencies(xrefs, gene_id, tax_num)

                # edge cases of id | symbol | chr | map_loc:
                # 263     AMD1P2    X|Y  with   Xq28 and Yq12
                # 438     ASMT      X|Y  with   Xp22.3 or Yp11.3    # in PAR
                # no idea why there's two bands listed - possibly 2 assemblies
                # 419     ART3      4    with   4q21.1|4p15.1-p14
                # 28227   PPP2R3B   X|Y  Xp22.33; Yp11.3            # in PAR
                # this is of "unknown" type == susceptibility
                # 619538  OMS     10|19|3 10q26.3;19q13.42-q13.43;3p25.3
                # unlocated scaffold
                # 101928066       LOC101928066    1|Un    -\
                # mouse --> 2C3
                # 11435   Chrna1  2       2 C3|2 43.76 cM
                # mouse --> 11B1.1
                # 11548   Adra1b  11      11 B1.1|11 25.81 cM
                # 11717   Ampd3   7       7 57.85 cM|7 E2-E3        # mouse
                # 14421   B4galnt1        10      10 D3|10 74.5 cM  # mouse
                # 323212  wu:fb92e12      19|20   -                 # fish
                # 323368  ints10  6|18    -                         # fish
                # 323666  wu:fc06e02      11|23   -                 # fish

                # feel that the chr placement can't be trusted in this table
                # when there is > 1 listed
                # with the exception of human X|Y,
                # we will only take those that align to one chr

                # FIXME remove the chr mapping below
                # when we pull in the genomic coords
                if str(chrom) != '-' and str(chrom) != '':
                    if re.search(r'\|', str(chrom)) and \
                            str(chrom) not in ['X|Y', 'X; Y']:
                        # means that there's uncertainty in the mapping.
                        # so skip it
                        # TODO we'll need to figure out how to deal with
                        # >1 loc mapping
                        logger.info(
                            '%s is non-uniquely mapped to %s.' +
                            ' Skipping for now.',
                            gene_id, str(chr))
                        continue
                        # X|Y	Xp22.33;Yp11.3

                    # if(not re.match(
                    #        r'(\d+|(MT)|[XY]|(Un)$',str(chr).strip())):
                    #    print('odd chr=',str(chr))
                    if str(chrom) == 'X; Y':
                        chrom = 'X|Y'  # rewrite the PAR regions for processing
                    # do this in a loop to allow PAR regions like X|Y
                    for c in re.split(r'\|', str(chrom)):
                        # assume that the chromosome label is added elsewhere
                        geno.addChromosomeClass(c, tax_id, None)
                        mychrom = makeChromID(c, tax_num, 'CHR')
                        # temporarily use taxnum for the disambiguating label
                        mychrom_syn = makeChromLabel(c, tax_num)
                        model.addSynonym(mychrom, mychrom_syn)
                        band_match = re.match(
                            r'[0-9A-Z]+[pq](\d+)?(\.\d+)?$', map_loc)
                        if band_match is not None and \
                                len(band_match.groups()) > 0:
                            # if tax_num != '9606':
                            #     continue
                            # this matches the regular kind of chrs,
                            # so make that kind of band
                            # not sure why this matches?
                            #   chrX|Y or 10090chr12|Un"
                            # TODO we probably need a different regex
                            # per organism
                            # the maploc_id already has the numeric chromosome
                            # in it, strip it first
                            bid = re.sub(r'^'+c, '', map_loc)
                            # the generic location (no coordinates)
                            maploc_id = makeChromID(c+bid, tax_num, 'CHR')
                            # print(map_loc,'-->',bid,'-->',maploc_id)
                            # Assume it's type will be added elsewhere
                            band = Feature(g, maploc_id, None, None)
                            band.addFeatureToGraph()
                            # add the band as the containing feature
                            g.addTriple(
                                gene_id,
                                Feature.object_properties['is_subsequence_of'],
                                maploc_id)
                        else:
                            # TODO handle these cases: examples are:
                            # 15q11-q22,Xp21.2-p11.23,15q22-qter,10q11.1-q24,
                            # 12p13.3-p13.2|12p13-p12,1p13.3|1p21.3-p13.1,
                            # 12cen-q21,22q13.3|22q13.3
                            logger.debug(
                                'not regular band pattern for %s: %s',
                                gene_id, map_loc)
                            # add the gene as a subsequence of the chromosome
                            g.addTriple(
                                gene_id,
                                Feature.object_properties['is_subsequence_of'],
                                mychrom)

                geno.addTaxon(tax_id, gene_id)

        return

    def _add_gene_equivalencies(self, xrefs, gene_id, taxon):
        """
        Add equivalentClass and sameAs relationships

        Uses external resource map located in
        /resources/clique_leader.yaml to determine
        if an ID space is a clique leader
        """

        clique_map = self.open_and_parse_yaml(self.resources['clique_leader'])

        if self.testMode:
            graph = self.testgraph
        else:
            graph = self.graph

        filter_out = ['Vega', 'IMGT/GENE-DB', 'Araport']
        taxon_spec_filters = {
            '10090': ['ENSEMBL']
        }
        if taxon in taxon_spec_filters:
            filter_out += taxon_spec_filters[taxon]

        model = Model(graph)
        # deal with the xrefs
        # MIM:614444|HGNC:HGNC:16851|Ensembl:ENSG00000136828|HPRD:11479|Vega:OTTHUMG00000020696
        for ref in xrefs.strip().split('|'):
            xref_curie = self._cleanup_id(ref)
            if xref_curie is not None and xref_curie.strip() != '':
                if re.match(r'HPRD', xref_curie):
                    # proteins are not == genes.
                    model.addTriple(
                        gene_id,
                        self.properties['has_gene_product'], xref_curie)
                    continue
                    # skip some of these for now
                if xref_curie.split(':')[0] in filter_out:
                    continue
                if re.match(r'^OMIM', xref_curie):
                    if DipperUtil.is_omim_disease(xref_curie):
                        continue
                try:
                    if self.class_or_indiv.get(gene_id) == 'C':
                        model.addEquivalentClass(
                            gene_id, xref_curie)
                        if int(taxon) in clique_map:
                            if clique_map[int(taxon)] == xref_curie.split(':')[0]:
                                model.makeLeader(xref_curie)
                            elif clique_map[int(taxon)] == gene_id.split(':')[0]:
                                model.makeLeader(gene_id)
                    else:
                        model.addSameIndividual(gene_id, xref_curie)
                except AssertionError as e:
                    logger.warn("Error parsing {0}: {1}".format(gene_id, e))
        return

    def _get_gene_history(self, limit):
        """
        Loops through the gene_history file and adds the old gene ids
        as deprecated classes, where the new gene id is the replacement for it.
        The old gene symbol is added as a synonym to the gene.
        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("Processing Gene records")
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files['gene_history']['file']))
        logger.info("FILE: %s", myfile)
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                # skip comments
                line = line.decode().strip()
                if re.match(r'^#', line):
                    continue
                (tax_num, gene_num, discontinued_num, discontinued_symbol,
                 discontinued_date) = line.split('\t')

                # set filter=None in init if you don't want to have a filter
                # if self.filter is not None:
                #     if ((self.filter == 'taxids' and \
                #          (int(tax_num) not in self.tax_ids))
                #             or (self.filter == 'geneids' and \
                #                 (int(gene_num) not in self.gene_ids))):
                #         continue
                #  end filter

                if gene_num == '-' or discontinued_num == '-':
                    continue

                if self.testMode and int(gene_num) not in self.gene_ids:
                    continue

                if not self.testMode and int(tax_num) not in self.tax_ids:
                    continue

                line_counter += 1
                gene_id = ':'.join(('NCBIGene', gene_num))
                discontinued_gene_id = ':'.join(('NCBIGene', discontinued_num))

                # add the two genes
                if self.class_or_indiv.get(gene_id) == 'C':
                    model.addClassToGraph(gene_id, None)
                    model.addClassToGraph(
                        discontinued_gene_id, discontinued_symbol)

                    # add the new gene id to replace the old gene id
                    model.addDeprecatedClass(discontinued_gene_id, [gene_id])
                else:
                    model.addIndividualToGraph(gene_id, None)
                    model.addIndividualToGraph(
                        discontinued_gene_id, discontinued_symbol)
                    model.addDeprecatedIndividual(
                        discontinued_gene_id, [gene_id])

                # also add the old symbol as a synonym of the new gene
                model.addSynonym(gene_id, discontinued_symbol)

                if (not self.testMode) and\
                        (limit is not None and line_counter > limit):
                    break

        return

    def _get_gene2pubmed(self, limit):
        """
        Loops through the gene2pubmed file and adds a simple triple to say
        that a given publication is_about a gene.
        Publications are added as NamedIndividuals.

        These are filtered on the taxon.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("Processing Gene records")
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files['gene2pubmed']['file']))
        logger.info("FILE: %s", myfile)
        assoc_counter = 0
        with gzip.open(myfile, 'rb') as f:
            for line in f:
                # skip comments
                line = line.decode().strip()
                if re.match(r'^#', line):
                    continue
                (tax_num, gene_num, pubmed_num) = line.split('\t')

                # ## set filter=None in init if you don't want to have a filter
                # if self.filter is not None:
                #     if ((self.filter == 'taxids' and \
                #          (int(tax_num) not in self.tax_ids))
                #        or (self.filter == 'geneids' and \
                #            (int(gene_num) not in self.gene_ids))):
                #         continue
                # #### end filter

                if self.testMode and int(gene_num) not in self.gene_ids:
                    continue

                if not self.testMode and int(tax_num) not in self.tax_ids:
                    continue

                if gene_num == '-' or pubmed_num == '-':
                    continue

                line_counter += 1
                gene_id = ':'.join(('NCBIGene', gene_num))
                pubmed_id = ':'.join(('PMID', pubmed_num))

                if self.class_or_indiv.get(gene_id) == 'C':
                    model.addClassToGraph(gene_id, None)
                else:
                    model.addIndividualToGraph(gene_id, None)
                # add the publication as a NamedIndividual
                # add type publication
                model.addIndividualToGraph(pubmed_id, None, None)
                reference = Reference(
                    g, pubmed_id, Reference.ref_types['journal_article'])
                reference.addRefToGraph()
                g.addTriple(
                    pubmed_id, model.object_properties['is_about'], gene_id)
                assoc_counter += 1
                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        logger.info(
            "Processed %d pub-gene associations", assoc_counter)

        return

    @staticmethod
    def map_type_of_gene(sotype):
        so_id = 'SO:0000110'
        type_to_so_map = {
            'ncRNA': 'SO:0001263',
            'other': 'SO:0000110',
            'protein-coding': 'SO:0001217',
            'pseudo': 'SO:0000336',
            'rRNA': 'SO:0001637',
            'snRNA': 'SO:0001268',
            'snoRNA': 'SO:0001267',
            'tRNA': 'SO:0001272',
            'unknown': 'SO:0000110',
            'scRNA': 'SO:0001266',
            # mature transcript - there is no good mapping
            'miscRNA': 'SO:0000233',
            'chromosome': 'SO:0000340',
            'chromosome_arm': 'SO:0000105',
            'chromosome_band': 'SO:0000341',
            'chromosome_part': 'SO:0000830'
        }

        if sotype in type_to_so_map:
            so_id = type_to_so_map.get(sotype)
        else:
            logger.warning(
                "unmapped code %s. Defaulting to 'SO:0000110', " +
                "sequence_feature.", sotype)

        return so_id

    @staticmethod
    def _cleanup_id(i):
        """
        Clean up messy id prefixes
        :param i:
        :return:
        """
        cleanid = i
        # MIM:123456 --> #OMIM:123456
        cleanid = re.sub(r'^MIM', 'OMIM', cleanid)

        # HGNC:HGNC --> HGNC
        cleanid = re.sub(r'^HGNC:HGNC', 'HGNC', cleanid)

        # Ensembl --> ENSEMBL
        cleanid = re.sub(r'^Ensembl', 'ENSEMBL', cleanid)

        # MGI:MGI --> MGI
        cleanid = re.sub(r'^MGI:MGI', 'MGI', cleanid)

        cleanid = re.sub(r'FLYBASE', 'FlyBase', cleanid)

        return cleanid

    def getTestSuite(self):
        import unittest
        from tests.test_ncbi import NCBITestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(NCBITestCase)

        return test_suite

    def add_orthologs_by_gene_group(self, graph, gene_ids):
        """
        This will get orthologies between human and other vertebrate genomes
        based on the gene_group annotation pipeline from NCBI.
        More information 9can be learned here:
        http://www.ncbi.nlm.nih.gov/news/03-13-2014-gene-provides-orthologs-regions/
        The method for associations is described in
        [PMCID:3882889](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882889/)
        == [PMID:24063302](http://www.ncbi.nlm.nih.gov/pubmed/24063302/).
        Because these are only between human and vertebrate genomes,
        they will certainly miss out on very distant orthologies,
        and should not be considered complete.

        We do not run this within the NCBI parser itself;
        rather it is a convenience function for others parsers to call.

        :param graph:
        :param gene_ids:  Gene ids to fetch the orthology
        :return:

        """

        logger.info("getting gene groups")
        line_counter = 0
        f = '/'.join((self.rawdir, self.files['gene_group']['file']))
        found_counter = 0
        # because many of the orthologous groups are grouped by human gene,
        # we need to do this by generating two-way hash

        # group_id => orthologs
        # ortholog id => group
        # this will be the fastest approach, though not memory-efficient.
        geno = Genotype(graph)
        model = Model(graph)
        group_to_orthology = {}
        gene_to_group = {}
        gene_to_taxon = {}

        with gzip.open(f, 'rb') as csvfile:
            filereader = csv.reader(
                io.TextIOWrapper(csvfile, newline=""),
                delimiter='\t',
                quotechar='\"')

            for row in filereader:
                # skip comment lines
                if re.match(r'\#', ''.join(row)):
                    continue
                line_counter += 1
                (tax_a, gene_a, rel, tax_b, gene_b) = row

                if rel != 'Ortholog':
                    continue

                if gene_a not in group_to_orthology:
                    group_to_orthology[gene_a] = set()
                group_to_orthology[gene_a].add(gene_b)

                if gene_b not in gene_to_group:
                    gene_to_group[gene_b] = set()
                gene_to_group[gene_b].add(gene_a)

                gene_to_taxon[gene_a] = tax_a
                gene_to_taxon[gene_b] = tax_b

                # also add the group lead as a member of the group
                group_to_orthology[gene_a].add(gene_a)

            # end loop through gene_group file
        logger.debug("Finished hashing gene groups")
        logger.debug("Making orthology associations")
        for gid in gene_ids:
            gene_num = re.sub(r'NCBIGene:', '', gid)
            group_nums = gene_to_group.get(gene_num)
            if group_nums is not None:
                for group_num in group_nums:
                    orthologs = group_to_orthology.get(group_num)
                    if orthologs is not None:
                        for o in orthologs:
                            oid = 'NCBIGene:'+str(o)
                            model.addClassToGraph(
                                oid, None, Genotype.genoparts['gene'])
                            otaxid = 'NCBITaxon:'+str(gene_to_taxon[o])
                            geno.addTaxon(otaxid, oid)
                            assoc = OrthologyAssoc(graph, self.name, gid, oid)
                            assoc.add_source('PMID:24063302')
                            assoc.add_association_to_graph()
                            # todo get gene label for orthologs -
                            # this could get expensive
                            found_counter += 1

            # finish loop through annotated genes
        logger.info(
            "Made %d orthology relationships for %d genes",
            found_counter, len(gene_ids))
        return
