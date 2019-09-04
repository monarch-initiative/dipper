import re
import gzip
import logging

from dipper.sources.OMIMSource import OMIMSource
from dipper.models.Model import Model
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Genotype import Genotype
from dipper.models.GenomicFeature import Feature, makeChromID, makeChromLabel
from dipper.models.Reference import Reference
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class NCBIGene(OMIMSource):
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
    The gene ids in test_ids.yaml will be used to subset the data when testing.

    All entries in the gene_history file are added as deprecated classes,
    and linked to the current gene id, with "replaced_by" relationships.

    Since we do not know much about the specific link in the gene2pubmed;
    we simply create a "mentions" relationship.

    """

    files = {
        'gene_info': {
            'file': 'gene_info.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz',
            'columns': [
                'tax_id',
                'GeneID',
                'Symbol',
                'LocusTag',
                'Synonyms',
                'dbXrefs',
                'chromosome',
                'map_location',
                'description',
                'type_of_gene',
                'Symbol_from_nomenclature_authority',
                'Full_name_from_nomenclature_authority',
                'Nomenclature_status',
                'Other_designations',
                'Modification_date',
                'Feature_type',
            ]
        },
        'gene_history': {
            'file': 'gene_history.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz',
            'columns': [
                'tax_id',
                'GeneID',
                'Discontinued_GeneID',
                'Discontinued_Symbol',
                'Discontinue_Date',
            ]
        },
        'gene2pubmed': {
            'file': 'gene2pubmed.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz',
            'columns': [
                'tax_id',
                'GeneID',
                'PubMed_ID',
            ]
        },
        'gene_group': {
            'file': 'gene_group.gz',
            'url': 'http://ftp.ncbi.nih.gov/gene/DATA/gene_group.gz',
            'columns': [
                'tax_id',
                'GeneID',
                'relationship',
                'Other_tax_id',
                'Other_GeneID',
            ]
        }
    }

    resources = {
        'clique_leader': '../../resources/clique_leader.yaml'
    }

    def __init__(
            self,
            graph_type,
            are_bnodes_skolemized,
            tax_ids=None,
            gene_ids=None
    ):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'ncbigene',
            ingest_title='National Center for Biotechnology Information',
            ingest_url='http://ncbi.nih.nlm.gov/gene',
            # ingest_desc=None,
            license_url='https://creativecommons.org/publicdomain/mark/1.0/',
            data_rights='https://www.ncbi.nlm.nih.gov/home/about/policies/'
            # file_handle=None
        )

        self.tax_ids = tax_ids
        self.gene_ids = gene_ids
        self.id_filter = 'taxids'   # 'geneids

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]

        self.tax_ids = [str(x) for x in self.tax_ids]

        LOG.info("Filtering on the following taxa: %s", tax_ids)

        if 'gene' in self.all_test_ids:
            self.gene_ids = self.all_test_ids['gene']
        else:
            LOG.warning("not configured with gene test ids.")
            self.gene_ids = []

        self.class_or_indiv = {}

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self._get_gene_info(limit)
        self._get_gene_history(limit)
        self._get_gene2pubmed(limit)

        LOG.info("Done parsing files.")

    def _get_gene_info(self, limit):
        """
        Currently loops through the gene_info file and
        creates the genes as classes, typed with SO.  It will add their label,
        any alternate labels as synonyms, alternate ids as equivalent classes.
        HPRDs get added as protein products.
        The chromosome and chr band get added as blank node regions,
        and the gene is faldo:located
        on the chr band.
        :param limit:
        :return:

        """
        src_key = 'gene_info'
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        geno = Genotype(graph)
        model = Model(graph)

        # not unzipping the file
        LOG.info("Processing 'Gene Info' records")
        line_counter = 0
        gene_info = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("FILE: %s", gene_info)
        # Add taxa and genome classes for those in our filter

        band_regex = re.compile(r'[0-9A-Z]+[pq](\d+)?(\.\d+)?$')
        for tax_num in self.tax_ids:
            tax_id = ':'.join(('NCBITaxon', tax_num))
            # tax label can get added elsewhere
            geno.addGenome(tax_id, tax_num)
            # label added elsewhere
            model.addClassToGraph(tax_id, None,
                                  class_category=blv.terms.OrganismTaxon.value)

        col = self.files[src_key]['columns']
        with gzip.open(gene_info, 'rb') as tsv:
            row = tsv.readline().decode().strip().split('\t')
            row[0] = row[0][1:]  # strip comment
            if not self.check_fileheader(col, row):
                pass

            for line in tsv:
                line = line.strip()
                line_counter += 1
                if line[0] == '#':  # skip comments
                    continue
                row = line.decode().strip().split('\t')

                # ##set filter=None in init if you don't want to have a filter
                # if self.id_filter is not None:
                #     if ((self.id_filter == 'taxids' and \
                #          (tax_num not in self.tax_ids))
                #           or (self.id_filter == 'geneids' and \
                #               (int(gene_num) not in self.gene_ids))):
                #         continue
                # #### end filter

                tax_num = row[col.index('tax_id')]
                gene_num = row[col.index('GeneID')]
                symbol = row[col.index('Symbol')]
                # = row[col.index('LocusTag')]
                synonyms = row[col.index('Synonyms')].strip()
                dbxrefs = row[col.index('dbXrefs')].strip()
                chrom = row[col.index('chromosome')].strip()
                map_loc = row[col.index('map_location')].strip()
                desc = row[col.index('description')]
                gtype = row[col.index('type_of_gene')].strip()
                # = row[col.index('Symbol_from_nomenclature_authority')]
                name = row[col.index('Full_name_from_nomenclature_authority')]
                # = row[col.index('Nomenclature_status')]
                other_designations = row[col.index('Other_designations')].strip()
                # = row[col.index('Modification_date')}
                # = row[col.index('Feature_type')]

                if self.test_mode and int(gene_num) not in self.gene_ids:
                    continue
                if not self.test_mode and tax_num not in self.tax_ids:
                    continue
                tax_id = ':'.join(('NCBITaxon', tax_num))
                gene_id = ':'.join(('NCBIGene', gene_num))

                gene_type_id = self.resolve(gtype)

                if symbol == 'NEWENTRY':
                    label = None
                else:
                    label = symbol
                # sequence feature, not a gene
                if gene_type_id == self.globaltt['sequence_feature']:
                    self.class_or_indiv[gene_id] = 'I'
                else:
                    self.class_or_indiv[gene_id] = 'C'

                if not self.test_mode and limit is not None and line_counter > limit:
                    continue

                if self.class_or_indiv[gene_id] == 'C':
                    model.addClassToGraph(gene_id, label, gene_type_id, desc,
                                          class_category=blv.terms.Gene.value)
                    # NCBI will be the default leader (for non mods),
                    # so we will not add the leader designation here.
                else:
                    model.addIndividualToGraph(gene_id, label, gene_type_id, desc,
                                               ind_category=blv.terms.GenomicEntity.value)
                    # in this case, they aren't genes.
                    # so we want someone else to be the leader

                if name != '-':
                    model.addSynonym(gene_id, name, class_category=blv.terms.Gene.value)
                if synonyms != '-':
                    for syn in synonyms.split('|'):
                        model.addSynonym(
                            gene_id, syn.strip(), model.globaltt['has_related_synonym'],
                            class_category=blv.terms.Gene.value)
                if other_designations != '-':
                    for syn in other_designations.split('|'):
                        model.addSynonym(
                            gene_id, syn.strip(), model.globaltt['has_related_synonym'],
                            class_category=blv.terms.Gene.value)
                if dbxrefs != '-':
                    self._add_gene_equivalencies(dbxrefs, gene_id, tax_id)

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
                if chrom != '-' and chrom != '':
                    if re.search(r'\|', chrom) and chrom not in ['X|Y', 'X; Y']:
                        # means that there's uncertainty in the mapping.
                        # so skip it
                        # TODO we'll need to figure out how to deal with
                        # >1 loc mapping
                        LOG.info(
                            '%s is non-uniquely mapped to %s. Skipping for now.',
                            gene_id, chrom)
                        continue
                        # X|Y	Xp22.33;Yp11.3

                    # if(not re.match(
                    #        r'(\d+|(MT)|[XY]|(Un)$',str(chr).strip())):
                    #    print('odd chr=',str(chr))
                    if chrom == 'X; Y':
                        chrom = 'X|Y'  # rewrite the PAR regions for processing
                    # do this in a loop to allow PAR regions like X|Y
                    for chromosome in re.split(r'\|', chrom):
                        # assume that the chromosome label is added elsewhere
                        geno.addChromosomeClass(chromosome, tax_id, None)
                        mychrom = makeChromID(chromosome, tax_num, 'CHR')
                        # temporarily use taxnum for the disambiguating label
                        mychrom_syn = makeChromLabel(chromosome, tax_num)
                        model.addSynonym(mychrom, mychrom_syn)

                        band_match = re.match(band_regex, map_loc)
                        if band_match is not None and len(band_match.groups()) > 0:
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
                            bid = re.sub(r'^' + chromosome, '', map_loc)
                            # the generic location (no coordinates)
                            maploc_id = makeChromID(chromosome + bid, tax_num, 'CHR')
                            # print(map_loc,'-->',bid,'-->',maploc_id)
                            # Assume it's type will be added elsewhere
                            band = Feature(graph, maploc_id, None, None)
                            band.addFeatureToGraph()
                            # add the band as the containing feature
                            graph.addTriple(
                                gene_id,
                                self.globaltt['is subsequence of'],
                                maploc_id)
                        else:
                            # TODO handle these cases: examples are:
                            # 15q11-q22,Xp21.2-p11.23,15q22-qter,10q11.1-q24,
                            # 12p13.3-p13.2|12p13-p12,1p13.3|1p21.3-p13.1,
                            # 12cen-q21,22q13.3|22q13.3
                            LOG.debug(
                                'not regular band pattern for %s: %s', gene_id, map_loc)
                            # add the gene as a subsequence of the chromosome
                            graph.addTriple(
                                gene_id, self.globaltt['is subsequence of'], mychrom)

                geno.addTaxon(tax_id, gene_id)

    def _add_gene_equivalencies(self, xrefs, gene_id, taxon):
        """
        Add equivalentClass and sameAs relationships

        Uses external resource map located in
        /resources/clique_leader.yaml to determine
        if an NCBITaxon ID space is a clique leader
        """

        clique_map = self.open_and_parse_yaml(self.resources['clique_leader'])

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        filter_out = ['Vega', 'IMGT/GENE-DB', 'Araport']

        # deal with the dbxrefs
        # MIM:614444|HGNC:HGNC:16851|Ensembl:ENSG00000136828|HPRD:11479|Vega:OTTHUMG00000020696

        for dbxref in xrefs.strip().split('|'):
            prefix = ':'.join(dbxref.split(':')[:-1]).strip()
            if prefix in self.localtt:
                prefix = self.localtt[prefix]
            dbxref_curie = ':'.join((prefix, dbxref.split(':')[-1]))

            if dbxref_curie is not None and prefix != '':
                if prefix == 'HPRD':  # proteins are not == genes.
                    model.addTriple(
                        gene_id, self.globaltt['has gene product'], dbxref_curie,
                        subject_category=blv.terms.Gene.value,
                        object_category=blv.terms.GeneProduct.value)
                    continue
                    # skip some of these for now based on curie prefix
                if prefix in filter_out:
                    continue

                if prefix == 'ENSEMBL':
                    model.addXref(gene_id, dbxref_curie, class_category=blv.terms.Gene.value)
                if prefix == 'OMIM':
                    if dbxref_curie in self.omim_replaced:
                        repl = self.omim_replaced[dbxref_curie]
                        for omim in repl:
                            if omim in self.omim_type and \
                                    self.omim_type[omim] == self.globaltt['gene']:
                                dbxref_curie = omim
                    if dbxref_curie in self.omim_type and \
                            self.omim_type[dbxref_curie] != self.globaltt['gene']:
                        continue
                try:
                    if self.class_or_indiv.get(gene_id) == 'C':
                        model.addEquivalentClass(gene_id, dbxref_curie,
                                                 subject_category=blv.terms.Gene.value,
                                                 object_category=blv.terms.Gene.value)
                        if taxon in clique_map:
                            if clique_map[taxon] == prefix:
                                model.makeLeader(dbxref_curie,
                                                 node_category=blv.terms.Gene.value)
                            elif clique_map[taxon] == gene_id.split(':')[0]:
                                model.makeLeader(gene_id,
                                                 node_category=blv.terms.Gene.value)
                    else:
                        model.addSameIndividual(gene_id, dbxref_curie,
                                                subject_category=blv.terms.Gene.value,
                                                object_category=blv.terms.Gene.value)
                except AssertionError as err:
                    LOG.warning("Error parsing %s: %s", gene_id, err)

    def _get_gene_history(self, limit):
        """
        Loops through the gene_history file and adds the old gene ids
        as deprecated classes, where the new gene id is the replacement for it.
        The old gene symbol is added as a synonym to the gene.
        :param limit:
        :return:

        """
        src_key = 'gene_history'
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        LOG.info("Processing Gene records")
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("FILE: %s", myfile)
        col = self.files[src_key]['columns']
        with gzip.open(myfile, 'rb') as tsv:
            row = tsv.readline().decode().strip().split('\t')
            row[0] = row[0][1:]  # strip comment
            if not self.check_fileheader(col, row):
                pass

            for line in tsv:
                # skip comments
                row = line.decode().strip().split('\t')
                if row[0][0] == '#':
                    continue

                tax_num = row[col.index('tax_id')].strip()
                gene_num = row[col.index('GeneID')].strip()
                discontinued_num = row[col.index('Discontinued_GeneID')].strip()
                discontinued_symbol = row[col.index('Discontinued_Symbol')].strip()
                # discontinued_date = row[col.index('Discontinue_Date')]

                # set filter=None in init if you don't want to have a filter
                # if self.id_filter is not None:
                #     if ((self.id_filter == 'taxids' and \
                #          (int(tax_num) not in self.tax_ids))
                #             or (self.id_filter == 'geneids' and \
                #                 (int(gene_num) not in self.gene_ids))):
                #         continue
                #  end filter

                if gene_num == '-' or discontinued_num == '-':
                    continue

                if self.test_mode and gene_num not in self.gene_ids:
                    continue

                if not self.test_mode and tax_num not in self.tax_ids:
                    continue

                line_counter += 1
                gene_id = ':'.join(('NCBIGene', gene_num))
                discontinued_gene_id = ':'.join(('NCBIGene', discontinued_num))
                # add the two genes
                if self.class_or_indiv.get(gene_id) == 'C':
                    model.addClassToGraph(gene_id, None,
                                          class_category=blv.terms.Gene.value)
                    model.addClassToGraph(discontinued_gene_id, discontinued_symbol,
                                          class_category=blv.terms.Gene.value)

                    # add the new gene id to replace the old gene id
                    model.addDeprecatedClass(discontinued_gene_id, [gene_id],
                                             old_id_category=blv.terms.Gene.value,
                                             new_ids_category=blv.terms.Gene.value)
                else:
                    model.addIndividualToGraph(gene_id, None,
                                               ind_category=blv.terms.Gene.value)
                    model.addIndividualToGraph(
                        discontinued_gene_id, discontinued_symbol,
                        ind_category=blv.terms.Gene.value)
                    model.addDeprecatedIndividual(discontinued_gene_id, [gene_id],
                                                  old_id_category=blv.terms.Gene.value,
                                                  new_ids_category=blv.terms.Gene.value)

                # also add the old symbol as a synonym of the new gene
                model.addSynonym(gene_id, discontinued_symbol,
                                 class_category=blv.terms.Gene.value)

                if not self.test_mode and (limit is not None and line_counter > limit):
                    break

    def _get_gene2pubmed(self, limit):
        """
        Loops through the gene2pubmed file and adds a simple triple to say
        that a given publication is_about a gene.
        Publications are added as NamedIndividuals.

        These are filtered on the taxon.

        :param limit:
        :return:

        """
        src_key = 'gene2pubmed'
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        LOG.info("Processing Gene records")
        line_counter = 0
        myfile = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("FILE: %s", myfile)
        assoc_counter = 0
        col = self.files[src_key]['columns']
        with gzip.open(myfile, 'rb') as tsv:
            row = tsv.readline().decode().strip().split('\t')
            row[0] = row[0][1:]  # strip comment
            if not self.check_fileheader(col, row):
                pass

            for line in tsv:
                line_counter += 1
                # skip comments
                row = line.decode().strip().split('\t')
                if row[0][0] == '#':
                    continue
                tax_num = row[col.index('tax_id')].strip()
                gene_num = row[col.index('GeneID')].strip()
                pubmed_num = row[col.index('PubMed_ID')].strip()

                # ## set id_filter=None in init if you don't want to have a filter
                # if self.id_filter is not None:
                #     if ((self.id_filter == 'taxids' and \
                #          (int(tax_num) not in self.tax_ids))
                #        or (self.id_filter == 'geneids' and \
                #            (int(gene_num) not in self.gene_ids))):
                #         continue
                # #### end filter

                if self.test_mode and int(gene_num) not in self.gene_ids:
                    continue

                if not self.test_mode and tax_num not in self.tax_ids:
                    continue

                if gene_num == '-' or pubmed_num == '-':
                    continue

                gene_id = ':'.join(('NCBIGene', gene_num))
                pubmed_id = ':'.join(('PMID', pubmed_num))

                if self.class_or_indiv.get(gene_id) == 'C':
                    model.addClassToGraph(gene_id, None,
                                          class_category=blv.terms.Gene.value)
                else:
                    model.addIndividualToGraph(gene_id, None,
                                               ind_category=blv.terms.Gene.value)
                # add the publication as a NamedIndividual
                # add type publication
                model.addIndividualToGraph(pubmed_id, None, None,
                                           ind_category=blv.terms.Publication.value)
                reference = Reference(
                    graph, pubmed_id, self.globaltt['journal article'])
                reference.addRefToGraph()
                graph.addTriple(
                    pubmed_id, self.globaltt['is_about'], gene_id,
                    subject_category=blv.terms.Publication.value,
                    object_category=blv.terms.Gene.value)
                assoc_counter += 1
                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        LOG.info("Processed %d pub-gene associations", assoc_counter)

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
        src_key = 'gene_group'
        LOG.info("getting gene groups")
        src_file = '/'.join((self.rawdir, self.files[src_key]['file']))
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
        col = self.files[src_key]['columns']

        with gzip.open(src_file, 'rb') as tsv:
            row = tsv.readline().decode().strip().split('\t')
            row[0] = row[0][1:]  # strip octothorp
            if not self.check_fileheader(col, row):
                pass
            for row in tsv:
                row = row.decode().strip().split('\t')
                tax_a = row[col.index('tax_id')]
                gene_a = row[col.index('GeneID')]
                rel = row[col.index('relationship')]
                tax_b = row[col.index('Other_tax_id')]
                gene_b = row[col.index('Other_GeneID')]

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
        LOG.debug("Finished hashing gene groups")
        LOG.debug("Making orthology associations")
        for gid in gene_ids:
            gene_num = re.sub(r'NCBIGene:', '', gid)
            group_nums = gene_to_group.get(gene_num)
            if group_nums is not None:
                for group_num in group_nums:
                    orthologs = group_to_orthology.get(group_num)
                    if orthologs is not None:
                        for orth in orthologs:
                            oid = 'NCBIGene:' + str(orth)
                            model.addClassToGraph(oid, None, self.globaltt['gene'],
                                                  class_category=blv.terms.Gene.value)
                            otaxid = 'NCBITaxon:' + str(gene_to_taxon[orth])
                            geno.addTaxon(otaxid, oid, genopart_category=blv.terms.Gene.value)
                            assoc = OrthologyAssoc(graph, self.name, gid, oid)
                            assoc.add_source('PMID:24063302')
                            assoc.add_association_to_graph()
                            # todo get gene label for orthologs -
                            # this could get expensive
                            found_counter += 1

            # finish loop through annotated genes
        LOG.info(
            "Made %d orthology relationships for %d genes",
            found_counter, len(gene_ids))
