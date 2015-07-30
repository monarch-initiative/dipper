import csv
import re
import os
from datetime import datetime
from stat import *
import logging

from dipper.utils.GraphUtils import GraphUtils
from dipper.sources.Source import Source
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Dataset import Dataset
from dipper.models.Reference import Reference
from dipper.models.Genotype import Genotype
from dipper import curie_map
from dipper import config


logger = logging.getLogger(__name__)


class MMRRC(Source):
    """

    """

    files = {
        'catalog': {'file': 'mmrrc_catalog_data.csv',
                   'url': 'https://www.mmrrc.org/about/mmrrc_catalog_data.csv'},
    }


    def __init__(self):
        Source.__init__(self, 'mmrrc')

        self.load_bindings()

        self.dataset = Dataset('mmrrc', 'Mutant Mouse Regional Resource Centers',
                               'https://www.mmrrc.org', None,
                               None)

        if 'test_ids' not in config.get_config() and 'gene' not in config.get_config()['test_ids']:
            logger.warn("not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['gene']

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)
        fname = '/'.join((self.rawdir,self.files['catalog']['file']))
        st = os.stat(fname)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        #TODO note can set the data version to what is in the header
        # first line like: This MMRRC catalog data file was generated on 2015-04-22

        self.dataset.setVersion(filedate)

        return


    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_phenotype_data( limit)

        logger.info("Finished parsing.")

        return


    def _process_phenotype_data(self, limit):
        '''

STRAIN/STOCK_ID,STRAIN/STOCK_DESIGNATION,STRAIN_TYPE,STATE,MGI_ALLELE_ACCESSION_ID,ALLELE_SYMBOL,ALLELE_NAME,MUTATION_TYPE,CHROMOSOME,MGI_GENE_ACCESSION_ID,GENE_SYMBOL,GENE_NAME,SDS_URL,ACCEPTED_DATE,MPT_IDS,PUBMED_IDS,RESEARCH_AREAS

        # NOTE: If a Strain carries more than one mutation, each Mutation description,
        i.e., the set: (Mutation Type - Chromosome - Gene Symbol - Gene Name - Allele Symbol - Allele Name)
        will require a separate line.


        :param raw:
        :param limit:
        :return:
        '''
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        gu = GraphUtils(curie_map.get())
        fname = '/'.join((self.rawdir,self.files['catalog']['file']))

        self.strain_hash = {}
        self.id_label_hash = {}
        mouse_taxon = 'NCBITaxon:10090'
        geno = Genotype(g)
        with open(fname, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for row in filereader:
                line_counter += 1
                # skip the first line, which is the version info
                # and skip the second line, which is the header
                if line_counter < 4:
                    continue

                (strain_id, strain_label, strain_type_symbol, strain_state,
                 mgi_allele_id, mgi_allele_symbol, mgi_allele_name, mutation_type, chrom,
                 mgi_gene_id, mgi_gene_symbol, mgi_gene_name, sds_url, accepted_date, mp_ids,
                pubmed_nums, research_areas) = row

                if self.testMode and (mgi_gene_id not in self.test_ids):
                    continue

                # get the variant or gene to save for later building of the genotype
                if strain_id not in self.strain_hash:
                    self.strain_hash[strain_id] = { 'variants': set(), 'genes': set() }

                if mgi_allele_id != '':
                    self.strain_hash[strain_id]['variants'].add(mgi_allele_id)
                    self.id_label_hash[mgi_allele_id] = mgi_allele_symbol
                    var_type = self._get_variant_type_from_abbrev(mutation_type)
                    # make a sequence alteration for this variant locus, and link the variation type to it
                    # sa_id = '_'+re.sub(':','',mgi_allele_id)+'SA'
                    # if self.nobnodes:
                    #     sa_id = ':'+sa_id
                    # gu.addIndividualToGraph(g, sa_id, None, var_type)
                    # geno.addSequenceAlterationToVariantLocus(sa_id, mgi_allele_id)
                if mgi_gene_id != '':
                    if re.match('Gene ID:', mgi_gene_id):
                        mgi_gene_id = re.sub('Gene ID:\s*', 'NCBIGene:', mgi_gene_id)
                    elif not re.match('MGI', mgi_gene_id):
                        logger.info("Gene id not recognized: %s", mgi_gene_id)
                    self.strain_hash[strain_id]['genes'].add(mgi_gene_id)
                    self.id_label_hash[mgi_gene_id] = mgi_gene_symbol

                # catch some errors - some things have gene labels, but no identifiers - report
                if mgi_gene_symbol.strip() != '' and mgi_gene_id == '':
                    logger.error("Found a gene label with no identifier for strain %s: %s", strain_id, mgi_gene_symbol)

                # split apart the mp ids
                # ataxia [MP:0001393] ,hypoactivity [MP:0001402] ...
                # mp_ids are now a comma delimited list with MP terms in brackets
                phenotype_ids = []
                if mp_ids != '':
                    for i in re.split(',', mp_ids):
                        i = i.strip()
                        mps = re.search('\[(.*)\]', i)
                        if mps is not None:
                            mp_id = mps.group(1).strip()
                            phenotype_ids.append(mp_id)

                #pubmed ids are space delimited
                pubmed_ids = []
                if pubmed_nums.strip() != '':
                    for i in re.split('\s+', pubmed_nums):
                        pmid = 'PMID:'+i.strip()
                        pubmed_ids.append(pmid)
                        r = Reference(pmid, Reference.ref_types['journal_article'])
                        r.addRefToGraph(g)

                # TODO the contents of the strain_num after the dash are the "holding center".  strip this out.
                # should map to https://www.mmrrc.org/catalog/sds.php?mmrrc_id=36643

                #https://www.mmrrc.org/catalog/sds.php?mmrrc_id=00001 is a good example of 4 genotype parts

                # TODO add strain_type
                gu.addClassToGraph(g, mouse_taxon, None)
                if research_areas.strip() == '':
                    research_areas = None
                else:
                    research_areas = 'Research Areas: '+research_areas
                strain_type = mouse_taxon
                if strain_state == 'ES':
                    stem_cell_class = 'CL:0000034'
                    strain_type = stem_cell_class
                gu.addIndividualToGraph(g, strain_id, strain_label, mouse_taxon, research_areas)  # an instance of mouse??

                for pid in phenotype_ids:
                    gu.addClassToGraph(g, pid, None)   #assume the phenotype label is in the ontology
                    assoc = G2PAssoc(self.name, strain_id, pid, gu.object_properties['has_phenotype'])
                    for p in pubmed_ids:
                        assoc.add_source(p)
                    assoc.add_association_to_graph(g)

                # TODO 1.  build a genotype from the alleles
                # TODO 2.  pull the background out of the strain


                if not self.testMode and (limit is not None and line_counter > limit):
                    break

            # now that we've collected all of the variant information, build it
            # we don't know their zygosities
            for s in self.strain_hash:
                h = self.strain_hash.get(s)
                vars = h['variants']
                genes = h['genes']
                vl_set = set()
                # make variant loci for each gene
                if len(vars) > 1:
                    logger.info("There's >1 variant for %s: %s", s, str(vars))
                elif len(vars) == 1:
                    vl_id = vars.pop()
                    vl_symbol = self.id_label_hash[vl_id]
                    if len(genes) > 0:
                        for gene in genes:
                            if len(vars) == 1:
                                self.id_label_hash[vl_id] = vl_symbol
                                geno.addAllele(vl_id, vl_symbol, geno.genoparts['variant_locus'])
                                geno.addAlleleOfGene(vl_id, gene, geno.object_properties['has_alternate_part'])
                                vl_set.add(vl_id)
                    else:
                        gu.addIndividualToGraph(g, vl_id, vl_symbol)  # is this a sequence alteration?
                        vl_set.add(vl_id)
                else:  #len(vars) == 0
                    for gene in genes:
                        vl_id = '_'+gene+'-VL'
                        if self.nobnodes:
                            vl_id = ':'+vl_id
                        vl_id = re.sub(':', '', vl_id)
                        vl_symbol = 'some variant of '+self.id_label_hash[gene]
                        self.id_label_hash[vl_id] = vl_symbol
                        geno.addAllele(vl_id, vl_symbol, geno.genoparts['variant_locus'])
                        geno.addAlleleOfGene(vl_id, gene, geno.object_properties['has_alternate_part'])
                        vl_set.add(vl_id)

                    # make the vslcs
                    vl_list = sorted(vl_set)
                    vslc_list = []
                    for vl in vl_list:
                        vslc_id = '_'+vl+'U'  # for unknown zygosity
                        if self.nobnodes:
                            vslc_id = ':'+vslc_id
                        vslc_id = re.sub(':', '', vslc_id)
                        vslc_label = self.id_label_hash[vl] + '/<?>'
                        self.id_label_hash[vslc_id] = vslc_label
                        vslc_list.append(vslc_id)
                        geno.addPartsToVSLC(vslc_id, vl, None, geno.zygosity['indeterminate'],
                                            geno.object_properties['has_alternate_part'], None)
                    if len(vslc_list) > 0:
                        gvc_id = '-'.join(vslc_list)
                        gvc_id = re.sub(':', '', gvc_id)
                        print('gvc=',gvc_id)
                        if self.nobnodes:
                            gvc_id = ':'+gvc_id
                        gvc_label = '; '.join(self.id_label_hash[v] for v in vslc_list)
                        print('gvc_label=',gvc_label)
                        gu.addIndividualToGraph(g, gvc_id, gvc_label, geno.genoparts['genomic_variation_complement'])
                        for vslc_id in vslc_list:
                            geno.addVSLCtoParent(vslc_id, gvc_id)

                        # todo if we get the appropriate background, then we can generate an appropriate genotype
                        # but for now, we will just link the strain to the gvc
                        gu.addTriple(g, s, geno.object_properties['has_genotype'], gvc_id)
                    else:
                        logger.info("Strain %s is not making a proper genotype.", s)


            gu.loadProperties(g, G2PAssoc.object_properties, G2PAssoc.OBJECTPROP)
            gu.loadProperties(g, G2PAssoc.datatype_properties, G2PAssoc.DATAPROP)
            gu.loadProperties(g, G2PAssoc.annotation_properties, G2PAssoc.ANNOTPROP)

        return

    def _get_variant_type_from_abbrev(self, abbrev):
        """
        All variants are generically typed as "sequence_alterations" unless otherwise stated.
        :param abbrev:
        :return:
        """
        variant_type = None

        var_dict = {
            'SM': 'SO:0001059',  # spontaneous mutation
            'TM': 'SO:0001059',  # targeted mutation
            'TG': 'SO:xxxxxxx',  # transgenic
            'GT': 'SO:0001059',  # gene trap
            'CI': 'SO:0001059', # chemically induced mutation
            'RAD': 'SO:0001059',  # radiation induced mutation
            'CH': 'SO:1000183',  # chromosomal aberration --> chromosomal structure variation
            'RB': 'SO:1000043',  # Robertsonian translocation
            'TL': 'SO:1000048',  #	reciprocal translocation
            'TP': 'SO:0000453', # transposition
            'INV': 'SO:1000036',  # inversion
            'INS': 'SO:0000667',  # insertion
            'DEL': 'SO:0000159',  # deletion
            'DP': 'SO:1000035',  #duplication
            'OTH': 'SO:0001059'  # other
        }
        if abbrev in var_dict:
            variant_type = var_dict[abbrev]
        else:
            logger.warn("Variant type not recognized: %s", abbrev)

        return variant_type

    def getTestSuite(self):
        import unittest
        from tests.test_mmrrc import MMRRCTestCase
        # TODO add G2PAssoc tests

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MMRRCTestCase)

        return test_suite
