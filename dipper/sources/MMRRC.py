import csv
import re
import os
from datetime import datetime
from stat import ST_CTIME
import logging

from dipper.sources.Source import Source
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Dataset import Dataset
from dipper.models.Reference import Reference
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model

logger = logging.getLogger(__name__)


class MMRRC(Source):
    """
    Here we process the Mutant Mouse Resource and Research Center
    (https://www.mmrrc.org) strain data,
    which includes:
    *  strains, their mutant alleles
    *  phenotypes of the alleles
    *  descriptions of the research uses of the strains

    Note that some gene identifiers are not included
    (for many of the transgenics with human genes) in the raw data.
    We do our best to process the links between the variant and
    the affected gene, but sometimes the mapping is not clear,
    and we do not include it.
    Many of these details will be solved by merging this source with
    the MGI data source, who has the variant-to-gene designations.

    Also note that even though the strain pages at the MMRRC site do list
    phenotypic differences in the context of the strain backgrounds,
    they do not provide that data to us,
    and thus we cannot supply that disambiguation here.
    """

    files = {
        'catalog': {
            'file': 'mmrrc_catalog_data.csv',
            'url': 'https://www.mmrrc.org/about/mmrrc_catalog_data.csv'},
    }

    test_ids = [
        'MMRRC:037507-MU', 'MMRRC:041175-UCD', 'MMRRC:036933-UNC',
        'MMRRC:037884-UCD', 'MMRRC:000255-MU', 'MMRRC:037372-UCD',
        'MMRRC:000001-UNC'
    ]

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'mmrrc')
        self.strain_hash = {}
        self.id_label_hash = {}
        self.dataset = Dataset(
            'mmrrc', 'Mutant Mouse Regional Resource Centers',
            'https://www.mmrrc.org', None,
            'https://www.mmrrc.org/about/data_download.php')

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)
        fname = '/'.join((self.rawdir, self.files['catalog']['file']))
        st = os.stat(fname)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        # TODO note: can set the data version to what is in the header
        # first line like:
        # This MMRRC catalog data file was generated on 2015-04-22

        self.dataset.setVersion(filedate)

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_phenotype_data(limit)

        logger.info("Finished parsing.")

        return

    def _process_phenotype_data(self, limit):
        """
        NOTE: If a Strain carries more than one mutation,
        then each Mutation description,
        i.e., the set: (
            Mutation Type - Chromosome - Gene Symbol -
            Gene Name - Allele Symbol - Allele Name)
        will require a separate line.

        Note that MMRRC curates phenotypes to alleles,
        even though they distribute only one file with the
        phenotypes appearing to be associated with a strain.

        So, here we process the allele-to-phenotype relationships separately
        from the strain-to-allele relationships.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        fname = '/'.join((self.rawdir, self.files['catalog']['file']))

        self.strain_hash = {}
        self.id_label_hash = {}
        genes_with_no_ids = set()
        stem_cell_class = 'CL:0000034'
        mouse_taxon = 'NCBITaxon:10090'
        geno = Genotype(g)
        with open(fname, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for row in filereader:
                line_counter += 1
                # skip the first 3 lines which are header, etc.
                if line_counter < 4:
                    continue

                (strain_id, strain_label, strain_type_symbol, strain_state,
                 mgi_allele_id, mgi_allele_symbol, mgi_allele_name,
                 mutation_type, chrom, mgi_gene_id, mgi_gene_symbol,
                 mgi_gene_name, sds_url, accepted_date, mp_ids, pubmed_nums,
                 research_areas) = row

                if self.testMode and (strain_id not in self.test_ids) \
                        or mgi_gene_name == 'withdrawn':
                    continue

                # strip off stuff after the dash -
                # is the holding center important?
                # MMRRC:00001-UNC --> MMRRC:00001
                strain_id = re.sub(r'-\w+$', '', strain_id)

                self.id_label_hash[strain_id] = strain_label

                # get the variant or gene to save for later building of
                # the genotype
                if strain_id not in self.strain_hash:
                    self.strain_hash[strain_id] = {'variants': set(),
                                                   'genes': set()}

                # clean up the bad one
                if mgi_allele_id == 'multiple mutation':
                    logger.error("Erroneous gene id: %s", mgi_allele_id)
                    mgi_allele_id = ''

                if mgi_allele_id != '':
                    self.strain_hash[strain_id]['variants'].add(mgi_allele_id)
                    self.id_label_hash[mgi_allele_id] = mgi_allele_symbol

                    # use the following if needing to add the
                    # sequence alteration types
                    # var_type =
                    #   self._get_variant_type_from_abbrev(mutation_type)
                    # make a sequence alteration for this variant locus,
                    # and link the variation type to it
                    # sa_id = '_'+re.sub(r':','',mgi_allele_id)+'SA'
                    # if self.nobnodes:
                    #     sa_id = ':'+sa_id
                    # gu.addIndividualToGraph(g, sa_id, None, var_type)
                    # geno.addSequenceAlterationToVariantLocus(sa_id,
                    #                                          mgi_allele_id)

                # scrub out any spaces
                mgi_gene_id = re.sub(r'\s+', '', mgi_gene_id)
                if mgi_gene_id.strip() != '':
                    if re.match(r'Gene\s*ID:', mgi_gene_id, re.I):
                        mgi_gene_id = re.sub(r'Gene\s*ID:\s*', 'NCBIGene:',
                                             mgi_gene_id)
                    elif not re.match(r'MGI', mgi_gene_id):
                        logger.info("Gene id not recognized: %s", mgi_gene_id)
                        if re.match(r'\d+$', mgi_gene_id):
                            # assume that if it's all numbers, then it's MGI
                            mgi_gene_id = 'MGI:'+str(mgi_gene_id)
                            logger.info("Assuming numerics are MGI.")
                    self.strain_hash[strain_id]['genes'].add(mgi_gene_id)
                    self.id_label_hash[mgi_gene_id] = mgi_gene_symbol

                # catch some errors -
                # some things have gene labels, but no identifiers - report
                if mgi_gene_symbol.strip() != '' and mgi_gene_id == '':
                    logger.error(
                        "Gene label with no identifier for strain %s: %s",
                        strain_id, mgi_gene_symbol)
                    genes_with_no_ids.add(mgi_gene_symbol.strip())
                    # make a temp id for genes that aren't identified
                    # tmp_gene_id = '_'+mgi_gene_symbol
                    # self.id_label_hash[tmp_gene_id] = mgi_gene_symbol
                    # self.strain_hash[strain_id]['genes'].add(tmp_gene_id)

                # split apart the mp ids
                # ataxia [MP:0001393] ,hypoactivity [MP:0001402] ...
                # mp_ids are now a comma delimited list
                # with MP terms in brackets
                phenotype_ids = []
                if mp_ids != '':
                    for i in re.split(r',', mp_ids):
                        i = i.strip()
                        mps = re.search(r'\[(.*)\]', i)
                        if mps is not None:
                            mp_id = mps.group(1).strip()
                            phenotype_ids.append(mp_id)

                # pubmed ids are space delimited
                pubmed_ids = []
                if pubmed_nums.strip() != '':
                    for i in re.split(r'\s+', pubmed_nums):
                        pmid = 'PMID:'+i.strip()
                        pubmed_ids.append(pmid)
                        r = Reference(g, pmid,
                                      Reference.ref_types['journal_article'])
                        r.addRefToGraph()

                # https://www.mmrrc.org/catalog/sds.php?mmrrc_id=00001
                # is a good example of 4 genotype parts

                model.addClassToGraph(mouse_taxon, None)
                if research_areas.strip() == '':
                    research_areas = None
                else:
                    research_areas = 'Research Areas: '+research_areas
                strain_type = mouse_taxon
                if strain_state == 'ES':
                    strain_type = stem_cell_class
                model.addIndividualToGraph(
                    strain_id, strain_label, strain_type,
                    research_areas)  # an inst of mouse??
                model.makeLeader(strain_id)

                # phenotypes are associated with the alleles
                for pid in phenotype_ids:
                    # assume the phenotype label is in the ontology
                    model.addClassToGraph(pid, None)
                    if mgi_allele_id is not None and mgi_allele_id != '':
                        assoc = G2PAssoc(g, self.name, mgi_allele_id, pid,
                                         model.object_properties['has_phenotype'])
                        for p in pubmed_ids:
                            assoc.add_source(p)
                        assoc.add_association_to_graph()
                    else:
                        logger.info("Phenotypes and no allele for %s",
                                    strain_id)

                if not self.testMode and (
                        limit is not None and line_counter > limit):
                    break

            # now that we've collected all of the variant information, build it
            # we don't know their zygosities
            for s in self.strain_hash:
                h = self.strain_hash.get(s)
                variants = h['variants']
                genes = h['genes']
                vl_set = set()
                # make variant loci for each gene
                if len(variants) > 0:
                    for v in variants:
                        vl_id = v
                        vl_symbol = self.id_label_hash[vl_id]
                        geno.addAllele(vl_id, vl_symbol,
                                       geno.genoparts['variant_locus'])
                        vl_set.add(vl_id)
                        if len(variants) == 1 and len(genes) == 1:
                            for gene in genes:
                                geno.addAlleleOfGene(vl_id, gene)
                        else:
                            geno.addAllele(vl_id, vl_symbol)
                else:  # len(vars) == 0
                    # it's just anonymous variants in some gene
                    for gene in genes:
                        vl_id = '_:' + re.sub(r':', '', gene) + '-VL'
                        vl_symbol = self.id_label_hash[gene]+'<?>'
                        self.id_label_hash[vl_id] = vl_symbol
                        geno.addAllele(vl_id, vl_symbol,
                                       geno.genoparts['variant_locus'])
                        geno.addGene(gene, self.id_label_hash[gene])
                        geno.addAlleleOfGene(vl_id, gene)
                        vl_set.add(vl_id)

                # make the vslcs
                vl_list = sorted(vl_set)
                vslc_list = []
                for vl in vl_list:
                    # for unknown zygosity
                    vslc_id = re.sub(r'^_', '', vl)+'U'
                    vslc_id = re.sub(r':', '', vslc_id)
                    vslc_id = '_:' + vslc_id
                    vslc_label = self.id_label_hash[vl] + '/?'
                    self.id_label_hash[vslc_id] = vslc_label
                    vslc_list.append(vslc_id)
                    geno.addPartsToVSLC(
                        vslc_id, vl, None, geno.zygosity['indeterminate'],
                        geno.object_properties['has_alternate_part'], None)
                    model.addIndividualToGraph(
                        vslc_id, vslc_label,
                        geno.genoparts['variant_single_locus_complement'])
                if len(vslc_list) > 0:
                    if len(vslc_list) > 1:
                        gvc_id = '-'.join(vslc_list)
                        gvc_id = re.sub(r'_|:', '', gvc_id)
                        gvc_id = '_:'+gvc_id
                        gvc_label = \
                            '; '.join(self.id_label_hash[v] for v in vslc_list)
                        model.addIndividualToGraph(
                            gvc_id, gvc_label,
                            geno.genoparts['genomic_variation_complement'])
                        for vslc_id in vslc_list:
                            geno.addVSLCtoParent(vslc_id, gvc_id)
                    else:
                        # the GVC == VSLC, so don't have to make an extra piece
                        gvc_id = vslc_list.pop()
                        gvc_label = self.id_label_hash[gvc_id]

                    genotype_label = gvc_label + ' [n.s.]'
                    bkgd_id = \
                        re.sub(r':', '', '-'.join(
                            (geno.genoparts['unspecified_genomic_background'],
                             s)))
                    genotype_id = '-'.join((gvc_id, bkgd_id))
                    bkgd_id = '_:'+bkgd_id
                    geno.addTaxon(mouse_taxon, bkgd_id)
                    geno.addGenomicBackground(
                        bkgd_id, 'unspecified ('+s+')',
                        geno.genoparts['unspecified_genomic_background'],
                        "A placeholder for the " +
                        "unspecified genetic background for "+s)
                    geno.addGenomicBackgroundToGenotype(
                        bkgd_id, genotype_id,
                        geno.genoparts['unspecified_genomic_background'])
                    geno.addParts(
                        gvc_id, genotype_id,
                        geno.object_properties['has_alternate_part'])
                    geno.addGenotype(genotype_id, genotype_label)
                    g.addTriple(
                        s, geno.object_properties['has_genotype'],
                        genotype_id)
                else:
                    # logger.debug(
                    #   "Strain %s is not making a proper genotype.", s)
                    pass

            logger.warning(
                "The following gene symbols did not list identifiers: %s",
                str(sorted(list(genes_with_no_ids))))

        return

    @staticmethod
    def _get_variant_type_from_abbrev(abbrev):
        """
        All variants are generically typed as "sequence_alterations"
        unless otherwise stated.
        :param abbrev:
        :return:

        """
        variant_type = None

        var_dict = {
            'SM': 'SO:0001059',  # spontaneous mutation
            'TM': 'SO:0001059',  # targeted mutation
            'TG': 'SO:xxxxxxx',  # transgenic
            'GT': 'SO:0001059',  # gene trap
            'CI': 'SO:0001059',  # chemically induced mutation
            'RAD': 'SO:0001059',  # radiation induced mutation
            # chromosomal aberration --> chromosomal structure variation
            'CH': 'SO:1000183',
            'RB': 'SO:1000043',  # Robertsonian translocation
            'TL': 'SO:1000048',  # reciprocal translocation
            'TP': 'SO:0000453',  # transposition
            'INV': 'SO:1000036',  # inversion
            'INS': 'SO:0000667',  # insertion
            'DEL': 'SO:0000159',  # deletion
            'DP': 'SO:1000035',  # duplication
            'OTH': 'SO:0001059'  # other
        }
        if abbrev in var_dict:
            variant_type = var_dict[abbrev]
        else:
            logger.warning("Variant type not recognized: %s", abbrev)

        return variant_type

    def getTestSuite(self):
        import unittest
        from tests.test_mmrrc import MMRRCTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MMRRCTestCase)

        return test_suite
