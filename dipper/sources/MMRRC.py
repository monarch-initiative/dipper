import csv
import re
import os
from datetime import datetime
from stat import ST_CTIME
import logging

from dipper.sources.Source import Source
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Reference import Reference
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model

LOG = logging.getLogger(__name__)


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
            'url': 'https://www.mmrrc.org/about/mmrrc_catalog_data.csv',
            'columns': [
                'STRAIN/STOCK_ID',
                'STRAIN/STOCK_DESIGNATION',
                'STRAIN_TYPE',
                'STATE',
                'MGI_ALLELE_ACCESSION_ID',
                'ALLELE_SYMBOL',
                'ALLELE_NAME',
                'MUTATION_TYPE',
                'CHROMOSOME',
                'MGI_GENE_ACCESSION_ID',
                'GENE_SYMBOL',
                'GENE_NAME',
                'SDS_URL',
                'ACCEPTED_DATE',
                'MPT_IDS',
                'PUBMED_IDS',
                'RESEARCH_AREAS',
            ]
        },
    }

    test_ids = [  # move to resources/?
        'MMRRC:037507-MU', 'MMRRC:041175-UCD', 'MMRRC:036933-UNC', 'MMRRC:037884-UCD',
        'MMRRC:000255-MU', 'MMRRC:037372-UCD', 'MMRRC:000001-UNC'
    ]

    def __init__(self, graph_type, are_bnodes_skolemized, skip_stats=False):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            skip_stats=skip_stats,
            name='mmrrc',
            ingest_title='Mutant Mouse Regional Resource Centers',
            ingest_url='https://www.mmrrc.org',
            ingest_logo='https://github.com/monarch-initiative/monarch-ui/blob/master/public/img/sources/source-mmrrc.png',
            # license_url=None,
            data_rights='https://www.mmrrc.org/about/data_download.php'
            # file_handle=None
        )
        self.strain_hash = {}
        self.id_label_hash = {}
        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)
        # TODO note: can set the data version to what is in the header
        # first line like:
        # This MMRRC catalog data file was generated on 2015-04-22

        return

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self._process_phenotype_data(limit)

        LOG.info("Finished parsing.")

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

        src_key = 'catalog'
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        fname = '/'.join((self.rawdir, self.files[src_key]['file']))

        self.strain_hash = {}
        self.id_label_hash = {}
        genes_with_no_ids = set()
        stem_cell_class = self.globaltt['stem cell']
        mouse_taxon = self.globaltt['Mus musculus']
        geno = Genotype(graph)
        with open(fname, 'r', encoding="utf8") as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            # This MMRRC catalog data file was generated on YYYY-MM-DD
            # insert or check date w/dataset
            line = next(reader)
            # gen_date = line[-10:]
            line = next(reader)
            col = self.files['catalog']['columns']
            if col != line:
                LOG.error(
                    '%s\nExpected Headers:\t%s\nRecived Headers:\t%s\n',
                    src_key, col, line)
                LOG.info(set(col) - set(line))

            line = next(reader)
            if line != []:
                LOG.warning('Expected third line to be blank. got "%s" instead', line)

            for row in reader:
                strain_id = row[col.index('STRAIN/STOCK_ID')].strip()
                strain_label = row[col.index('STRAIN/STOCK_DESIGNATION')]
                # strain_type_symbol = row[col.index('STRAIN_TYPE')]
                strain_state = row[col.index('STATE')]
                mgi_allele_id = row[col.index('MGI_ALLELE_ACCESSION_ID')].strip()
                mgi_allele_symbol = row[col.index('ALLELE_SYMBOL')]
                # mgi_allele_name = row[col.index('ALLELE_NAME')]
                # mutation_type = row[col.index('MUTATION_TYPE')]
                # chrom = row[col.index('CHROMOSOME')]
                mgi_gene_id = row[col.index('MGI_GENE_ACCESSION_ID')].strip()
                mgi_gene_symbol = row[col.index('GENE_SYMBOL')].strip()
                mgi_gene_name = row[col.index('GENE_NAME')]
                # sds_url = row[col.index('SDS_URL')]
                # accepted_date = row[col.index('ACCEPTED_DATE')]
                mpt_ids = row[col.index('MPT_IDS')].strip()
                pubmed_nums = row[col.index('PUBMED_IDS')].strip()
                research_areas = row[col.index('RESEARCH_AREAS')].strip()

                if self.test_mode and (strain_id not in self.test_ids) \
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
                    self.strain_hash[strain_id] = {
                        'variants': set(), 'genes': set()}

                # flag bad ones
                if mgi_allele_id[:4] != 'MGI:' and mgi_allele_id != '':
                    LOG.error("Erroneous MGI allele id: %s", mgi_allele_id)
                    if mgi_allele_id[:3] == 'MG:':
                        mgi_allele_id = 'MGI:' + mgi_allele_id[3:]
                    else:
                        mgi_allele_id = ''

                if mgi_allele_id != '':
                    self.strain_hash[strain_id]['variants'].add(mgi_allele_id)
                    self.id_label_hash[mgi_allele_id] = mgi_allele_symbol

                    # use the following if needing to add the sequence alteration types
                    # var_type = self.localtt[mutation_type]
                    # make a sequence alteration for this variant locus,
                    # and link the variation type to it
                    # sa_id = '_'+re.sub(r':','',mgi_allele_id)+'SA'
                    # if self.nobnodes:
                    #     sa_id = ':'+sa_id
                    # gu.addIndividualToGraph(g, sa_id, None, var_type)
                    # geno.addSequenceAlterationToVariantLocus(sa_id, mgi_allele_id)

                # scrub out any spaces, fix known issues
                mgi_gene_id = re.sub(r'\s+', '', mgi_gene_id)
                if mgi_gene_id == 'NULL':
                    mgi_gene_id = ''
                elif mgi_gene_id[:7] == 'GeneID:':
                    mgi_gene_id = 'NCBIGene:' + mgi_gene_id[7:]

                if mgi_gene_id != '':
                    [curie, localid] = mgi_gene_id.split(':')
                    if curie not in ['MGI', 'NCBIGene']:
                        LOG.info("MGI Gene id not recognized: %s", mgi_gene_id)
                    self.strain_hash[strain_id]['genes'].add(mgi_gene_id)
                    self.id_label_hash[mgi_gene_id] = mgi_gene_symbol

                # catch some errors - too many. report summary at the end
                # some things have gene labels, but no identifiers - report
                if mgi_gene_symbol != '' and mgi_gene_id == '':
                    # LOG.error(
                    #    "Gene label with no MGI identifier for strain %s: %s",
                    #    strain_id, mgi_gene_symbol)
                    genes_with_no_ids.add(mgi_gene_symbol)
                    # make a temp id for genes that aren't identified ... err wow.
                    # tmp_gene_id = '_' + mgi_gene_symbol
                    # self.id_label_hash[tmp_gene_id.strip()] = mgi_gene_symbol
                    # self.strain_hash[strain_id]['genes'].add(tmp_gene_id)

                # split apart the mp ids
                # ataxia [MP:0001393] ,hypoactivity [MP:0001402] ...
                # mpt_ids are a comma delimited list
                # labels with MP terms following in brackets
                phenotype_ids = []
                if mpt_ids != '':
                    for lb_mp in mpt_ids.split(r','):
                        lb_mp = lb_mp.strip()
                        if lb_mp[-1:] == ']' and lb_mp[-12:-8] == '[MP:':
                            phenotype_ids.append(lb_mp[-11:-2])

                # pubmed ids are space delimited
                pubmed_ids = []
                if pubmed_nums != '':
                    for pm_num in re.split(r'\s+', pubmed_nums):
                        pmid = 'PMID:' + pm_num.strip()
                        pubmed_ids.append(pmid)
                        ref = Reference(graph, pmid, self.globaltt['journal article'])
                        ref.addRefToGraph()

                # https://www.mmrrc.org/catalog/sds.php?mmrrc_id=00001
                # is a good example of 4 genotype parts

                model.addClassToGraph(mouse_taxon, None)
                if research_areas == '':
                    research_areas = None
                else:
                    research_areas = 'Research Areas: ' + research_areas
                strain_type = mouse_taxon
                if strain_state == 'ES':
                    strain_type = stem_cell_class
                model.addIndividualToGraph(   # an inst of mouse??
                    strain_id, strain_label, strain_type, research_areas)
                model.makeLeader(strain_id)

                # phenotypes are associated with the alleles
                for pid in phenotype_ids:
                    # assume the phenotype label is in some ontology
                    model.addClassToGraph(pid, None)
                    if mgi_allele_id is not None and mgi_allele_id != '':
                        assoc = G2PAssoc(
                            graph, self.name, mgi_allele_id, pid,
                            self.globaltt['has phenotype'])
                        for p in pubmed_ids:
                            assoc.add_source(p)
                        assoc.add_association_to_graph()
                    else:
                        LOG.info("Phenotypes and no allele for %s", strain_id)

                if not self.test_mode and (
                        limit is not None and reader.line_num > limit):
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
                    for var in variants:
                        vl_id = var.strip()
                        vl_symbol = self.id_label_hash[vl_id]
                        geno.addAllele(
                            vl_id, vl_symbol, self.globaltt['variant_locus'])
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
                        geno.addAllele(
                            vl_id, vl_symbol, self.globaltt['variant_locus'])
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
                        vslc_id, vl, None, self.globaltt['indeterminate'],
                        self.globaltt['has_variant_part'], None)
                    model.addIndividualToGraph(
                        vslc_id, vslc_label,
                        self.globaltt['variant single locus complement'])
                if len(vslc_list) > 0:
                    if len(vslc_list) > 1:
                        gvc_id = '-'.join(vslc_list)
                        gvc_id = re.sub(r'_|:', '', gvc_id)
                        gvc_id = '_:'+gvc_id
                        gvc_label = '; '.join(self.id_label_hash[v] for v in vslc_list)
                        model.addIndividualToGraph(
                            gvc_id, gvc_label,
                            self.globaltt['genomic_variation_complement'])
                        for vslc_id in vslc_list:
                            geno.addVSLCtoParent(vslc_id, gvc_id)
                    else:
                        # the GVC == VSLC, so don't have to make an extra piece
                        gvc_id = vslc_list.pop()
                        gvc_label = self.id_label_hash[gvc_id]

                    genotype_label = gvc_label + ' [n.s.]'
                    bkgd_id = re.sub(
                        r':', '', '-'.join((
                            self.globaltt['unspecified_genomic_background'], s)))
                    genotype_id = '-'.join((gvc_id, bkgd_id))
                    bkgd_id = '_:' + bkgd_id
                    geno.addTaxon(mouse_taxon, bkgd_id)
                    geno.addGenomicBackground(
                        bkgd_id, 'unspecified (' + s + ')',
                        self.globaltt['unspecified_genomic_background'],
                        "A placeholder for the unspecified genetic background for " + s)
                    geno.addGenomicBackgroundToGenotype(
                        bkgd_id, genotype_id,
                        self.globaltt['unspecified_genomic_background'])
                    geno.addParts(
                        gvc_id, genotype_id, self.globaltt['has_variant_part'])
                    geno.addGenotype(genotype_id, genotype_label)
                    graph.addTriple(
                        s, self.globaltt['has_genotype'], genotype_id)
                else:
                    # LOG.debug(
                    #   "Strain %s is not making a proper genotype.", s)
                    pass

            LOG.warning(
                "The following gene symbols did not list identifiers: %s",
                str(sorted(list(genes_with_no_ids))))
            LOG.error(
                '%i symbols given are missing their gene identifiers',
                len(genes_with_no_ids))

        return

    def getTestSuite(self):
        import unittest
        from tests.test_mmrrc import MMRRCTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MMRRCTestCase)

        return test_suite
