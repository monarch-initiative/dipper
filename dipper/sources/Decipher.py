import logging
import re
import io
import csv
from zipfile import ZipFile

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Model import Model
from dipper.models.Reference import Reference
from dipper import config
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.sources.HGNC import HGNC
from dipper.models.Genotype import Genotype

logger = logging.getLogger(__name__)


class Decipher(Source):
    """
    The Decipher group curates and assembles the Development Disorder Genotype
    Phenotype Database (DDG2P) which is a curated list of genes reported to be
    associated with developmental disorders, compiled by clinicians as part of
    the DDD study to facilitate clinical feedback of likely causal variants.

    Beware that the redistribution of this data is a bit unclear from the
    [license](https://decipher.sanger.ac.uk/legal).
    If you intend to distribute this data, be sure to have the appropriate
    licenses in place.

    """

    files = {
        'annot': {
            'file': 'ddg2p.zip',
            'url': 'https://decipher.sanger.ac.uk/files/ddd/ddg2p.zip'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'decipher')

        self.dataset = Dataset(
            'decipher', 'Development Disorder Genotype – Phenotype Database',
            'https://decipher.sanger.ac.uk/', None,
            'https://decipher.sanger.ac.uk/legal')

        if 'test_ids' not in config.get_config() \
                or 'disease' not in config.get_config()['test_ids']:
            logger.warning("not configured with disease test ids.")
            self.test_ids = []
        else:
            self.test_ids = config.get_config()['test_ids']['disease']

        self.g = self.graph
        self.geno = Genotype(self.g)
        self.model = Model(self.g)

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)

        # since there's a dependency on HGNC files; fetch those too

        hgnc = HGNC()
        hgnc.fetch(is_dl_forced)

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True
            self.g = self.testgraph
        else:
            self.g = self.graph

        self.geno = Genotype(self.g)

        # rare disease-phenotype associations
        self._process_ddg2p_annotations(limit)

        logger.info("Finished parsing.")

        return

    def _process_ddg2p_annotations(self, limit):
        """
        The ddg2p annotations associate a gene symbol to an omim disease,
        along with some HPO ids and pubs. The gene symbols come from gencode,
        which in turn come from HGNC official gene symbols.  Therefore,
        we use the HGNC source class to get the id/symbol mapping for
        use in our annotations here.

        According to http://www.gencodegenes.org/faq.html,
        "Gene names are usually HGNC or MGI-approved gene symbols mapped
        to the GENCODE genes by the Ensembl xref pipeline. Sometimes,
        when there is no official gene symbol, the Havana clone-based
        name is used."

        The kind of variation that is linked to a disease is indicated
        (LOF, GOF, CNV, etc) in the source data.
        Here, we create an anonymous variant of the specified gene of
        the indicated type (mapped to the sequence ontology (SO)).

        :param limit:
        :return:

        """

        line_counter = 0
        if self.g is not None:
            g = self.g
        else:
            g = self.graph

        # in order for this to work, we need to map the HGNC id-symbol;
        hgnc = HGNC()
        hgnc_symbol_id_map = hgnc.get_symbol_id_map()

        myzip = ZipFile(
            '/'.join((self.rawdir, self.files['annot']['file'])), 'r')

        # use the ddg2p.txt file
        fname = 'ddg2p.txt'

        unmapped_omim_counter = 0
        unmapped_gene_count = 0
        with myzip.open(fname, 'r') as f:
            f = io.TextIOWrapper(f)
            reader = csv.reader(f, delimiter='\t', quotechar='\"')
            # score_means_by_measure = {}
            # strain_scores_by_measure = {}   # TODO theseare unused
            for row in reader:
                line_counter += 1
                if re.match(r'#', row[0]):   # skip comments
                    continue

                (gencode_gene_name, mode, category, consequence, disease, omim,
                 ddg2p_id, pubmed_ids, hpo_codes) = row

                hgnc_id = hgnc_symbol_id_map.get(gencode_gene_name.strip())
                if hgnc_id is None:
                    logger.error(
                        "Couldn't map the gene symbol %s to HGNC.",
                        gencode_gene_name)
                    unmapped_gene_count += 1
                    continue
                # add the gene
                self.model.addClassToGraph(hgnc_id, gencode_gene_name)

                # TODO make VSLC with the variation
                #   to associate with the disorder
                # TODO use the Inheritance and Mutation consequence
                #   to classify the VSLCs

                allele_id = self.make_allele_by_consequence(
                    consequence, hgnc_id, gencode_gene_name)

                if omim.strip() != '':
                    omim_id = 'OMIM:'+str(omim.strip())
                    # assume this is declared elsewhere in ontology
                    self.model.addClassToGraph(omim_id, None)

                    if category.strip() == 'Confirmed DD gene':
                        rel = self.model.object_properties['has_phenotype']
                    elif category.strip() == 'Probable DD gene':
                        rel = self.model.object_properties['has_phenotype']
                    elif category.strip() == 'Possible DD gene':
                        rel = self.model.object_properties['contributes_to']
                    elif category.strip() == 'Not DD gene':
                        # TODO negative annotation
                        continue
                    assoc = G2PAssoc(g, self.name, allele_id, omim_id)
                    # TODO 'rel' is assigned to but never used

                    for p in re.split(r';', pubmed_ids):
                        p = p.strip()
                        if p != '':
                            pmid = 'PMID:'+str(p)
                            r = Reference(
                                g, pmid, Reference.ref_types['journal_article'])
                            r.addRefToGraph()
                            assoc.add_source(pmid)

                    assoc.add_association_to_graph()
                else:
                    # these are unmapped to a disease id.
                    # note that some match OMIM disease labels
                    # but the identifiers are just not included.
                    # TODO consider mapping to OMIM or DOIDs in other ways
                    logger.warning(
                        "No omim id on line %d\n%s", line_counter, str(row))
                    unmapped_omim_counter += 1

                # TODO hpo phenotypes
                # since the DDG2P file is not documented,
                # I don't know what the HPO annotations are actually about
                # are they about the gene?  the omim disease?  something else?
                # So, we wont create associations until this is clarified

                if not self.testMode and limit is not None \
                        and line_counter > limit:
                    break

        myzip.close()
        logger.warning(
            "gene-disorder associations with no omim id: %d",
            unmapped_omim_counter)
        logger.warning("unmapped gene count: %d", unmapped_gene_count)

        return

    def make_allele_by_consequence(self, consequence, gene_id, gene_symbol):
        """
        Given a "consequence" label that describes a variation type,
        create an anonymous variant of the specified gene as an instance of
        that consequence type.

        :param consequence:
        :param gene_id:
        :param gene_symbol:
        :return: allele_id
        """

        allele_id = None

        # Loss of function : Nonsense, frame-shifting indel,
        #   essential splice site mutation, whole gene deletion or any other
        #   mutation where functional analysis demonstrates clear reduction
        #   or loss of function
        # All missense/in frame : Where all the mutations described in the data
        #   source are either missense or in frame deletions and there is no
        #   evidence favoring either loss-of-function, activating or
        #   dominant negative effect
        # Dominant negative : Mutation within one allele of a gene that creates
        #   a significantly greater deleterious effect on gene product
        #   function than a monoallelic loss of function mutation
        # Activating : Mutation, usually missense that results in
        #   a constitutive functional activation of the gene product
        # Increased gene dosage : Copy number variation that increases
        #   the functional dosage of the gene
        # Cis-regulatory or promotor mutation : Mutation in cis-regulatory
        #   elements that lies outwith the known transcription unit and
        #   promotor of the controlled gene
        # Uncertain : Where the exact nature of the mutation is unclear or
        #   not recorded
        so_type = {                                 # type of variant
            'Loss of function': 'SO:0002054',       # loss of function
            'All missense/in frame': 'SO:0001583',  # missense
            'Dominant negative': 'SO:0002052',      # dominant negative
            'Activating': 'SO:0002053',             # gain of function
            'Increased gene dosage': 'SO:0001742',  # copy number gain
            # regulatory region
            'Cis-regulatory or promotor mutation': 'SO:0001566',
            'Uncertain': 'SO:0001060',              # generic sequence
            '5 or 3UTR mutation': 'SO:0001622',     # UTR
        }

        type_id = so_type.get(consequence)
        if type_id is None:
            logger.warning("Consequence type unmapped: %s", str(consequence))
            type_id = 'SO:0001060'  # sequence variant

        # make the allele
        allele_id = ''.join((gene_id, type_id))
        allele_id = re.sub(r':', '', allele_id)
        allele_id = '_:'+allele_id  # make this a BNode
        allele_label = ' '.join((consequence, 'allele in', gene_symbol))

        self.model.addIndividualToGraph(allele_id, allele_label, type_id)
        self.geno.addAlleleOfGene(allele_id, gene_id)

        return allele_id

    # def getTestSuite(self):
    #     # TODO
    #     import unittest
    #     from tests.test_decipher import DecipherTestCase
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(DecipherTestCase)
    #
    #     return test_suite
