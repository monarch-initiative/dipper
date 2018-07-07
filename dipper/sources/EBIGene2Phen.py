import logging
import csv
import gzip
import requests
from contextlib import closing
from typing import List, Optional

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Model import Model
from dipper.models.Genotype import Genotype


logger = logging.getLogger(__name__)


class EBIGene2Phen(Source):
    """
    From EBI:
    The gene2phenotype dataset (G2P) integrates data on genes,
    variants and phenotypes for example relating to developmental
    disorders. It is constructed entirely from published literature,
    and is primarily an inclusion list to allow targeted filtering of
    genome-wide data for diagnostic purposes. The dataset was compiled
    with respect to published genes, and annotated with types of disease-
    causing gene variants. Each row of the dataset associates a gene with
    a disease phenotype via an evidence level, inheritance mechanism and
    mutation consequence. Some genes therefore appear in the database
    more than once, where different genetic mechanisms result in
    different phenotypes.

    Disclaimer: https://www.ebi.ac.uk/gene2phenotype/disclaimer
    Terms of Use: https://www.ebi.ac.uk/about/terms-of-use#general
    Documentation: https://www.ebi.ac.uk/gene2phenotype/documentation
                   https://www.clinicalgenome.org/site/assets/files/
                   2757/fitzpatrick_ddg2p.pdf

    This script operates on the Developmental Disorders (DDG2P.csv) file
    In the future we may update to include the cancer gene disease pairs
    in the  CancerG2P.csv file
    """

    EBI_BASE = "https://www.ebi.ac.uk/gene2phenotype/downloads/"

    files = {
        'developmental_disorders': {
            'file': 'DDG2P.csv.gz',
            'url': EBI_BASE + 'DDG2P.csv.gz'}
    }

    map_files = {
        # map for unmapped disease labels to MONDO ids
        'mondo_map': 'https://data.monarchinitiative.org/dipper/'
                     'cache/unmapped_ebi_diseases.tsv'
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'ebi_g2p')

        self.dataset = Dataset(
            'ebi',
            'EBI',
            'https://www.ebi.ac.uk/gene2phenotype',
            None,
            None,
            'https://www.ebi.ac.uk/about/terms-of-use#general',
        )

        # Load mondo map file
        self.mondo_map = {}
        map_url = self.map_files['mondo_map']
        with closing(requests.get(map_url, stream=True)) as map_fh:
            for line in map_fh.iter_lines():
                row = line.decode('utf-8').split('\t')
                self.mondo_map[row[0]] = row[1]

        self.global_terms = self.open_and_parse_yaml('../../translationtable/global_terms.yaml')
        self.translation_table = self.open_and_parse_yaml('../../translationtable/ebi.yaml')

        return

    def fetch(self, is_dl_forced: bool=False):
        """
        Fetch DDG2P.csv.gz and check headers to see
        if it has been updated

        :param is_dl_forced: {bool}
        :return: None
        """
        self.get_files(is_dl_forced)

        return

    def parse(self, limit: Optional[int]=None):
        """
        Here we parse each row of the gene to phenotype file

        We create anonymous variants along with their attributes
        (allelic requirement, functional consequence)
        and connect these to genes and diseases

        genes are connected to variants via
        global_terms['has_affected_locus']

        variants are connected to attributes via:
        global_terms['has_allelic_requirement']
        global_terms['has_functional_consequence']

        variants are connected to disease based on
        mappings to the DDD category column,
        see the translationtable specific to this source
        for mappings

        For cases where there are no disease OMIM id,
        we either use a disease cache file with mappings
        to MONDO that has been manually curated

        :param limit: {int} number of rows to parse
        :return: None
        """
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")
        row_count = 0
        file_path = '/'.join((self.rawdir,
                              self.files['developmental_disorders']['file']))

        with gzip.open(file_path, 'rt') as csvfile:
            reader = csv.reader(csvfile)
            next(reader) # header
            for row in reader:
                if limit is None or row_count <= (limit + 1):
                    self._add_gene_disease(row)
                else:
                    break
                row_count += 1

        logger.info("Done parsing.")

        return

    def _add_gene_disease(self, row:List):
        """
        Parse and add gene variant disease model
        Model building happens in _build_gene_disease_model

        :param row {List}: single row from DDG2P.csv
        :return: None
        """
        if len(row) != 14:
            raise ValueError("Unexpected number of fields "
                             "for row {}".format(row))

        (
            gene_symbol,
            gene_omim_id,
            disease_label,
            disease_omim_id,
            g2p_relation_label,
            allelic_requirement,
            mutation_consequence,
            phenotypes,
            organ_specificity_list,
            pmids,
            panel,
            prev_symbols,
            hgnc_id,
            entry_date,
        ) = row
        variant_label = "variant of {}".format(gene_symbol)
        if disease_omim_id == 'No disease mim':
            # check if we've manually curated
            if disease_label in self.mondo_map:
                disease_id = self.mondo_map[disease_label]
            else:
                return # sorry for this
        else:
            disease_id = 'OMIM:'+disease_omim_id

        hgnc_curie = 'HGNC:' + hgnc_id
        relation_global_term = self.translation_table[g2p_relation_label]
        relation_curie = self.global_terms[relation_global_term]
        consequence_predicate = None

        if mutation_consequence != 'uncertain' and mutation_consequence != '':
            consequence = self._get_consequence_predicate(mutation_consequence)
            consequence_relation = self.global_terms[consequence]
            consequence_curie = self.global_terms[
                self.translation_table[mutation_consequence]
            ]
            variant_label = "{} {}".format(mutation_consequence, variant_label)
        else:
            consequence_relation = None
            consequence_curie = None

        if allelic_requirement != '':
            requirement_curie = self.global_terms[
                self.translation_table[allelic_requirement]
            ]
            variant_label = "{} {}".format(allelic_requirement, variant_label)
        else:
            requirement_curie = None

        if pmids != '':
            pmid_list = ['PMID:'+pmid for pmid in pmids.split(';')]
        else:
            pmid_list = []

        # build the model
        # Should we build a reusable object and/or tuple that
        # could be passed to a more general model builder for
        # this and orphanet (and maybe clinvar)
        self._build_gene_disease_model(
            hgnc_curie,
            relation_curie,
            disease_id,
            variant_label,
            consequence_relation,
            consequence_curie,
            requirement_curie,
            pmid_list
        )
        return

    def _build_gene_disease_model(
            self,
            gene_id,
            relation_id,
            disease_id,
            variant_label,
            consequence_predicate=None,
            consequence_id=None,
            allelic_requirement=None,
            pmids=None):
        """
        Builds gene variant disease model

        :return: None
        """
        model = Model(self.graph)
        geno = Genotype(self.graph)

        pmids = [] if pmids is None else pmids

        is_variant = False
        variant_or_gene = gene_id

        variant_id_string = variant_label + 'EBI'
        if allelic_requirement is not None:
            variant_id_string = variant_label + allelic_requirement + 'EBI'

        variant_bnode = self.make_id(variant_id_string, "_")

        if consequence_predicate is not None \
                and consequence_id is not None:
            is_variant = True
            model.addTriple(variant_bnode,
                            consequence_predicate,
                            consequence_id)
            # Hack to add labels to terms that
            # don't exist in an ontology
            if consequence_id.startswith(':'):
                model.addLabel(consequence_id,
                               consequence_id.strip(':').replace('_', ' '))

        if allelic_requirement is not None:
            is_variant = True
            model.addTriple(variant_bnode,
                            self.global_terms['has_allelic_requirement'],
                            allelic_requirement)
            if allelic_requirement.startswith(':'):
                model.addLabel(allelic_requirement,
                               allelic_requirement.strip(':').replace('_', ' '))
        if is_variant:
            variant_or_gene = variant_bnode
            # Typically we would type the variant using the
            # molecular consequence, but these are not specific
            # enough for us to make mappings (see translation table)
            model.addIndividualToGraph(variant_bnode,
                                       variant_label,
                                       geno.genoparts['variant_locus'])
            geno.addAffectedLocus(variant_bnode, gene_id)
            model.addBlankNodeAnnotation(variant_bnode)

        assoc = G2PAssoc(
            self.graph, self.name, variant_or_gene, disease_id, relation_id)
        assoc.source = pmids
        assoc.add_association_to_graph()

        return

    @staticmethod
    def _get_consequence_predicate(consequence):
        consequence_map = {
            'has_molecular_consequence' : [
                '5_prime or 3_prime UTR mutation',
                'all missense/in frame',
                'cis-regulatory or promotor mutation',
                'part of contiguous gene duplication'
            ],
            'has_functional_consequence': [
                'activating',
                'dominant negative',
                'increased gene dosage',
                'loss of function'
            ]
        }
        consequence_type = 'uncertain'
        for typ, typ_list in consequence_map.items():
            if consequence in typ_list:
                consequence_type = typ

        return consequence_type
