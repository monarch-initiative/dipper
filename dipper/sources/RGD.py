from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
from dipper.models.Reference import Reference
import csv
import logging

__author__ = 'timputman'

logger = logging.getLogger(__name__)


class RGD(Source):
    """
    Ingest of Rat Genome Database gene to mammalian phenotype gaf file

    """
    RGD_BASE = 'ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/'
    files = {
        'rat_gene2mammalian_phenotype': {
            'file': 'rattus_genes_mp',
            'url': RGD_BASE + 'rattus_genes_mp'},
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'rgd')
        self.dataset = Dataset(
            'rgd', 'RGD', 'http://rgd.mcw.edu/', None, None)

        self.global_terms = Source.open_and_parse_yaml('../../translationtable/global_terms.yaml')

    def fetch(self, is_dl_forced=False):
        """
        Override Source.fetch()
        Fetches resources from rat_genome_database using the rat_genome_database ftp site
        Args:
            :param is_dl_forced (bool): Force download
        Returns:
            :return None
        """
        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        """
        Override Source.parse()
        Args:
            :param limit (int, optional) limit the number of rows processed
        Returns:
            :return None
        """
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        rgd_file = '/'.join((self.rawdir, self.files['rat_gene2mammalian_phenotype']['file']))
        with open(rgd_file, 'r') as file:
            rgd_reader = csv.reader(file, delimiter='\t')
            for x in range(3):
                next(rgd_reader)
            headers = [
                'db',
                'db_object_id',
                'db_object_symbol',
                'qualifier',
                'ontology_id',
                'db_reference',
                'evidence_code',
                'with_or_from',
                'aspect',
                'db_object_name',
                'db_object_synonym',
                'db_object_type',
                'taxon',
                'date',
                'assigned_by',
                'annotation_extension',
                'gene_product_form_id',
            ]

            for i, row in enumerate(rgd_reader):
                assoc = dict(zip(headers, row))
                self.make_association(record=assoc)
                if limit is not None and i > limit:
                    break
        return

    def make_association(self, record):
        """
        contstruct the association
        :param record:
        :return: modeled association of  genotype to mammalian phenotype
        """
        model = Model(self.graph)

        # define the triple
        gene = '{0}:{1}'.format(record['db'], str(record['db_object_id']))
        relation = 'RO:0002200'
        phenotype = record['ontology_id']

        # instantiate the association
        g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=phenotype, pred=relation)

        # add the references
        references = record['db_reference'].split('|')
        # created RGDRef prefix in curie map to route to proper reference URL in RGD
        references = [x.replace('RGD', 'RGDRef') if 'PMID' not in x else x for x in references]

        if len(references) > 0:
            # make first ref in list the source
            g2p_assoc.add_source(identifier=references[0])
            ref_model = Reference(
                self.graph, references[0],
                Reference.ref_types['publication']
            )
            ref_model.addRefToGraph()

        if len(references) > 1:
            # create equivalent source for any other refs in list
            # This seems to be specific to this source and there could be non-equivalent references in this list
            for ref in references[1:]:
                model.addSameIndividual(sub=references[0], obj=ref)

        # add the date created on
        g2p_assoc.add_date(date=record['date'])
        g2p_assoc.add_evidence(self.global_terms[record['evidence_code']])
        g2p_assoc.add_association_to_graph()

        return
