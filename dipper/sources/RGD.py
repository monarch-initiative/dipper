from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.models.Provenance import Provenance
from dipper.models.Dataset import Dataset
from ontobio.io.gafparser import GafParser
from pprint import pprint
import logging

__author__ = 'timputman'

logger = logging.getLogger(__name__)


class RGD(Source):
    RGD_BASE = 'ftp://ftp.rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/'
    files = {
        'rat_gene2mammalian_phenotype': {
            'file': 'rattus_genes_mp',
            'url': RGD_BASE + 'rattus_genes_mp'},
    }

    map_files = {
        'eco_map': {'IEA': 'ECO:0000501',
                    'IAGP': 'ECO:0005613',
                    'IDA': 'ECO:0000314',
                    'IMP': 'ECO:0000315',
                    'IED': 'ECO:0005611',
                    'TAS': 'ECO:0000304',
                    'QTM': 'ECO:0000061',
                    'ISS': 'ECO:0000250',
                    'NAS': 'ECO:0000303',
                    'IPM': 'ECO:0005612',
                    'IEP': 'ECO:0000270'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'rat_genome_database')
        self.dataset = Dataset(
            'rat_genome_database', 'Rat_Genome_Database', 'http://rgd.mcw.edu/', None,
            None)

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
        p = GafParser()
        assocs = p.parse(open(rgd_file, "r"))

        for i, assoc in enumerate(assocs):
            if assoc['date'] is None:
                print(assoc)
            if assoc['relation']['id'] is None:
                assoc['relation']['id'] = 'RO:0002200'
            self.make_association(assoc)
            if limit is not None and i > limit:
                break

        return

    def make_association(self, record):
        model = Model(self.graph)
        provenance_model = Provenance(self.graph)
        redate = record['date'].replace('-', '')
        assertion_bnode = self.make_id("{0}{1}{2}".format(record['subject']['label'],
                                                          record['subject']['id'],
                                                          record['object']['id']
                                                          ), '_')

        provenance_model.add_date_created(prov_type=assertion_bnode, date=redate)

        model.addIndividualToGraph(
            assertion_bnode, None,
            provenance_model.provenance_types['assertion'])

        gene = record['subject']['id']
        relation = record['relation']['id']
        phenotype = record['object']['id']

        g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=phenotype, pred=relation)
        references = record['evidence']['has_supporting_reference']

        if len(references) > 0:
            g2p_assoc.add_source(identifier=references[0])
        if len(references) > 1:
            for ref in references[1:]:
                model.addSameIndividual(sub=references[0], obj=ref)

        g2p_assoc.add_evidence(RGD.map_files['eco_map'][record['evidence']['type']])
        g2p_assoc.add_association_to_graph()

        return


example_assoc = {'aspect': 'N',
                 'date': '2006-10-17',
                 'evidence': {'has_supporting_reference': ['RGD:1581602', 'PMID:16368876'],
                              'type': 'IAGP',
                              'with_support_from': []},
                 'negated': False,
                 'object': {'id': 'MP:0005211', 'taxon': 'NCBITaxon:10116'},
                 'provided_by': 'RGD',
                 'qualifiers': [],
                 'relation': {'id': None},
                 'source_line': 'RGD\t621503\tKcnq1\t\tMP:0005211\tRGD:1581602|PMID:16368876\t'
                                'IAGP\t\tN\tpotassium voltage-gated channel subfamily Q member '
                                '1\t\tgene\ttaxon:10116\t20061017\tRGD\t\t\n',
                 'subject': {'fullname': 'potassium voltage-gated channel subfamily Q member 1',
                             'id': 'RGD:621503',
                             'label': 'Kcnq1',
                             'synonyms': [],
                             'taxon': {'id': 'NCBITaxon:10116'},
                             'type': 'gene'},
                 'subject_extensions': [{'filler': '\n', 'property': 'isoform'}]}
