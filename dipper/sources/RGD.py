from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.models.Provenance import Provenance
from dipper.models.Dataset import Dataset
from ontobio.io.gafparser import GafParser
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
        super().__init__(graph_type, are_bnodes_skolemized, 'rat_genome_database')
        self.dataset = Dataset(
            'rat_genome_database', 'Rat_Genome_Database', 'http://rgd.mcw.edu/', None,
            None)

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

        # ontobio gafparser implemented here
        p = GafParser()
        assocs = p.parse(open(rgd_file, "r"))

        for i, assoc in enumerate(assocs):
            assoc['relation']['id'] = 'RO:0002200'
            self.make_association(assoc)
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
        provenance_model = Provenance(self.graph)
        redate = record['date'].replace('-', '')

        # date created is currently modeled as assertion but this is up for review
        assertion_bnode = self.make_id("{0}{1}{2}".format(record['subject']['label'],
                                                          record['subject']['id'],
                                                          record['object']['id']
                                                          ), '_')

        provenance_model.add_date_created(prov_type=assertion_bnode, date=redate)

        model.addIndividualToGraph(
            assertion_bnode, None,
            provenance_model.provenance_types['assertion'])

        # define the triple
        gene = record['subject']['id']
        relation = record['relation']['id']
        phenotype = record['object']['id']

        g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=phenotype, pred=relation)
        references = record['evidence']['has_supporting_reference']

        if len(references) > 0:
            # make first ref in list the source
            g2p_assoc.add_source(identifier=references[0])
        if len(references) > 1:
            # create equivalent source for any other refs in list
            # This seems to be specific to this source and there could be non-equivalent references in this list
            for ref in references[1:]:
                model.addSameIndividual(sub=references[0], obj=ref)

        g2p_assoc.add_evidence(self.global_terms[record['evidence']['type']])
        g2p_assoc.add_association_to_graph()

        return


