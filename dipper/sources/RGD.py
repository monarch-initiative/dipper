import logging

from ontobio.io.gafparser import GafParser
from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.models.Reference import Reference

__author__ = 'timputman'

LOG = logging.getLogger(__name__)


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

    def __init__(self, graph_type, are_bnodes_skolemized, skip_stats=False):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            skip_stats=skip_stats,
            name='rgd',
            ingest_title='Rat Genome Database',
            ingest_url='http://rgd.mcw.edu/',
            ingest_logo='https://github.com/monarch-initiative/monarch-ui/blob/master/public/img/sources/source-rgd.png',
            license_url=None,
            data_rights='https://rgd.mcw.edu/wg/disclaimer/',
            # file_handle=None
        )
        self.dataset.set_citation('https://rgd.mcw.edu/wg/citing-rgd/')

    def fetch(self, is_dl_forced=False):
        """
        Override Source.fetch()
        Fetches resources from rat_genome_database via the rat_genome_database ftp site
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
            LOG.info("Only parsing first %d rows", limit)

        rgd_file = '/'.join(
            (self.rawdir, self.files['rat_gene2mammalian_phenotype']['file']))
        # ontobio gafparser implemented here
        p = GafParser()
        assocs = p.parse(open(rgd_file, "r"))

        for i, assoc in enumerate(assocs):
            if 'relation' in assoc.keys():
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

        record['relation']['id'] = self.resolve("has phenotype")
        # define the triple
        gene = record['subject']['id']
        relation = record['relation']['id']
        phenotype = record['object']['id']

        # instantiate the association
        g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=phenotype, pred=relation)

        # add the references
        references = record['evidence']['has_supporting_reference']
        # created RGDRef prefix in curie map to route to proper reference URL in RGD
        references = [
            x.replace('RGD', 'RGDRef') if 'PMID' not in x else x for x in references]

        if len(references) > 0:
            # make first ref in list the source
            g2p_assoc.add_source(identifier=references[0])
            ref_model = Reference(
                self.graph, references[0],
                self.globaltt['publication']
            )
            ref_model.addRefToGraph()

        if len(references) > 1:
            # create equivalent source for any other refs in list
            # This seems to be specific to this source and
            # there could be non-equivalent references in this list
            for ref in references[1:]:
                model.addSameIndividual(sub=references[0], obj=ref)

        # add the date created on
        g2p_assoc.add_date(date=record['date'])
        g2p_assoc.add_evidence(self.resolve(record['evidence']['type']))  # ?set where?
        g2p_assoc.add_association_to_graph()

        return
