from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Pathway import Pathway
from dipper.models.Dataset import Dataset
import logging
import os
import yaml
import csv

logger = logging.getLogger(__name__)


class Reactome(Source):
    """
    Reactome is a free, open-source, curated and peer reviewed pathway database.
    (http://reactome.org/)
    """
    PANTHER_BASE = "http://www.reactome.org/download/current/"
    files = {
        'ensembl2pathway': {
            'file': 'Ensembl2Reactome.txt',
            'url': PANTHER_BASE+'Ensembl2Reactome.txt'},
    }
    map_files = {
        'eco_map': '../../resources/eco_map.yaml',
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'reactome')
        self.dataset = Dataset(
            'reactome', 'Reactome', 'http://reactome.org/', None,
            'http://reactome.org/pages/about/license-agreement/')


    def fetch(self, is_dl_forced=False):
        """
        Override Source.fetch()
        Fetches resources from reactome using the Reactome.files dictionary
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
        self.load_bindings()
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        file = '/'.join((self.rawdir, self.files['ensembl2pathway']['file']))
        self._parse_ensembl_reactome(file, limit)

        return

    def _parse_ensembl_reactome(self, file=None, limit=None):
        """
        Parse ensembl gene to reactome pathway file
        :param file: file path (not handle)
        :param limit: limit (int, optional) limit the number of rows processed
        :return: None
        """
        eco_map = self._get_eco_map()
        count = 0
        with open(file, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                (gene, pathway_id, pathway_iri, pathway_label,
                 go_ecode, species_name) = row
                count += 1
                self._add_gene_pathway_association(
                    eco_map, gene, pathway_id, pathway_label, go_ecode)

                if limit is not None and count >= limit:
                    break

        return

    def _get_eco_map(self):
        """
        :return: dict
        """
        eco_map = {}
        if os.path.exists(os.path.join(os.path.dirname(__file__),
                                       self.map_files['eco_map'])):
            map_file = open(os.path.join(
                os.path.dirname(__file__), self.map_files['eco_map']), 'r')
            eco_map = yaml.load(map_file)
            map_file.close()
        else:
            logger.warn("IMPC map file not found")

        return eco_map

    def _add_gene_pathway_association(
            self, eco_map, gene, pathway_id, pathway_label, go_ecode):
        pathway = Pathway(self.graph)

        pathway_curie = "REACT:" + pathway_id
        gene_curie = "ENSEMBL:" + gene.strip()
        eco_curie = eco_map[go_ecode]
        pathway.addPathway(pathway_curie, pathway_label)
        pathway.addComponentToPathway(gene_curie, pathway_curie)

        association = Assoc(self.name)
        association.sub = gene_curie
        association.rel = pathway.object_properties['involved_in']
        association.obj = pathway_curie
        association.set_association_id()
        association.add_evidence(eco_curie)
        association.add_association_to_graph(self.graph)
        return




