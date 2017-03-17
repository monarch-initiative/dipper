from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Pathway import Pathway
from dipper.models.Dataset import Dataset
import logging
import csv

logger = logging.getLogger(__name__)


class Reactome(Source):
    """
    Reactome is a free, open-source, curated and peer reviewed pathway database.
    (http://reactome.org/)
    """
    REACTOME_BASE = "http://www.reactome.org/download/current/"
    files = {
        'ensembl2pathway': {
            'file': 'Ensembl2Reactome.txt',
            'url': REACTOME_BASE+'Ensembl2Reactome.txt'},
        'chebi2pathway': {
            'file': 'ChEBI2Reactome.txt',
            'url': REACTOME_BASE + 'ChEBI2Reactome.txt'},
    }

    map_files = {
        'eco_map': 'http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt',
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
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        ensembl_file = '/'.join((self.rawdir, self.files['ensembl2pathway']['file']))
        self._parse_reactome_association_file(
            ensembl_file, limit, subject_prefix='ENSEMBL', object_prefix='REACT')
        chebi_file = '/'.join((self.rawdir, self.files['chebi2pathway']['file']))
        self._parse_reactome_association_file(
            chebi_file, limit, subject_prefix='CHEBI', object_prefix='REACT')

        return

    def _parse_reactome_association_file(
            self, file, limit=None, subject_prefix=None, object_prefix=None):
        """
        Parse ensembl gene to reactome pathway file
        :param file: file path (not handle)
        :param limit: limit (int, optional) limit the number of rows processed
        :return: None
        """
        eco_map = Reactome.get_eco_map(Reactome.map_files['eco_map'])
        count = 0
        with open(file, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                (component, pathway_id, pathway_iri, pathway_label,
                 go_ecode, species_name) = row
                count += 1
                self._add_component_pathway_association(
                    eco_map, component, subject_prefix, pathway_id, object_prefix,
                    pathway_label, go_ecode)

                if limit is not None and count >= limit:
                    break

        return

    def _add_component_pathway_association(
            self, eco_map, component, component_prefix, pathway_id,
            pathway_prefix, pathway_label, go_ecode):
        pathway = Pathway(self.graph)

        pathway_curie = "{}:{}".format(pathway_prefix, pathway_id)
        gene_curie = "{}:{}".format(component_prefix, component.strip())
        eco_curie = eco_map[go_ecode]
        pathway.addPathway(pathway_curie, pathway_label)
        pathway.addComponentToPathway(gene_curie, pathway_curie)

        association = Assoc(self.graph, self.name)
        association.sub = gene_curie
        association.rel = pathway.object_properties['involved_in']
        association.obj = pathway_curie
        association.set_association_id()
        association.add_evidence(eco_curie)
        association.add_association_to_graph()
        return




