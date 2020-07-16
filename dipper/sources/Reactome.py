import logging
import csv
import yaml
from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Pathway import Pathway

LOG = logging.getLogger(__name__)


class Reactome(Source):
    """
    Reactome is a free, open-source, curated and peer reviewed pathway database.
    (http://reactome.org/)
    """
    REACTOME_BASE = "http://www.reactome.org/download/current/"
    files = {
        'ensembl2pathway': {
            'file': 'Ensembl2Reactome.txt',
            'url': REACTOME_BASE + 'Ensembl2Reactome.txt',
            'columns': [          # e.g.
                'component',      # 1:C34F11.9e
                'pathway_id',     # 2:R-CEL-201688
                'pathway_iri',    # 3:https://reactome.org/PathwayBrowser/#/R-CEL-123
                'pathway_label',  # 4:WNT mediated activation of DVL
                'go_ecode',       # 5:IEA
                'species_nam',    # 6:Caenorhabditis elegans
            ]
        },
        'chebi2pathway': {
            'file': 'ChEBI2Reactome.txt',
            'url': REACTOME_BASE + 'ChEBI2Reactome.txt',
            'columns': [          # e.g.
                'component',      # 1:10033
                'pathway_id',     # 2:R-BTA-6806664
                'pathway_iri',    # 3:https://reactome.org/PathwayBrowser/#/R-BTA-123
                'pathway_label',  # 4:Metabolism of vitamin K
                'go_ecode',       # 5:IEA
                'species_nam',    # 6:Bos taurus
            ]
        },
        'gaf-eco-mapping': {
            'file': 'gaf-eco-mapping.yaml',
            'url': '/'.join((Source.DIPPERCACHE, 'reactome', 'gaf-eco-mapping.yaml')),
        }
    }

    def __init__(
            self,
            graph_type,
            are_bnodes_skolemized,
            data_release_version=None
    ):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='reactome',
            ingest_title='Reactome',
            ingest_url='http://reactome.org/',
            ingest_logo='source-reactome.png',
            license_url=None,
            data_rights='https://reactome.org/license/'
            # file_handle=None
        )
        # gaf evidence code mapping is built in parse(), after the file is fetched.
        self.gaf_eco = {}

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

    def parse(self, limit=None):
        """
        Override Source.parse()
        Args:
            :param limit (int, optional) limit the number of rows processed
        Returns:
            :return None
        """
        src_key = 'gaf-eco-mapping'
        yamlfile = '/'.join((self.rawdir, self.files[src_key]['file']))
        with open(yamlfile, 'r') as yfh:
            self.gaf_eco = yaml.safe_load(yfh)

        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        self._parse_reactome_association_file(
            'ensembl2pathway', limit, subject_prefix='ENSEMBL', object_prefix='REACT'
        )

        self._parse_reactome_association_file(
            'chebi2pathway', limit, subject_prefix='CHEBI', object_prefix='REACT'
        )

    def _parse_reactome_association_file(
            self,
            src_key,
            limit=None,
            subject_prefix=None,
            object_prefix=None
    ):
        """
        Parse ensembl gene to reactome pathway file
        :param file: file path (not handle)
        :param limit: limit (int, optional) limit the number of rows processed
        :return: None
        """
        src_file = '/'.join((self.rawdir, self.files[src_key]['file']))
        col = self.files[src_key]['columns']

        with open(src_file, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                component = row[col.index('component')].strip()
                pathway_id = row[col.index('pathway_id')].strip()
                # pathway_iri = row[col.index('pathway_iri')]
                pathway_label = row[col.index('pathway_label')].strip()
                go_ecode = row[col.index('go_ecode')].strip()
                # species_name = row[col.index('species_name')]

                gene_curie = ':'.join((subject_prefix, component))
                pathway_curie = ':'.join((object_prefix, pathway_id))

                eco_curie = None
                if go_ecode in self.gaf_eco:
                    eco_curie = self.gaf_eco[go_ecode]
                else:
                    LOG.error(
                        'Evidence code %s not found in %s', go_ecode, str(self.gaf_eco))

                self._add_component_pathway_association(
                   gene_curie, pathway_curie, pathway_label, eco_curie)

                if limit is not None and reader.line_num >= limit:
                    break

    def _add_component_pathway_association(
            self, gene_curie, pathway_curie, pathway_label, eco_curie
    ):

        pathway = Pathway(self.graph)
        pathway.addPathway(pathway_curie, pathway_label)
        pathway.addComponentToPathway(gene_curie, pathway_curie)

        association = Assoc(self.graph, self.name)
        association.sub = gene_curie
        association.rel = self.globaltt['involved in']
        association.obj = pathway_curie
        association.set_association_id()
        association.add_evidence(eco_curie)
        association.add_association_to_graph()
