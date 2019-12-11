import logging


from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.assoc.Association import Assoc

LOG = logging.getLogger(__name__)


class Template(Source):
    """

    """

    files = {
        'src_key': {
            'file': '',
            'url': '',
            'path': '',

            'columns': [
               'one',
               'two',
               'three',
            ]
        }
    }

    def __init__(
            self,
            graph_type='streamed_graph',
            are_bnodes_skolemized,
            data_release_version=None,
            tax_ids=None,
            version=None):
        """
        :param tax_ids: [str,], List of NCBI taxon  identifiers
        :return:
        """
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='bgee',
            ingest_title='Bgee Gene expression data in animals',
            ingest_url='http://bgee.org/',
            ingest_logo='source-bgee.png',
            # license_url=None,
            data_rights='https://bgee.org/?page=about'
            # file_handle=None
        )
        self.default_taxa = {
            x: self.globaltt[x].split(':')[-1] for x in self.default_species}
        # names for logging
        self.txid_name = {v: k for k, v in self.default_taxa.items()}

        if tax_ids is None:
            tax_ids = self.default_taxa.values()
        self.tax_ids = [str(x) for x in tax_ids]  # incase they were passed in
        LOG.info(
            "Filtering on tax_ids %s",
            [{t: self.txid_name[t]} for t in self.tax_ids])

        if version is None:
            self.version = 'current'
        else:
            self.version = version

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced: boolean, force download
        :return:
        """


    def parse(self):
        """
        :return: None
        """
