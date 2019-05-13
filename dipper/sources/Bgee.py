import logging
import os
import re
import time
import ftplib
import gzip
from datetime import datetime
from stat import ST_SIZE

import pandas as pd
from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.assoc.Association import Assoc


LOG = logging.getLogger(__name__)
BGEE_FTP = 'ftp.bgee.org'


class Bgee(Source):
    """
    Bgee is a database to retrieve and compare gene expression
    patterns between animal species.

    Bgee first maps heterogeneous expression data (currently RNA-Seq,
    Affymetrix, in situ hybridization, and EST data) to anatomy and
    development of different species.

    Then, in order to perform automated cross species comparisons,
    homology relationships across anatomies, and comparison criteria
    between developmental stages, are designed.
    """

    default_species = [
        "Cavia porcellus",                  # guinea pig
        "Mus musculus",                     # mouse
        "Rattus norvegicus",                # rat
        "Monodelphis domestica",            # gray short-tailed opossum
        "Anolis carolinensis",              # lizard
        "Caenorhabditis elegans",           # worm
        "Drosophila melanogaster",          # fly
        "Danio rerio",                      # zebrafish
        "Xenopus (Silurana) tropicalis",    # frog
        "Gallus gallus",                    # chicken
        "Ornithorhynchus anatinus",         # platypus
        "Erinaceus europaeus",              # hedgehog
        "Macaca mulatta",                   # monkey (rhesus macaque)
        "Gorilla gorilla",                  # gorilla
        "Pan paniscus",                     # bonobo
        "Pan troglodytes",                  # chimp
        "Homo sapiens",                     # corporal people
        "Canis lupus familiaris",           # dog  "Canis familiaris"
        "Felis catus",                      # cat
        "Equus caballus",                   # horse
        "Sus scrofa",                       # pig
        "Bos taurus",                       # cow
        "Oryctolagus cuniculus",            # rabbit
        # 7217	Drosophila_ananassae
        # 7230	Drosophila_mojavensis
        # 7237	Drosophila_pseudoobscura
        # 7240	Drosophila_simulans
        # 7244	Drosophila_virilis
        # 7245	Drosophila_yakuba
    ]
    files = {
        'anat_entity': {
            'path': '/download/ranks/anat_entity/',
            'pattern': re.compile(r'^[0-9]+_anat_entity_all_data_.*.tsv.gz'),
            'columns': [
                'Ensembl gene ID',
                'gene name',
                'anatomical entity ID',     # uberon
                'anatomical entity name',
                'rank score',
                'XRefs to BTO',
            ]
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None, version=None):
        """
        :param tax_ids: [str,], List of NCBI taxon  identifiers
        :return:
        """
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'bgee',
            ingest_title='Bgee Gene expression data in animals',
            ingest_url='http://bgee.org/',
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

        (files_to_download, ftp) = self._get_file_list(
            self.files['anat_entity']['path'],
            self.files['anat_entity']['pattern'])

        LOG.info(
            'Will Check \n%s\nfrom %s',
            '\n'.join(list(files_to_download)), ftp.getwelcome())

        for dlname in files_to_download:
            localfile = '/'.join((self.rawdir, dlname))
            info = ftp.sendcmd("MLST {}".format(dlname))   # fetch remote file stats
            info = info.split('\n')[1].strip()              # drop pre & post script
            info = info.split(';')                          # partition fields
            info = [item.strip() for item in info[:-1]]     # cleanup an drop final name
            info = [item.split('=') for item in info]       # make pairs
            info = {item[0]: item[1] for item in info}      # transform list to dict

            LOG.info(
                '%s\n'
                'Remote File Size: %i\n'
                'Remote timestamp: %s',
                dlname, int(info['size']),
                self._convert_ftp_time_to_iso(info['modify']))

            if not os.path.exists(localfile) or is_dl_forced or \
                    self.checkIfRemoteIsNewer(
                            localfile, int(info['size']), info['modify']):

                LOG.info("Fetching %s", dlname)
                LOG.info("Writing to %s", localfile)

                ftp.retrbinary('RETR {}'.format(dlname), open(localfile, 'wb').write)
                remote_dt = Bgee._convert_ftp_time_to_iso(info['modify'])
                os.utime(
                    localfile,
                    (time.mktime(remote_dt.timetuple()),
                     time.mktime(remote_dt.timetuple())))
        ftp.quit()

    def parse(self, limit=None):
        """
        Given the input taxa, expects files in the raw directory
        with the name {tax_id}_anat_entity_all_data_Pan_troglodytes.tsv.zip

        :param limit: int Limit to top ranked anatomy associations per group
        :return: None
        """

        files_to_download, ftp = self._get_file_list(
            self.files['anat_entity']['path'],
            self.files['anat_entity']['pattern'], None)
        for dlname in files_to_download:
            localfile = '/'.join((self.rawdir, dlname))
            with gzip.open(localfile, 'rt', encoding='ISO-8859-1') as fh:
                LOG.info("Processing %s", localfile)
                self._parse_gene_anatomy(fh, limit)

    def _parse_gene_anatomy(self, fh, limit):
        """
        Process anat_entity files with columns:
        Ensembl gene ID,gene name, anatomical entity ID,
        anatomical entity name, rank score, XRefs to BTO

        :param fh: filehandle
        :param limit: int, limit per group
        :return: None
        """
        dataframe = pd.read_csv(fh, sep='\t')
        col = self.files['anat_entity']['columns']
        if not self.check_fileheader(col, list(dataframe)):
            pass

        gene_groups = dataframe.sort_values(
            'rank score', ascending=False).groupby('Ensembl gene ID')

        if limit is None:
            limit = 20
        gene_groups = gene_groups.head(limit).groupby('Ensembl gene ID')

        for gene, group in gene_groups:
            for index, row in group.iterrows():
                self._add_gene_anatomy_association(
                    row['Ensembl gene ID'].strip(),
                    row['anatomical entity ID'].strip(),
                    row['rank score']
                )
                # uberon <==> bto equivelance?

    def _add_gene_anatomy_association(self, gene_id, anatomy_curie, rank):
        """
        :param gene_id: str Non curified ID
        :param gene_label: str Gene symbol
        :param anatomy_curie: str curified anatomy term
        :param rank: str rank
        :return: None
        """
        g2a_association = Assoc(self.graph, self.name)
        model = Model(self.graph)
        gene_curie = "ENSEMBL:{}".format(gene_id)

        rank = re.sub(r',', '', str(rank))  # ? can't do RE on a float ...
        model.addIndividualToGraph(gene_curie, None)
        g2a_association.sub = gene_curie
        g2a_association.obj = anatomy_curie
        g2a_association.rel = self.globaltt['expressed in']
        g2a_association.add_association_to_graph()
        g2a_association.add_predicate_object(
            self.globaltt['has_quantifier'], float(rank), 'Literal', 'xsd:float')

    # Override
    def checkIfRemoteIsNewer(self, localfile, remote_size, remote_modify):
        """
        Overrides checkIfRemoteIsNewer in Source class

        :param localfile: str file path
        :param remote_size: str bytes
        :param remote_modify: str last modify date in the form 20160705042714
        :return: boolean True if remote file is newer else False
        """
        is_remote_newer = False
        status = os.stat(localfile)
        LOG.info(
            "\nLocal file size: %i"
            "\nLocal Timestamp: %s",
            status[ST_SIZE], datetime.fromtimestamp(status.st_mtime))
        remote_dt = Bgee._convert_ftp_time_to_iso(remote_modify)

        if remote_dt != datetime.fromtimestamp(status.st_mtime) or \
                status[ST_SIZE] != int(remote_size):
            is_remote_newer = True
            LOG.info(
                "Object on server is has different size %i and/or date %s",
                remote_size, remote_dt)

        return is_remote_newer

    @staticmethod
    def _convert_ftp_time_to_iso(ftp_time):
        """
        Convert datetime in the format 20160705042714 to a datetime object

        :return: datetime object
        """
        date_time = datetime(
            int(ftp_time[:4]), int(ftp_time[4:6]), int(ftp_time[6:8]),
            int(ftp_time[8:10]), int(ftp_time[10:12]), int(ftp_time[12:14]))
        return date_time

    def _get_file_list(self, working_dir, file_regex=re.compile(r'.*'), ftp=None):
        """
        Get file list from ftp server filtered by taxon
        :return: Tuple of (Generator object with Tuple(
            file name, info object), ftp object)
        """

        if ftp is None:
            ftp = ftplib.FTP(BGEE_FTP)
            ftp.login("anonymous", "info@monarchinitiative.org")

        working_dir = "{}{}".format(self.version, working_dir)
        LOG.info('Looking for remote files in %s', working_dir)

        ftp.cwd(working_dir)
        remote_files = ftp.nlst()

        # LOG.info('All remote files \n%s', '\n'.join(remote_files))
        files_to_download = [
            dnload for dnload in remote_files if re.match(file_regex, dnload) and
            re.findall(r'^\d+', dnload)[0] in self.tax_ids]
        # LOG.info('Choosing remote files \n%s', '\n'.join(list(files_to_download)))

        return files_to_download, ftp
