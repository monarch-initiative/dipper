import logging
import os
from stat import ST_SIZE
import re
from zipfile import ZipFile
from datetime import datetime
import time
from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
from dipper.models.assoc.Association import Assoc
from dipper.models.Genotype import Genotype
import ftplib
import pandas as pd

logger = logging.getLogger(__name__)


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

    BGEE_FTP = 'ftp.bgee.org'
    DEFAULT_TAXA = [10090, 10116, 13616, 28377, 6239,
                    7227, 7955, 8364, 9031, 9258,
                    9544, 9593, 9597, 9598, 9606,
                    9823, 9913]
    files = {
        'anat_entity': {
            'path': '/download/ranks/anat_entity/',
            'pattern': re.compile(r'.*_all_data_.*')
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None, version=None):
        """
        :param tax_ids: [int,], List of taxa
        :return:
        """
        super().__init__(graph_type, are_bnodes_skolemized, 'bgee')
        if tax_ids is None:
            self.tax_ids = Bgee.DEFAULT_TAXA
        else:
            logger.info("Filtering on taxa {}".format(tax_ids))
            self.tax_ids = tax_ids

        self.dataset = Dataset(
            'Bgee', 'Bgee Gene expression data in animals',
            'http://bgee.org/')

        if version is None:
            self.version = 'current'
        else:
            self.version = version

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced: boolean, force download
        :return:
        """
        for group in self.files:
            files_to_download, ftp = \
                self._get_file_list(self.files[group]['path'],
                                    self.files[group]['pattern'])
            for name, info in files_to_download:
                localfile = '/'.join((self.rawdir, name))
                if not os.path.exists(localfile)\
                        or is_dl_forced\
                        or self.checkIfRemoteIsNewer(localfile, info['size'],
                                                     info['modify']):
                    logger.info("Fetching {}".format(name))
                    logger.info("Writing to {}".format(localfile))
                    ftp.retrbinary('RETR {}'.format(name), open(localfile, 'wb').write)
                    remote_dt = Bgee._convert_ftp_time_to_iso(info['modify'])
                    os.utime(localfile, (time.mktime(remote_dt.timetuple()),
                                         time.mktime(remote_dt.timetuple())))

        ftp.quit()

        return

    def parse(self, limit=None):
        """
        Given the input taxa, expects files in the raw directory
        with the name {tax_id}_anat_entity_all_data_Pan_troglodytes.tsv.zip
        :param limit: int Limit to top ranked anatomy associations per group
        :return: None
        """
        files_to_download, ftp = \
            self._get_file_list(self.files['anat_entity']['path'],
                                self.files['anat_entity']['pattern'])
        for name, info in files_to_download:
            localfile = '/'.join((self.rawdir, name))
            with ZipFile(localfile, 'r') as zip_file:
                fh = zip_file.open(re.sub(r'\.zip$', '', name))
                logger.info("Processing {}".format(name))
                self._parse_gene_anatomy(fh, limit)

        return

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
        gene_groups = dataframe.sort_values('rank score', ascending=False)\
                               .groupby('Ensembl gene ID')

        if limit is not None:
            gene_groups = gene_groups.head(limit).groupby('Ensembl gene ID')

        for gene, group in gene_groups:
            for index, row in group.iterrows():
                self._add_gene_anatomy_association(
                    row['Ensembl gene ID'], row['anatomical entity ID'],
                    row['rank score']
                )
        return

    def _add_gene_anatomy_association(self, gene_id, anatomy_curie, rank):
        """
        :param gene_id: str Non curified ID
        :param gene_label: str Gene symbol
        :param anatomy_curie: str curified anatomy term
        :param rank: str rank
        :return: None
        """
        g2a_association = Assoc(self.graph, self.name)
        genotype = Genotype(self.graph)
        model = Model(self.graph)
        gene_curie = "ENSEMBL:{}".format(gene_id)
        rank = re.sub(r',', '', rank)
        model.addIndividualToGraph(ind_id=gene_curie, label=None,
                                   ind_type=genotype.genoparts['gene'])
        g2a_association.sub = gene_curie
        g2a_association.obj = anatomy_curie
        g2a_association.rel = Assoc.object_properties['expressed_in']
        g2a_association.add_association_to_graph()
        g2a_association.add_predicate_object(
            Assoc.datatype_properties['has_quantifier'],
            float(rank), 'Literal', 'xsd:float')
        return

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
        logger.info("Local file date: {0}, size: {1}".format(
                    datetime.fromtimestamp(status.st_mtime),
                    status[ST_SIZE]))
        remote_dt = Bgee._convert_ftp_time_to_iso(remote_modify)

        if remote_dt != datetime.fromtimestamp(status.st_mtime) \
                or status[ST_SIZE] != int(remote_size):
            is_remote_newer = True
            logger.info(
                "Object on server is has different size {0} and/or "
                "date {1}".format(remote_size, remote_dt))

        return is_remote_newer

    @staticmethod
    def _convert_ftp_time_to_iso(ftp_time):
        """
        Convert datetime in the format 20160705042714 to a datetime object
        :return: datetime object
        """
        date_time = datetime(int(ftp_time[:4]), int(ftp_time[4:6]),
                             int(ftp_time[6:8]), int(ftp_time[8:10]),
                             int(ftp_time[10:12]), int(ftp_time[12:14]))
        return date_time

    def _get_file_list(self, working_dir, file_regex=re.compile(r'.*'), ftp=None):
        """
        Get file list from ftp server filtered by taxon
        :return: Tuple of (Generator object with Tuple(file name, info object), ftp object)
        """
        if ftp is None:
            ftp = ftplib.FTP(Bgee.BGEE_FTP)
            ftp.login("anonymous", "info@monarchinitiative.edu")

        working_dir = "{}{}".format(self.version, working_dir)

        ftp.cwd(working_dir)

        directory = ftp.mlsd()
        files = (value for value in directory if value[1]['type'] == 'file')
        files_to_download = (value for value in files
                             if re.match(file_regex, value[0])
                             and int(re.findall(r'^\d+', value[0])[0]) in self.tax_ids)

        return files_to_download, ftp
