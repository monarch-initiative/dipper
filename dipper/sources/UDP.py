from dipper.sources.Source import Source
from dipper import curie_map
from dipper.models.Genotype import Genotype
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Dataset import Dataset
import logging
import csv

logger = logging.getLogger(__name__)


class UDP(Source):
    """
    The National Institutes of Health (NIH) Undiagnosed Diseases Program (UDP)
    is part of the Undiagnosed Disease Network (UDN),
    an NIH Common Fund initiative that focuses on the most puzzling medical cases
    referred to the NIH Clinical Center in Bethesda, Maryland.
    from https://www.genome.gov/27544402/the-undiagnosed-diseases-program/

    Data is available by request for access via the NHGRI collaboration server:
    https://udplims-collab.nhgri.nih.gov/api

    Note this source class does not include a fetch method since the data is private
    The parser works generically when two tsv files are present in the raw directory
    /raw/udp with the structure

    udp_variants.tsv
    'Patient', 'Family', 'Chr', 'Build', 'Chromosome Position',
    'Reference Allele', 'Variant Allele', 'Parent of origin',
    'Allele Type', 'Mutation Type', 'Gene', 'Transcript', 'Original Amino Acid',
    'Variant Amino Acid', 'Amino Acid Change', 'Segregates with',
    'Position', 'Exon', 'Inheritance model', 'Zygosity', 'dbSNP ID', '1K Frequency',
    'Number of Alleles'

    udp_phenotypes.tsv
    'Patient', 'HPID', 'Present'

    The script also utilizes two mapping files
    udp_gene_map.tsv -  generated from scripts/fetch-gene-ids.py,
                        gene symbols from udp_variants
    udp_chr_rs.tsv - rsid(s) per coordinate greped from hg19 dbsnp file,
                     then disambiguated with eutils, see scripts/dbsnp/dbsnp.py

    """
    files = {
        'patient_phenotypes': {
            'file': 'udp_phenotypes.tsv'
        },
        'patient_variants': {
            'file': 'udp_variants.tsv'
        }
    }
    map_files = {
        'gene_map': '../../resources/udp/udp_gene_map',
        'dbsnp_map': '../../resources/udp/dbsnp'
    }

    def __init__(self):
        super().__init__('udp')
        self.dataset = Dataset(
            'udp', 'UDP', 'https://rarediseases.info.nih.gov/')

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

        phenotype_file = '/'.join((self.rawdir, self.files['patient_phenotypes']['file']))
        variant_file = '/'.join((self.rawdir, self.files['patient_variants']['file']))

        self._parse_patient_phenotypes(phenotype_file, limit)
        self._parse_patient_variants(phenotype_file, limit)

        return

    def _parse_patient_variants(self, file, limit):
        """
        :param file: file path
        :param limit: limit (int, optional) limit the number of rows processed
        :return:
        """
        return

    def _parse_patient_phenotypes(self, file, limit):
        """
        :param file: file path
        :param limit: limit (int, optional) limit the number of rows processed
        :return:
        """
        genotype_util = Genotype(self.graph)
        graph_util = GraphUtils(curie_map.get())
        line_counter = 0
        with open(file, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                (patient_id, hpo_curie, present) = row
                patient_curie = ':{0}'.format(patient_id)
                graph_util.addPerson(self.graph, patient_curie, patient_id)

                # Add Genotype
                intrinsic_geno_bnode = self.make_id(
                    "{0}-intrinsic-genotype".format(patient_id), "_")
                genotype_label = "{0} genotype".format(patient_id)
                genotype_util.addGenotype(intrinsic_geno_bnode,
                                          genotype_label,
                                          genotype_util.genoparts['intrinsic_genotype'])

                graph_util.addTriple(self.graph, patient_curie,
                                     genotype_util.object_properties['has_genotype'],
                                     intrinsic_geno_bnode)

                disease_bnode = self.make_id("{0}-disease".format(patient_id), "_")
                disease_label = "{0} disease".format(patient_id)
                graph_util.addIndividualToGraph(self.graph, disease_bnode, disease_label, "DOID:4")
                graph_util.addTriple(self.graph, patient_curie,
                                     graph_util.object_properties['has_phenotype'],
                                     disease_bnode)
                if present == 'yes':
                    graph_util.addTriple(self.graph, patient_curie,
                                         graph_util.object_properties['has_phenotype'],
                                         hpo_curie)

                line_counter += 1
                if not self.testMode and limit is not None \
                        and line_counter >= limit:
                    break
