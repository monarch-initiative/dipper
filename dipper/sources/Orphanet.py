import logging
import xml.etree.ElementTree as ET

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model
from dipper import config

logger = logging.getLogger(__name__)


class Orphanet(Source):
    """
    Orphanet’s aim is to help improve the diagnosis,
    care and treatment of patients with rare diseases.
    For Orphanet, we are currently only parsing the disease-gene associations.

     Note that ???

    """

    files = {
        'disease-gene': {
            'file': 'en_product6.xml',
            'url': 'http://www.orphadata.org/data/xml/en_product6.xml'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'orphanet')

        self.dataset = Dataset(
            'orphanet', 'Orphanet', 'http://www.orpha.net', None,
            'http://creativecommons.org/licenses/by-nd/3.0/',
            'http://omim.org/help/agreement')

        # check to see if there's any ids configured in the config;
        # otherwise, warn
        if 'test_ids' not in config.get_config() or \
                'disease' not in config.get_config()['test_ids']:
            logger.warning("not configured with disease test ids.")

        return

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced:
        :return:
        """
        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_diseasegene(limit)

        logger.info("Done parsing.")

        return

    def _process_diseasegene(self, limit):
        """
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        model = Model(g)

        myfile = '/'.join((self.rawdir, self.files['disease-gene']['file']))

        # PYLINT complains iterparse deprecated,
        # but as of py 3.4 only the optional & unsupplied parse arg is.
        for event, elem in ET.iterparse(myfile):
            if elem.tag == 'Disorder':
                # get the element name and id, ignoreS element name
                # id = elem.get('id') # some internal identifier
                disorder_num = elem.find('OrphaNumber').text

                disorder_id = 'Orphanet:'+str(disorder_num)

                if self.testMode and \
                        disorder_id not in \
                        config.get_config()['test_ids']['disease']:
                    continue

                disorder_label = elem.find('Name').text

                # make a hash of internal gene id to type for later lookup
                gene_iid_to_type = {}
                gene_list = elem.find('GeneList')
                for gene in gene_list.findall('Gene'):
                    gene_iid = gene.get('id')
                    gene_type = gene.find('GeneType').get('id')
                    gene_iid_to_type[gene_iid] = gene_type

                # assuming that these are in the ontology
                model.addClassToGraph(disorder_id, disorder_label)

                assoc_list = elem.find('DisorderGeneAssociationList')
                for a in assoc_list.findall('DisorderGeneAssociation'):
                    gene_iid = a.find('.//Gene').get('id')
                    gene_name = a.find('.//Gene/Name').text
                    gene_symbol = a.find('.//Gene/Symbol').text
                    gene_num = a.find('./Gene/OrphaNumber').text
                    gene_id = 'Orphanet:'+str(gene_num)
                    gene_type_id = \
                        self._map_gene_type_id(gene_iid_to_type[gene_iid])
                    model.addClassToGraph(
                        gene_id, gene_symbol, gene_type_id, gene_name)
                    syn_list = a.find('./Gene/SynonymList')
                    if int(syn_list.get('count')) > 0:
                        for s in syn_list.findall('./Synonym'):
                            model.addSynonym(gene_id, s.text)

                    dgtype = a.find('DisorderGeneAssociationType').get('id')
                    rel_id = self._map_rel_id(dgtype)
                    dg_label = \
                        a.find('./DisorderGeneAssociationType/Name').text
                    if rel_id is None:
                        logger.warning(
                            "Cannot map association type (%s) to RO " +
                            "for association (%s | %s).  Skipping.",
                            dg_label, disorder_label, gene_symbol)
                        continue

                    alt_locus_id = '_:'+gene_num+'-'+disorder_num+'VL'
                    alt_label = \
                        ' '.join(('some variant of', gene_symbol.strip(),
                                  'that is a', dg_label.lower(),
                                  disorder_label))

                    model.addIndividualToGraph(alt_locus_id, alt_label,
                                               geno.genoparts['variant_locus'])
                    geno.addAffectedLocus(alt_locus_id, gene_id)
                    model.addBlankNodeAnnotation(alt_locus_id)

                    # consider typing the gain/loss-of-function variants like:
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002054
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002053

                    # use "assessed" status to issue an evidence code
                    # FIXME I think that these codes are sub-optimal
                    status_code = \
                        a.find('DisorderGeneAssociationStatus').get('id')
                    # imported automatically asserted information
                    # used in automatic assertion
                    eco_id = 'ECO:0000323'
                    # Assessed
                    # TODO are these internal ids stable between releases?
                    if status_code == '17991':
                        # imported manually asserted information
                        # used in automatic assertion
                        eco_id = 'ECO:0000322'
                    # Non-traceable author statement ECO_0000034
                    # imported information in automatic assertion ECO_0000313

                    assoc = G2PAssoc(g, self.name, alt_locus_id,
                                     disorder_id, rel_id)
                    assoc.add_evidence(eco_id)
                    assoc.add_association_to_graph()

                    rlist = a.find('./Gene/ExternalReferenceList')
                    eqid = None

                    for r in rlist.findall('ExternalReference'):
                        if r.find('Source').text == 'Ensembl':
                            eqid = 'ENSEMBL:'+r.find('Reference').text
                        elif r.find('Source').text == 'HGNC':
                            eqid = 'HGNC:'+r.find('Reference').text
                        elif r.find('Source').text == 'OMIM':
                            eqid = 'OMIM:'+r.find('Reference').text
                        else:
                            pass  # skip the others for now
                        if eqid is not None:
                            model.addClassToGraph(eqid, None)
                            model.addEquivalentClass(gene_id, eqid)
                elem.clear()  # empty the element

            if self.testMode and limit is not None and line_counter > limit:
                return

        return

    # TODO rebuilding the map every time the function is called is
    # not at all efficent ...(try decorator staticmethod or remove from class?)
    # @staticmethod
    def _map_rel_id(self, orphanet_rel_id):
        # To check on  the "DisorderGeneAssociationType" to map
        # xmlstarlet sel -t -m \
        # './JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/\
        # DisorderGeneAssociation/DisorderGeneAssociationType'\
        # -v './@id' -o '    ' -v './Name' -n en_product6.xml |\
        # sort | uniq -c | sort -nr
        rel_id = None
        model = Model(self.graph)
        id_map = {
            # Disease-causing germline mutation(s) in (RO:0002200)
            '17949': model.object_properties['has_phenotype'],
            # Disease-causing somatic mutation(s) in
            '17955': model.object_properties['has_phenotype'],
            # Major susceptibility factor in  ('RO:0002326')
            '17961': model.object_properties['contributes_to'],
            # Modifying germline mutation in
            '17967': model.object_properties['contributes_to'],
            # Modifying somatic mutation in
            '17973': model.object_properties['contributes_to'],
            # Part of a fusion gene in
            '17979': model.object_properties['contributes_to'],
            # Role in the phenotype of
            '17985': model.object_properties['contributes_to'],
            '18273': None,  # Candidate gene tested in  FIXME?
            # Disease-causing germline mutation(s) (gain of function) in
            '25979': model.object_properties['has_phenotype'],
            # Disease-causing germline mutation(s) (loss of function) in
            '25972': model.object_properties['has_phenotype'],
        }

        if orphanet_rel_id in id_map:
            rel_id = id_map[orphanet_rel_id]
        else:
            logger.error(
                'Disease-gene association type (%s) not mapped.',
                orphanet_rel_id)

        return rel_id

    # To check the "GeneType" mappings
    # xmlstarlet sel -t -v \
    #   './JDBOR/DisorderList/Disorder//Gene/GeneType/Name'
    #   en_product6.xml | sort | uniq -c | sort -nr
    @staticmethod
    def _map_gene_type_id(orphanet_type_id):
        type_id = Genotype.genoparts['sequence_feature']
        id_map = {
            # locus
            '25986': Genotype.genoparts['sequence_feature'],
            # gene with protein product
            '25993': Genotype.genoparts['protein_coding_gene'],
            # Non-coding RNA
            '26046': Genotype.genoparts['ncRNA_gene']
        }

        if orphanet_type_id in id_map:
            type_id = id_map[orphanet_type_id]
        else:
            logger.error(
                'Gene type (%s) not mapped. Defaulting to SO:sequence_feature',
                orphanet_type_id)

        return type_id

    def getTestSuite(self):
        import unittest
        # TODO PYLINT Unable to import 'tests.test_orphanet
        from tests.test_orphanet import OrphanetTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(OrphanetTestCase)

        return test_suite
