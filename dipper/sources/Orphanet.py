import logging
import xml.etree.ElementTree as ET

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.utils.GraphUtils import GraphUtils
from dipper import config
from dipper import curie_map


logger = logging.getLogger(__name__)


class Orphanet(Source):
    """
    Orphanet’s aim is to help improve the diagnosis, care and treatment of patients with rare diseases.
    For Orphanet, we are currently only parsing the disease-gene associations.


     Note that
    """

    files = {
        'disease-gene': {
            'file': 'en_product6.xml',
            'url': 'http://www.orphadata.org/data/xml/en_product6.xml'}
    }

    def __init__(self):
        Source.__init__(self, 'orphanet')

        self.load_bindings()

        self.dataset = Dataset('orphanet', 'Orphanet', 'http://www.orpha.net',
                               None,
                               'http://creativecommons.org/licenses/by-nd/3.0/',
                               'http://omim.org/help/agreement')

        # check to see if there's any ids configured in the config; otherwise, warn
        if 'test_ids' not in config.get_config() or 'disease' not in config.get_config()['test_ids']:
            logger.warn("not configured with disease test ids.")

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

        self.load_core_bindings()
        self.load_bindings()

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
        gu = GraphUtils(curie_map.get())

        myfile = '/'.join((self.rawdir, self.files['disease-gene']['file']))

        for event, elem in ET.iterparse(myfile):
            if elem.tag == 'Disorder':
                # get the element name and id
                # id = elem.get('id') # some internal identifier
                disorder_num = elem.find('OrphaNumber').text

                disorder_id = 'Orphanet:'+str(disorder_num)

                if self.testMode and disorder_id not in config.get_config()['test_ids']['disease']:
                    continue

                disorder_label = elem.find('Name').text

                # make a hash of internal gene id to type for later lookup
                gene_iid_to_type = {}
                gene_list = elem.find('GeneList')
                for gene in gene_list.findall('Gene'):
                    gene_iid = gene.get('id')
                    gene_type = gene.find('GeneType').get('id')
                    gene_iid_to_type[gene_iid] = gene_type

                gu.addClassToGraph(g, disorder_id, disorder_label)   # assuming that these are in the ontology

                assoc_list = elem.find('DisorderGeneAssociationList')
                for a in assoc_list.findall('DisorderGeneAssociation'):
                    gene_iid = a.find('.//Gene').get('id')
                    gene_name = a.find('.//Gene/Name').text
                    gene_symbol = a.find('.//Gene/Symbol').text
                    gene_num = a.find('./Gene/OrphaNumber').text
                    gene_id = 'Orphanet:'+str(gene_num)
                    gene_type_id = self._map_gene_type_id(gene_iid_to_type[gene_iid])
                    gu.addClassToGraph(g, gene_id, gene_symbol, gene_type_id, gene_name)
                    syn_list = a.find('./Gene/SynonymList')
                    if int(syn_list.get('count')) > 0:
                        for s in syn_list.findall('./Synonym'):
                            gu.addSynonym(g, gene_id, s.text)

                    dgtype = a.find('DisorderGeneAssociationType').get('id')
                    rel_id = self._map_rel_id(dgtype)
                    dg_label = a.find('./DisorderGeneAssociationType/Name').text
                    if rel_id is None:
                        logger.warn("Cannot map association type (%s) to RO for association (%s | %s).  Skipping.",
                                    dg_label, disorder_label, gene_symbol)
                        continue

                    alt_locus_id = '_'+gene_num+'-'+disorder_num+'VL'
                    alt_label = ' '.join(('some variant of', gene_symbol.strip(),
                                          'that is a', dg_label.lower(), disorder_label))
                    if self.nobnodes:
                        alt_locus_id = ':'+alt_locus_id
                    gu.addIndividualToGraph(g, alt_locus_id, alt_label, geno.genoparts['variant_locus'])
                    geno.addAlleleOfGene(alt_locus_id, gene_id)

                    # consider typing the gain/loss-of-function variants like:
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002054
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002053

                    # use "assessed" status to issue an evidence code
                    # FIXME I think that these codes are sub-optimal
                    status_code = a.find('DisorderGeneAssociationStatus').get('id')
                    eco_id = 'ECO:0000323'  # imported automatically asserted information used in automatic assertion
                    if status_code == '17991':  # Assessed  # TODO are these internal ids stable between releases?
                        eco_id = 'ECO:0000322'  # imported manually asserted information used in automatic assertion
                    # Non-traceable author statement ECO_0000034
                    # imported information in automatic assertion ECO_0000313

                    assoc = G2PAssoc(self.name, alt_locus_id, disorder_id, rel_id)
                    assoc.add_evidence(eco_id)
                    assoc.add_association_to_graph(g)

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
                            gu.addClassToGraph(g, eqid, None)
                            gu.addEquivalentClass(g, gene_id, eqid)
                            pass
                elem.clear()  # discard the element

            if self.testMode and limit is not None and line_counter > limit:
                return

        gu.loadProperties(g, G2PAssoc.annotation_properties, G2PAssoc.ANNOTPROP)
        gu.loadProperties(g, G2PAssoc.datatype_properties, G2PAssoc.DATAPROP)
        gu.loadProperties(g, G2PAssoc.object_properties, G2PAssoc.OBJECTPROP)
        gu.loadAllProperties(g)

        return

    @staticmethod
    def _map_rel_id(orphanet_rel_id):
        # TODO check if these ids are stable for mapping
        rel_id = None
        gu = GraphUtils(curie_map.get())
        id_map = {
            '17949': gu.object_properties['has_phenotype'],  # Disease-causing germline mutation(s) in
            '17955': gu.object_properties['has_phenotype'],  # Disease-causing somatic mutation(s) in
            '17961': gu.object_properties['contributes_to'],  # Major susceptibility factor in
            '17967': gu.object_properties['contributes_to'],  # Modifying germline mutation in
            '17973': gu.object_properties['contributes_to'],  # Modifying somatic mutation in
            '17979': gu.object_properties['contributes_to'],  # Part of a fusion gene in
            '17985': gu.object_properties['contributes_to'],  # Role in the phenotype of
            '18273': None,  # Candidate gene tested in  FIXME?
            '25972': gu.object_properties['has_phenotype'],  # Disease-causing germline mutation(s) (loss of function) in
            '25979': gu.object_properties['has_phenotype']   # Disease-causing germline mutation(s) (gain of function) in
        }

        if orphanet_rel_id in id_map:
            rel_id = id_map[orphanet_rel_id]
        else:
            logger.error('Disease-gene association type (%s) not mapped.', orphanet_rel_id)

        return rel_id

    @staticmethod
    def _map_gene_type_id(orphanet_type_id):
        # TODO check if these ids are stable for mapping

        type_id = Genotype.genoparts['sequence_feature']
        id_map = {
            '25986': Genotype.genoparts['sequence_feature'],  # locus
            '25993': Genotype.genoparts['protein_coding_gene'],  # gene with protein product
            '26046': Genotype.genoparts['ncRNA_gene']  # Non-coding RNA
        }

        if orphanet_type_id in id_map:
            type_id = id_map[orphanet_type_id]
        else:
            logger.error('Gene type (%s) not mapped.  Defaulting to SO:sequence_feature', orphanet_type_id)

        return type_id

    def getTestSuite(self):
        import unittest
        from tests.test_orphanet import OrphanetTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(OrphanetTestCase)

        return test_suite
