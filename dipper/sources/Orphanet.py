import logging
import xml.etree.ElementTree as ET
from dipper.sources.Source import Source
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model

LOG = logging.getLogger(__name__)


class Orphanet(Source):
    """
    Orphanet’s aim is to help improve the diagnosis,
    care and treatment of patients with rare diseases.
    For Orphanet, we are currently only parsing the disease-gene associations.
    """

    """ 
    Some useful code:
    xmlstarlet sel -t  -v "/JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/
    """
    
    files = {
        'disease-gene': {
            'file': 'en_product6.xml',
            'url': 'http://www.orphadata.org/data/xml/en_product6.xml'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'orphanet',
            ingest_title='Orphanet',
            ingest_url='http://www.orpha.net',
            license_url=None,
            data_rights='http://www.orphadata.org/cgi-bin/index.php'
            # file_handle=None
        )

        if 'disease' not in self.all_test_ids:
            LOG.warning("not configured with disease test ids.")
            self.test_ids = []
        else:
            self.test_ids = self.all_test_ids['disease']

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
            LOG.info("Only parsing first %d rows", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self._process_diseasegene(limit)

        LOG.info("Done parsing.")

        return

    def _process_diseasegene(self, limit):
        """
        :param limit:
        :return:
        """
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        line_counter = 0

        model = Model(graph)

        myfile = '/'.join((self.rawdir, self.files['disease-gene']['file']))

        for event, elem in ET.iterparse(myfile):
            if elem.tag == 'Disorder':
                # get the element name and id, ignore element name
                # id = elem.get('id') # some internal identifier
                disorder_num = elem.find('OrphaNumber').text
                disorder_id = 'ORPHA:' + str(disorder_num)

                if self.test_mode and disorder_id not in self.all_test_ids['disease']:
                    continue
                disorder_label = elem.find('Name').text

                # assuming that these are in the ontology (...any particular one?)
                model.addClassToGraph(disorder_id, disorder_label)
                assoc_list = elem.find('DisorderGeneAssociationList')
                expected_genes = assoc_list.get('count')
                LOG.info(
                    'Expecting %s genes associated with disorder %s.',
                    expected_genes, disorder_id)
                processed_genes = 0
                for assoc in assoc_list.findall('DisorderGeneAssociation'):
                    processed_genes += 1
                    gene = assoc.find('Gene')

                    # get gene's curie  HGNC or Ensembl ...

                    lclid = gene.find('OrphaNumber').text
                    gene_curie = 'ORPHA:' + lclid
                    gene_set = {'ORPHA': lclid}
                    for gene_ref in gene.findall(
                            './ExternalReferenceList/ExternalReference'):
                        gene_set[gene_ref.find('Source').text] = \
                            gene_ref.find('Reference').text

                    # set priority (clique leader if available) but default to OPRHA
                    for pfx in ('HGNC', 'Ensembl', 'SwissProt'):
                        if pfx in gene_set:
                            if pfx in self.localtt:
                                pfx = self.localtt[pfx]
                            gene_curie = pfx + ':' + gene_set[pfx]
                            gene_set.pop(pfx)
                            model.addClassToGraph(gene_curie, None)
                            break

                    # TEC have reservations w.r.t aggerator links being gene classes
                    for prefix in gene_set:
                        lclid = gene_set[prefix]
                        if prefix in self.localtt:
                            prefix = self.localtt[prefix]

                        dbxref = prefix + ':' + lclid

                        if gene_curie != dbxref:
                            model.addClassToGraph(dbxref, None)
                            model.addEquivalentClass(gene_curie, dbxref)

                    # TEC. would prefer this not happen here. let HGNC handle it
                    # except there are some w/o explicit external links ...

                    gene_symbol = gene.find('Symbol').text

                    syn_list = gene.find('./SynonymList')
                    if int(syn_list.get('count')) > 0:
                        for syn in syn_list.findall('./Synonym'):
                            model.addSynonym(gene_curie, syn.text)

                    dg_label = assoc.find('./DisorderGeneAssociationType/Name').text

                    # use dg association status to issue an evidence code
                    # FIXME I think that these codes are sub-optimal
                    eco_id = self.resolve(
                        assoc.find('DisorderGeneAssociationStatus/Name').text)

                    rel_id = self.resolve(dg_label)
                    
                    g2p_assoc = G2PAssoc(self.graph, self.name, gene_curie, disorder_id, rel_id)
                    g2p_assoc.add_evidence(eco_id)
                    g2p_assoc.add_association_to_graph()

                elem.clear()  # empty the element
                if int(expected_genes) != processed_genes:
                    LOG.warning(
                        '% expected %s associated genes but we processed %i',
                        disorder_id, expected_genes, processed_genes)

            if self.test_mode and limit is not None and line_counter > limit:
                return

        return

