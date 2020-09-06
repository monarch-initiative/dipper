import logging
import xml.etree.ElementTree as ET
from dipper.sources.Source import Source
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class Orphanet(Source):
    """
    Orphanet’s aim is to help improve the diagnosis,
    care and treatment of patients with rare diseases.
    For Orphanet, we are currently only parsing the disease-gene associations.
    """

    """
    see: dipper/resources/orphanet/

    for finding what you as most likely working on
    xmlstarlet sel -t -c \
    "JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/DisorderGeneAssociation" \
        en_product6.xml

    """

    files = {
        'disease-gene': {
            'file': 'en_product6.xml',
            'url': 'http://www.orphadata.org/data/xml/en_product6.xml'}
    }

    def __init__(
            self, graph_type, are_bnodes_skolemized, data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='orphanet',
            ingest_title='Orphanet',
            ingest_url='http://www.orpha.net',
            ingest_logo='source-orphanet.png',
            license_url=None,
            data_rights='http://www.orphadata.org/cgi-bin/index.php'
            # file_handle=None
        )

        if 'disease' not in self.all_test_ids:
            LOG.warning("not configured with disease test ids.")
            self.test_ids = []
        else:
            self.test_ids = self.all_test_ids['disease']

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced:
        :return:
        """
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        self._process_diseasegene(limit)

        LOG.info("Done parsing.")

    def _process_diseasegene(self, limit):
        """
        :param limit:
        :return:
        """
        src_key = 'disease-gene'

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        model = Model(graph)

        xmlfile = '/'.join((self.rawdir, self.files[src_key]['file']))

        for event, elem in ET.iterparse(xmlfile):
            if elem.tag != 'Disorder':
                continue
            orphanumber = elem.find('OrphaCode').text
            disorder_curie = 'ORPHA:' + str(orphanumber)

            if self.test_mode and disorder_curie not in self.all_test_ids['disease']:
                continue

            # Orphanet mappings are expected to be in Mondo
            # these free-text disorder names become synonyms
            disorder_label = elem.find('Name').text
            model.addClassToGraph(
                disorder_curie, disorder_label, class_category=blv.terms['Disease'])

            assoc_list = elem.find('DisorderGeneAssociationList')
            expected_genes = assoc_list.get('count')
            if expected_genes == 0:
                LOG.info("%s has no genes.", disorder_curie)
                continue
            # LOG.info(  # too chatty in the logs
            #    'Expecting %s genes associated with disorder %s.',
            #    expected_genes, disorder_curie)
            processed_genes = 0
            for assoc in assoc_list.findall('DisorderGeneAssociation'):
                processed_genes += 1
                gene = assoc.find('Gene')

                # new as of 2020-May
                # (on ORPHA association)  geno? type protein coding gene
                #   7678 gene with protein product
                #   84 Non-coding RNA
                #   31 Disorder-associated locus
                # gene_type = gene.find('GeneType/Name').text

                # todo monochrom mapping?
                # gene_cyto = gene.find('LocusList/Locus/GeneLocus')

                # get gene's curie as a map {'prefix': 'localid'}
                # circa 2020-May
                #   7787 HGNC
                #   7745 Ensembl
                #   7711 OMIM
                #   7709 SwissProt
                #   7578 Genatlas
                #   6523 Reactome
                #   1101 IUPHAR
                gene_clique = {}
                gene_curie = None
                for gene_ref in gene.findall(
                        './ExternalReferenceList/ExternalReference'):
                    prefix = gene_ref.find('Source').text
                    if prefix in self.localtt:
                        prefix = self.localtt[prefix]
                    gene_clique[prefix] = gene_ref.find('Reference').text

                if len(gene_clique) == 0:
                    LOG.error('No gene at all for %s', disorder_curie)
                    break

                # gene representation to prefer
                for prefix in ('HGNC', 'ENSEMBL', 'SwissProt', 'OMIM'):
                    if prefix in gene_clique:
                        gene_curie = prefix + ':' + gene_clique[prefix]
                        gene_clique.pop(prefix)
                        model.addClassToGraph(gene_curie, None)
                        break  # one shot

                if gene_curie is None:
                    # settle for whatever
                    LOG.warning("No prefered gene have\n\t%s", gene_clique)
                    for prefix in gene_clique:
                        gene_curie = prefix + ':' + gene_clique[prefix]
                        gene_clique.pop(prefix)
                        model.addClassToGraph(gene_curie, None)
                        break  # one shot

                for prefix in gene_clique:
                    lclid = gene_clique[prefix]
                    dbxref = ':'.join((prefix, lclid))
                    if gene_curie != dbxref:
                        model.addClassToGraph(
                            dbxref, None, class_category=blv.terms['Gene']
                        )
                        model.addEquivalentClass(gene_curie, dbxref)
                        # note: when dbxref is a Genatlas family,
                        # the curie is/was mistaken as a literal
                        # due to the trailing '@' on the symbol

                syn_list = gene.find('./SynonymList')
                if int(syn_list.get('count')) > 0 and gene_curie is not None:
                    for syn in syn_list.findall('./Synonym'):
                        if syn is not None and syn.text != '':
                            model.addSynonym(gene_curie, syn.text)

                dg_label = assoc.find('./DisorderGeneAssociationType/Name').text
                # circa 2020-May
                #   4771 Disease-causing germline mutation(s) in
                #   1137 Disease-causing germline mutation(s) (loss of function) in
                #   576 Major susceptibility factor in
                #   359 Candidate gene tested in
                #   232 Role in the phenotype of
                #   232 Part of a fusion gene in
                #   211 Disease-causing germline mutation(s) (gain of function) in
                #   186 Disease-causing somatic mutation(s) in
                #   45 Biomarker tested in
                #   44 Modifying germline mutation in
                rel_curie = self.resolve(dg_label)

                # genotype_curie = self.resolve(" | ".join((dg_label, gene_type)))

                # use dg association status to issue an evidence code
                # FIXME these codes may be sub-optimal (there are only two)
                # maybe just attach a "pending" to the minority that need it.
                eco_id = self.resolve(
                    assoc.find('DisorderGeneAssociationStatus/Name').text)

                g2p_assoc = G2PAssoc(
                    self.graph,      # graph
                    self.name,       # definedby
                    gene_curie,      # entity_id
                    disorder_curie,  # phenotype_id
                    rel_curie        # rel=None
                )

                g2p_assoc.add_evidence(eco_id)
                g2p_assoc.add_association_to_graph()

            elem.clear()  # empty the element
            if int(expected_genes) != processed_genes:
                LOG.warning(
                    '%s expected %i associated genes but we processed %i',
                    disorder_curie, int(expected_genes), processed_genes)
