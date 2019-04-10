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
                        #       'OMIM', 'Genatlas','Reactome', 'IUPHAR'):
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

                    # gene_name = gene.find('Name').text
                    gene_symbol = gene.find('Symbol').text
                    # gene_iid = assoc.find('DisorderGeneAssociationType').get('id')
                    # gene_type_id = self.resolve(gene_iid)
                    # don't  know the 'type' of the gene for this class anymore
                    # model.addClassToGraph(
                    #    gene_curie, gene_symbol, gene_type_id, gene_name)

                    syn_list = gene.find('./SynonymList')
                    if int(syn_list.get('count')) > 0:
                        for syn in syn_list.findall('./Synonym'):
                            model.addSynonym(gene_curie, syn.text)

                    dg_label = assoc.find('./DisorderGeneAssociationType/Name').text
                    # rel_id = self.resolve(dg_label)

                    # alt_locus_id = '_:' + gene_num + '-' + disorder_num + 'VL'
                    # alt_label = ' '.join((
                    #    'some variant of', gene_symbol.strip(), disorder_label))
                    # model.addIndividualToGraph(
                    #    alt_locus_id, alt_label, self.globaltt['variant_locus'])
                    # geno.addAffectedLocus(alt_locus_id, gene_id)
                    # model.addBlankNodeAnnotation(alt_locus_id)
                    # consider typing the gain/loss-of-function variants like:
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002054
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002053

                    # use dg association status to issue an evidence code
                    # FIXME I think that these codes are sub-optimal
                    eco_id = self.resolve(
                        assoc.find('DisorderGeneAssociationStatus/Name').text)

                    # assoc = G2PAssoc(
                    #    graph, self.name, alt_locus_id, disorder_id, rel_id)
                    # assoc.add_evidence(eco_id)
                    # assoc.add_association_to_graph()

                    self.add_gene_to_disease(
                        dg_label, gene_curie, gene_symbol, disorder_id, eco_id)

                elem.clear()  # empty the element
                if int(expected_genes) != processed_genes:
                    LOG.warning(
                        '% expected %s associated genes but we processed %i',
                        disorder_id, expected_genes, processed_genes)

            if self.test_mode and limit is not None and line_counter > limit:
                return

        return

    def add_gene_to_disease(
            self,
            association_type,
            gene_id,
            gene_symbol,
            disease_id,
            eco_id):
        """
        Composes triples based on the DisorderGeneAssociationType element:
        AND the suffixes:

            - "gene phenotype"
            - "function consequence"
            - "cell origin"

        xmlstarlet sel -t  -v "/JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/
            DisorderGeneAssociation/DisorderGeneAssociationType/Name" en_product6.xml  \
            | sort -u

        Biomarker tested in
        Candidate gene tested in
        Disease-causing germline mutation(s) (gain of function) in
        Disease-causing germline mutation(s) in
        Disease-causing germline mutation(s) (loss of function) in
        Disease-causing somatic mutation(s) in
        Major susceptibility factor in
        Modifying germline mutation in
        Part of a fusion gene in
        Role in the phenotype of

        These labels are a composition of terms, we map:
        gene-disease predicate (has phenotype, contributes_to)
        variant-origin (germline, somatic)
        variant-functional consequence (loss, gain)

        To check on  the "DisorderGeneAssociationType" to id-label map
        xmlstarlet sel -t -m \
        './JDBOR/DisorderList/Disorder/DisorderGeneAssociationList/\
        DisorderGeneAssociation/DisorderGeneAssociationType'\
        -v './@id' -o '    ' -v './Name' -n en_product6.xml |\
        sort | uniq -c | sort -nr

        Although the id-label pairs appear to be stable after
        a few years, we map to the label instead of the id in
        case Orphanet changes their IDs

        :param association_type: {str} DisorderGeneAssociationType/Name,
                                       eg Role in the phenotype of

        :param gene_id: {str} gene id as curie
        :param gene_symbol: {str} HGVS gene symbol
        :param disease_id: {str} disease id as curie
        :param eco_id: {str} eco code as curie

        :return: None
        """

        model = Model(self.graph)
        geno = Genotype(self.graph)
        gene_or_variant = ""

        # If we know something about the variant such as functional consequence or
        # cellular origin make a blank node and attach the attributes
        is_variant = False
        variant_id_string = "{}{}".format(gene_id, disease_id)
        functional_consequence = None
        cell_origin = None

        # hard fail for no mappings/new terms, otherwise they go unnoticed
        if "{}|gene phenotype".format(association_type) not in self.localtt:
            raise ValueError(
                'Disease-gene association type {} not mapped'.format(association_type)
            )

        g2p_relation = self.resolve("|".join([association_type, "gene phenotype"]))

        # Variant attributes
        if "|".join([association_type, "function consequence"]) in self.localtt:
            is_variant = True
            local_key = "|".join([association_type, "function consequence"])
            functional_consequence = self.resolve(local_key)
            functional_consequence_lbl = self.localtt[local_key]
        if "|".join([association_type, "cell origin"]) in self.localtt:
            is_variant = True
            local_key = "|".join([association_type, "cell origin"])
            cell_origin = self.resolve(local_key)
            cell_origin_lbl = self.localtt[local_key]

        if is_variant:
            variant_label = "of {}".format(gene_symbol)
            if functional_consequence:
                variant_label = "{} {}".format(
                    functional_consequence_lbl.replace('_', ' '), variant_label
                )
                variant_id_string += functional_consequence_lbl
            else:
                variant_label = "variant {}".format(variant_label)

            if cell_origin:
                variant_label = "{} {}".format(cell_origin_lbl, variant_label)
                variant_id_string += cell_origin_lbl

            variant_bnode = self.make_id(variant_id_string, "_")
            model.addIndividualToGraph(
                variant_bnode, variant_label, self.globaltt['variant_locus'])
            geno.addAffectedLocus(variant_bnode, gene_id)
            model.addBlankNodeAnnotation(variant_bnode)

            self._add_variant_attributes(
                variant_bnode, functional_consequence, cell_origin)

            gene_or_variant = variant_bnode

        else:
            gene_or_variant = gene_id

        assoc = G2PAssoc(
            self.graph, self.name, gene_or_variant, disease_id, g2p_relation)
        assoc.add_evidence(eco_id)
        assoc.add_association_to_graph()

        return

    def _add_variant_attributes(
            self,
            variant_id,
            functional_consequence=None,
            cell_origin=None):
        """
        Add attributes to variant

        :param variant_id:
        :param functional_consequence_gt: {str} global term for
                                          functional consequence
        :param cell_origin_gt: global term for cell origin
        :return: None
        """
        model = Model(self.graph)

        if functional_consequence is not None:
            predicate = self.globaltt['has_functional_consequence']
            model.addTriple(variant_id, predicate, functional_consequence)

        if cell_origin is not None:
            predicate = self.globaltt['has_cell_origin']
            model.addTriple(variant_id, predicate, cell_origin)

        return

    # def getTestSuite(self):
    #    import unittest
    #    # TODO PYLINT Unable to import 'tests.test_orphanet  (there is none)
    #    from tests.test_orphanet import OrphanetTestCase
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(OrphanetTestCase)
    #
    #    return test_suite
