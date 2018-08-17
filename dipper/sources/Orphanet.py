import logging
import xml.etree.ElementTree as ET
from dipper.sources.Source import Source
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
            license_url='http://creativecommons.org/licenses/by-nd/3.0/',
            data_rights='http://omim.org/help/agreement'
            # file_handle=None
        )

        # check to see if there's any ids configured in the config;
        # otherwise, warn
        # TODO remove
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
            graph = self.testgraph
        else:
            graph = self.graph
        line_counter = 0

        geno = Genotype(graph)
        model = Model(graph)


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
                    gene_id = 'Orphanet:' + str(gene_num)
                    gene_type_id = self.resolved(gene_iid_to_type[gene_iid])
                    model.addClassToGraph(
                        gene_id, gene_symbol, gene_type_id, gene_name)
                    syn_list = a.find('./Gene/SynonymList')
                    if int(syn_list.get('count')) > 0:
                        for s in syn_list.findall('./Synonym'):
                            model.addSynonym(gene_id, s.text)

                    # IDs appear stable but removing for now  KS
                    #dgtype = a.find('DisorderGeneAssociationType').get('id')
                    #rel_id = self.resolve(dgtype)
                    dg_label = a.find('./DisorderGeneAssociationType/Name').text
                    #if rel_id is None:
                    #    logger.warning(
                    #        "Cannot map association type (%s) to RO " +
                    #        "for association (%s | %s).  Skipping.",
                    #        dg_label, disorder_label, gene_symbol)
                    #    continue

                    #alt_locus_id = '_:' + gene_num + '-' + disorder_num + 'VL'
                    #alt_label = ' '.join((
                    #    'some variant of', gene_symbol.strip(), disorder_label))

                    #model.addIndividualToGraph(
                    #    alt_locus_id, alt_label, self.globaltt['variant_locus'])
                    #geno.addAffectedLocus(alt_locus_id, gene_id)
                    #model.addBlankNodeAnnotation(alt_locus_id)

                    # consider typing the gain/loss-of-function variants like:
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002054
                    # http://sequenceontology.org/browser/current_svn/term/SO:0002053


                    # use "assessed" status to issue an evidence code
                    # FIXME I think that these codes are sub-optimal
                    status_code = a.find('DisorderGeneAssociationStatus').get('id')
                    # imported automatically asserted information
                    # used in automatic assertion
                    eco_id = self.globaltt[
                        'imported automatically asserted information used in automatic assertion']
                    # Assessed
                    # TODO are these internal ids stable between releases?
                    if status_code == '17991':
                        # imported manually asserted information
                        # used in automatic assertion
                        eco_id = self.globaltt[
                            'imported manually asserted information used in automatic assertion']
                    # Non-traceable author statement ECO_0000034
                    # imported information in automatic assertion ECO_0000313

                    #assoc = G2PAssoc(
                    #    graph, self.name, alt_locus_id, disorder_id, rel_id)
                    # assoc.add_evidence(eco_id)
                    # assoc.add_association_to_graph()

                    self.add_gene_to_disease(
                        dg_label, gene_id, gene_symbol, disorder_id, eco_id)

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


    def add_gene_to_disease(
            self,
            association_type,
            gene_id,
            gene_symbol,
            disease_id,
            eco_id):
        """
        Composes triples based on the DisorderGeneAssociationType
        element:

        xmlstarlet sel -t  -v "/JDBOR/DisorderList/Disorder/
        DisorderGeneAssociationList/DisorderGeneAssociation/
        DisorderGeneAssociationType/Name" en_product6.xml | sort -u

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
        gene-disease predicate (has_phenotype, contributes_to)
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

        # If we know something about the variant
        # such as functional consequence or cellular origin
        # make a blank node and attach the attributes
        is_variant = False
        variant_id_string = "{}{}".format(gene_id, disease_id)
        functional_consequence = None
        cell_origin = None

        # hard fail for no mappings/new terms,
        # otherwise they go unnoticed
        if "{}|gene phenotype".format(association_type) not in self.localtt:
            raise ValueError(
                'Disease-gene association type {} not mapped'
                    .format(association_type)
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
                variant_label = "{} {}"\
                    .format(
                        functional_consequence_lbl.replace('_', ' '),
                        variant_label
                )
                variant_id_string += functional_consequence_lbl
            else:
                variant_label = "variant {}".format(variant_label)

            if cell_origin:
                variant_label = "{} {}".format(cell_origin_lbl, variant_label)
                variant_id_string += cell_origin_lbl

            variant_bnode = self.make_id(variant_id_string, "_")
            model.addIndividualToGraph(variant_bnode, variant_label,
                                       geno.genoparts['variant_locus'])
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

    @staticmethod
    def _map_gene_type_id(orphanet_type_id):
        """
        To check the "GeneType" mappings
        xmlstarlet sel -t -v \
        './JDBOR/DisorderList/Disorder//Gene/GeneType/Name'
         en_product6.xml | sort | uniq -c | sort -nr

        :param orphanet_type_id:
        :return:
        """
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

        test_suite = unittest.TestLoader().loadTestsFromTestCase(OrphanetTestCase)

        return test_suite
