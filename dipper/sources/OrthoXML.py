import logging
import itertools
import lxml.etree
import os

from dipper.sources.Source import Source
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Model import Model
from dipper import config

__author__ = "Adrian Altenhoff"

logger = logging.getLogger(__name__)


class OrthoXMLParser(object):
    def __init__(self, xml):
        if not (hasattr(xml, 'docinfo') and xml.docinfo.root_name == "orthoXML"):
            raise ValueError('Expecting an lxml.etree object of an orthoXML file as input')
        self.obj = xml
        self.relations = []
        self.gene_mapping = self._build_gene_mapping()

    def _build_gene_mapping(self):
        gene_map = {}
        for gene in self.obj.findall(".//{http://orthoXML.org/2011/}gene"):
            gene_map[gene.get('id')] = gene
        return gene_map

    def extract_pairwise_relations(self, node=None):
        if node is None:
            root_nodes = self.default_node_list()
        else:
            root_nodes = [node]

        for group in root_nodes:
            self.relations = []
            self._extract_pw(group)
            for rel in self.relations:
                yield rel

    def _extract_pw(self, node):
        if self.is_leaf(node):
            return {self.leaf_label(node)}
        elif self.is_internal_node(node):
            nodes_of_children = [self._extract_pw(child) for child in self.get_children(node)]
            rel = lxml.etree.QName(node).localname
            for child1, child2 in itertools.combinations(nodes_of_children, 2):
                for gId1, gId2 in itertools.product(child1, child2):
                    self.relations.append((gId1, gId2, rel))
            nodes = set.union(*nodes_of_children)
            return nodes
        else:
            return set([])

    def default_node_list(self):
        return self.obj.find(".//{http://orthoXML.org/2011/}groups")

    def is_internal_node(self, node):
        return lxml.etree.QName(node).localname in ('orthologGroup', 'paralogGroup')

    def leaf_label(self, leaf):
        return leaf.get('id')

    def is_leaf(self, node):
        return lxml.etree.QName(node).localname == "geneRef"

    def get_children(self, node):
        return node


class OrthoXML(Source):
    """
    Extract the induced pairwise relations from an OrthoXML file.

    This base class is primarily intended to extract the orthologous
    and paralogous relations from a file in OrthoXML file containing
    the QfO reference species data set.

    A concreate method should subclass this class and overwrite the
    constructor method to provide the information about the dataset
    and a method name.
    """
    files = {}

    def __init__(self, graph_type, are_bnodes_skolemized, method, tax_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, method)
        self.tax_ids = tax_ids
        self._map_orthology_code_to_RO = {
            'orthologGroup': OrthologyAssoc.ortho_rel['orthologous'],
            'paralogGroup': OrthologyAssoc.ortho_rel['paralogous']}

        if 'test_ids' not in config.get_config() \
                or 'protein' not in config.get_config()['test_ids']:
            logger.warning("not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['protein']

        return

    def fetch(self, is_dl_forced=False):
        """
        :return: None

        """
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        """
        :return: None
        """

        if self.testOnly:
            self.testMode = True

        if self.tax_ids is None:
            logger.info(
                "No taxon filter set; Dumping all orthologous associations.")
        else:
            logger.info(
                "Only the following taxa will be dumped: %s",
                str(self.tax_ids))

        self._get_orthologs(limit)

        return

    def _get_orthologs(self, limit):
        """
        This will process each of the specified pairwise orthology files,
        creating orthology associations based on the specified orthology code.
        this currently assumes that each of the orthology files is identically
        formatted. Relationships are made between genes here.

        There is also a nominal amount of identifier re-formatting:
        MGI:MGI --> MGI
        Ensembl --> ENSEMBL

        we skip any genes where we don't know how to map the gene identifiers.
        For example, Gene:Huwe1 for RAT is not an identifier, so we skip any
        mappings to this identifier.  Often, the there are two entries for the
        same gene (base on equivalent Uniprot id), and so we are not actually
        losing any information.

        We presently have a hard-coded filter to select only orthology
        relationships where one of the pair is in our species of interest
        (Mouse and Human, for the moment).
        This will be added as a configurable parameter in the future.

        Genes are also added to a grouping class defined with a PANTHER id.

        Triples:
        <gene1_id> RO:othologous <gene2_id>
        <assoc_id> :hasSubject <gene1_id>
        <assoc_id> :hasObject <gene2_id>
        <assoc_id> :hasPredicate <RO:orthologous>
        <assoc_id> dc:evidence ECO:phylogenetic_evidence

        <panther_id> a DATA:gene_family
        <panther_id> RO:has_member <gene1_id>
        <panther_id> RO:has_member <gene2_id>

        :param limit:
        :return:

        """
        logger.info("getting orthologs and paralog relations")

        g = self.testgraph if self.testMode else self.graph
        model = Model(g)

        for k in self.files.keys():
            f = os.path.join(self.rawdir, self.files[k]['file'])
            matchcounter = 0
            logger.info("Parsing %s", f)

            xml = lxml.etree.parse(f)
            parser = OrthoXMLParser(xml)

            for gene_nr_a, gene_nr_b, rel_type in parser.extract_pairwise_relations():
                gene_a = parser.gene_mapping[gene_nr_a]
                gene_b = parser.gene_mapping[gene_nr_b]

                gene_id_a = gene_a.get('protId')
                gene_id_b = gene_b.get('protId')

                if self.testMode and not \
                        (gene_id_a in self.test_ids or gene_id_b in self.test_ids):
                    continue

                matchcounter += 1
                taxon_a = self.extract_taxon_info(gene_a)
                taxon_b = self.extract_taxon_info(gene_b)

                rel = self._map_orthology_code_to_RO[rel_type]
                evidence_id = 'ECO:0000080'  # phylogenetic evidence

                # add genes to graph;
                # assume labels will be taken care of elsewhere
                gene_id_a = "UniProtKB:" + gene_id_a
                gene_id_b = "UniProtKB:" + gene_id_b
                model.addClassToGraph(gene_id_a, None)
                model.addClassToGraph(gene_id_b, None)

                # add the association and relevant nodes to graph
                assoc = OrthologyAssoc(g, self.name, gene_id_a, gene_id_b, rel)
                assoc.add_evidence(evidence_id)

                # might as well add the taxon info for completeness
                g.addTriple(
                    gene_id_a, model.object_properties['in_taxon'], taxon_a)
                g.addTriple(
                    gene_id_b, model.object_properties['in_taxon'], taxon_b)

                assoc.add_association_to_graph()

                if not self.testMode \
                        and limit is not None and matchcounter > limit:
                    logger.warning("reached limit of relations to extract. Stopping early...")
                    break
                    # make report on unprocessed_gene_ids

            logger.info("finished processing %s", f)
        return

    def extract_taxon_info(self, gene_node):
        """extract the ncbi taxon id from a gene_node

        default implementation goes up to the species node in the xml
        and extracts the id from the attribute at that node.
        """
        species_node = gene_node.getparent().getparent().getparent()
        return "NCBITaxon:{}".format(species_node.get('NCBITaxId'))

    def clean_protein_id(self, protein_id):
        """makes sure protein_id is properly prefixed"""
        if protein_id.find(':') > 0:
            return protein_id

        if protein_id.startswith('ENS'):
            return "ENSEMBL:" + protein_id
        elif protein_id.startswith('FB'):
            return "FlyBase:" + protein_id
        else:
            return "UniProtKB:" + protein_id
