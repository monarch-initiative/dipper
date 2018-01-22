import logging
import itertools
import re
import functools
import lxml.etree
import os

from dipper.sources.Source import Source
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Model import Model
from dipper.models.Genotype import Genotype
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
    # regex to match UniProtKB Identifiers. taken from https://www.uniprot.org/help/accession_numbers
    _up_re = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')

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

        self._get_relations(limit)

        return

    def _get_relations(self, limit):
        """
        This will process each of the specified orthoxml files, and
        extracting the induced orthology and paralogy associations
        based on the specified xml group nodes.

        The specs for orthoxml can be found here: http://orthoxml.org

        We currently extract tripples for orthologous relations,
        paralogous relations and in_taxon relations to NCBITaxonId
        attributes, e.g.

        Triples:
        <protein1_id> RO:othologous <protein2_id>
        <assoc_id> :hasSubject <protein1_id>
        <assoc_id> :hasObject <protein2_id>
        <assoc_id> :hasPredicate <RO:orthologous>
        <assoc_id> dc:evidence ECO:phylogenetic_evidence

        :param limit: limit the number of induced pairwise relations
        :return: None

        """
        logger.info("getting ortholog and paralog relations")

        g = self.testgraph if self.testMode else self.graph
        model = Model(g)

        for k in self.files.keys():
            f = os.path.join(self.rawdir, self.files[k]['file'])
            matchcounter = 0
            logger.info("Parsing %s", f)

            xml = lxml.etree.parse(f)
            parser = OrthoXMLParser(xml)

            for protein_nr_a, protein_nr_b, rel_type in parser.extract_pairwise_relations():
                protein_a = parser.gene_mapping[protein_nr_a]
                protein_b = parser.gene_mapping[protein_nr_b]

                protein_id_a = protein_a.get('protId')
                protein_id_b = protein_b.get('protId')

                if self.testMode and not \
                        (protein_id_a in self.test_ids or protein_id_b in self.test_ids):
                    continue

                matchcounter += 1
                taxon_a = self.extract_taxon_info(protein_a)
                taxon_b = self.extract_taxon_info(protein_b)

                # check if both protein belong to taxa that are selected
                if (self.tax_ids is not None and
                    (int(re.sub(r'NCBITaxon:', '', taxon_a.rstrip()))
                            not in self.tax_ids) and
                    (int(re.sub(r'NCBITaxon:', '', taxon_b.rstrip()))
                            not in self.tax_ids)):
                        continue

                protein_id_a = self.clean_protein_id(protein_id_a)
                protein_id_b = self.clean_protein_id(protein_id_b)
                # add genes to graph if needed;
                # assume labels will be taken care of elsewhere
                self.add_protein_to_graph(protein_id_a, taxon_a, model)
                self.add_protein_to_graph(protein_id_b, taxon_b, model)

                rel = self._map_orthology_code_to_RO[rel_type]
                evidence_id = 'ECO:0000080'  # phylogenetic evidence
                # add the association and relevant nodes to graph
                assoc = OrthologyAssoc(g, self.name, protein_id_a, protein_id_b, rel)
                assoc.add_evidence(evidence_id)
                assoc.add_association_to_graph()

                if not self.testMode \
                        and limit is not None and matchcounter > limit:
                    logger.warning("reached limit of relations to extract. Stopping early...")
                    break
                    # make report on unprocessed_gene_ids

            logger.info("finished processing %s", f)
        return

    @functools.lru_cache(2**15)
    def add_protein_to_graph(self, protein_id, taxon, model):
        """adds protein nodes to the graph and adds a "in_taxon" triple.

        for efficency reasons, we cache which proteins we have already
        added using a least recently used cache."""
        model.addClassToGraph(protein_id, None, Genotype.genoparts['polypeptide'])
        model.graph.addTriple(
            protein_id, model.object_properties['in_taxon'], taxon)

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
        elif self._up_re.match(protein_id):
            return "UniProtKB:" + protein_id
        else:
            logger.warning("namespace of protein id is not known: "+protein_id)
            return "UNKNOWN:" + protein_id
