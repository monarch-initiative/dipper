import tarfile
import re
import logging

from dipper.sources.Source import Source
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
from dipper import config
from dipper import curie_map

__author__ = 'nicole'

logger = logging.getLogger(__name__)


class Panther(Source):
    """
    The pairwise orthology calls from Panther DB:
    http://pantherdb.org/ encompass 22 species,
    from the RefGenome and HCOP projects.
    Here, we map the orthology classes to RO homology relationships
    This resource may be extended in the future with additional species.

    This currently makes a graph of orthologous relationships between genes,
    with the assumption that gene metadata (labels, equivalent ids) are
    provided from other sources.

    Gene families are nominally created from the orthology files,
    though these are incomplete with no hierarchical (subfamily) information.
    This will get updated from the HMM files in the future.

    Note that there is a fair amount of identifier cleanup performed to align
    with our standard CURIE prefixes.

    The test graph of data is output based on configured
    "protein" identifiers in conf.json.

    By default, this will produce a file with ALL orthologous relationships.
    IF YOU WANT ONLY A SUBSET, YOU NEED TO PROVIDE A FILTER UPON CALLING THIS
    WITH THE TAXON IDS

    """
    PNTHDL = 'ftp://ftp.pantherdb.org/ortholog/current_release'
    files = {
        'refgenome': {
            'file': 'RefGenomeOrthologs.tar.gz',
            'url': PNTHDL+'/RefGenomeOrthologs.tar.gz'},
        'hcop': {
            'file': 'Orthologs_HCOP.tar.gz',
            'url': PNTHDL+'/Orthologs_HCOP.tar.gz'}
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'panther')
        self.tax_ids = tax_ids

        self.dataset = Dataset(
            'panther', 'Protein ANalysis THrough Evolutionary Relationships',
            'http://pantherdb.org/', None,
            'http://www.pantherdb.org/terms/disclaimer.jsp')

        # # Defaults
        # if self.tax_ids is None:
        #     self.tax_ids = [9606, 10090, 7955]

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
        # TODO the version number is tricky to get
        # we can't get it from redirects of the url
        # TODO use the remote timestamp of the file?

        return

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
        logger.info("getting orthologs")

        if self.testMode:
            g = self.testgraph

        else:
            g = self.graph
        model = Model(g)
        unprocessed_gene_ids = set()

        for k in self.files.keys():
            f = '/'.join((self.rawdir, self.files[k]['file']))
            matchcounter = 0
            mytar = tarfile.open(f, 'r:gz')

            # assume that the first entry is the item
            fname = mytar.getmembers()[0]
            logger.info("Parsing %s", fname.name)
            line_counter = 0
            with mytar.extractfile(fname) as csvfile:
                for line in csvfile:
                    # skip comment lines
                    if re.match(r'^#', line.decode()):
                        logger.info("Skipping header line")
                        continue
                    line_counter += 1

                    # a little feedback to the user since there's so many
                    if line_counter % 1000000 == 0:
                        logger.info(
                            "Processed %d lines from %s",
                            line_counter, fname.name)

                    line = line.decode().strip()

                    # parse each row. ancestor_taxon is unused
                    # HUMAN|Ensembl=ENSG00000184730|UniProtKB=Q0VD83
                    #   	MOUSE|MGI=MGI=2176230|UniProtKB=Q8VBT6
                    #       	LDO	Euarchontoglires	PTHR15964
                    (a, b, orthology_class, ancestor_taxon,
                     panther_id) = line.split('\t')
                    (species_a, gene_a, protein_a) = a.split('|')
                    (species_b, gene_b, protein_b) = b.split('|')

                    # skip the entries that don't have homolog relationships
                    # with the test ids
                    if self.testMode and not (
                            re.sub(r'UniProtKB=', '',
                                   protein_a) in self.test_ids or
                            re.sub(r'UniProtKB=', '', protein_b)
                            in self.test_ids):
                        continue

                    # map the taxon abbreviations to ncbi taxon ids
                    taxon_a = self._map_taxon_abbr_to_id(species_a)
                    taxon_b = self._map_taxon_abbr_to_id(species_b)

                    # ###uncomment the following code block
                    # if you want to filter based on taxid of favorite animals
                    # taxids = [9606,10090,10116,7227,7955,6239,8355]
                    # taxids = [9606] #human only
                    # retain only those orthologous relationships to genes
                    # in the specified taxids
                    # using AND will get you only those associations where
                    # gene1 AND gene2 are in the taxid list (most-filter)
                    # using OR will get you any associations where
                    # gene1 OR gene2 are in the taxid list (some-filter)
                    if (
                        self.tax_ids is not None and
                        (int(re.sub(r'NCBITaxon:', '', taxon_a.rstrip()))
                            not in self.tax_ids) and
                        (int(re.sub(
                            r'NCBITaxon:', '', taxon_b.rstrip())) not in
                            self.tax_ids)):
                        continue
                    else:
                        matchcounter += 1
                        if limit is not None and matchcounter > limit:
                            break

                    # ### end code block for filtering on taxon

                    # fix the gene identifiers
                    gene_a = re.sub(r'=', ':', gene_a)
                    gene_b = re.sub(r'=', ':', gene_b)

                    clean_gene = self._clean_up_gene_id(gene_a, species_a)
                    if clean_gene is None:
                        unprocessed_gene_ids.add(gene_a)
                    gene_a = clean_gene
                    clean_gene = self._clean_up_gene_id(gene_b, species_b)
                    if clean_gene is None:
                        unprocessed_gene_ids.add(gene_b)
                    gene_b = clean_gene

                    # a special case here; mostly some rat genes
                    # they use symbols instead of identifiers.  will skip
                    if gene_a is None or gene_b is None:
                        continue

                    rel = self._map_orthology_code_to_RO(orthology_class)

                    evidence_id = 'ECO:0000080'  # phylogenetic evidence

                    # add the association and relevant nodes to graph
                    assoc = OrthologyAssoc(g, self.name, gene_a, gene_b, rel)
                    assoc.add_evidence(evidence_id)

                    # add genes to graph;
                    # assume labels will be taken care of elsewhere
                    model.addClassToGraph(gene_a, None)
                    model.addClassToGraph(gene_b, None)

                    # might as well add the taxon info for completeness
                    g.addTriple(
                        gene_a, model.object_properties['in_taxon'], taxon_a)
                    g.addTriple(
                        gene_b, model.object_properties['in_taxon'], taxon_b)

                    assoc.add_association_to_graph()

                    # note this is incomplete...
                    # it won't construct the full family hierarchy,
                    # just the top-grouping
                    assoc.add_gene_family_to_graph(
                        ':'.join(('PANTHER', panther_id)))

                    if not self.testMode \
                            and limit is not None and line_counter > limit:
                        break

            logger.info("finished processing %s", f)
            logger.warning(
                "The following gene ids were unable to be processed: %s",
                str(unprocessed_gene_ids))

        return

    @staticmethod
    def _map_taxon_abbr_to_id(ptax):
        """
        Will map the panther-specific taxon abbreviations to NCBI taxon numbers
        :param ptax:
        :return: NCBITaxon id
        """
        taxid = None
        ptax_to_taxid_map = {
            'ANOCA': 28377,  # green lizard
            'ARATH': 3702,   # arabadopsis
            'BOVIN': 9913,   # cow
            'CAEEL': 6239,   # worm
            'CANFA': 9615,   # dog
            'CHICK': 9031,   # chicken
            'DANRE': 7955,   # zebrafish
            'DICDI': 44689,  # discodium
            'DROME': 7227,   # drosophila melanogaster
            'ECOLI': 562,
            'HORSE': 9796,   # horses
            'HUMAN': 9606,   # humans
            'MACMU': 9544,   # macaque
            'MONDO': 13616,  # opossum
            'MOUSE': 10090,  # mouse
            'ORNAN': 9258,   # orangutan
            'PANTR': 9598,   # chimp
            'PIG': 9823,
            'RAT': 10116,
            'SCHPO': 4896,   # pombe yeast
            'TAKRU': 31033,  # pufferfish
            'XENTR': 8364,   # xenopus
            'YEAST': 4932,   # yeast
        }

        if ptax in ptax_to_taxid_map:
            taxid = ':'.join(('NCBITaxon', str(ptax_to_taxid_map.get(ptax))))
        else:
            logger.error("unmapped taxon code %s", ptax)

        return taxid

    @staticmethod
    def _map_orthology_code_to_RO(ortho):
        """
        Map the panther-specific orthology code (P,O,LDO,X,LDX)
        to relationship-ontology
        identifiers.
        :param ortho: orthology code
        :return: RO identifier
        """
        ortho_rel = OrthologyAssoc.ortho_rel
        ro_id = ortho_rel['orthologous']  # in orthology relationship with
        ortho_to_ro_map = {
            'P': ortho_rel['paralogous'],
            'O': ortho_rel['orthologous'],
            'LDO': ortho_rel['least_diverged_orthologous'],
            'X': ortho_rel['xenologous'],
            'LDX': ortho_rel['xenologous']
        }

        if ortho in ortho_to_ro_map:
            ro_id = ortho_to_ro_map.get(ortho)
        else:
            logger.warning(
                "unmapped orthology code %s. Defaulting to 'orthology'", ortho)

        return ro_id

    @staticmethod
    def _clean_up_gene_id(geneid, sp):
        """
        A series of identifier rewriting to conform with
        standard gene identifiers.
        :param geneid:
        :param sp:
        :return:
        """
        # special case for MGI
        geneid = re.sub(r'MGI:MGI:', 'MGI:', geneid)

        # rewrite Ensembl --> ENSEMBL
        geneid = re.sub(r'Ensembl', 'ENSEMBL', geneid)

        # rewrite Gene:CELE --> WormBase
        # these are old-school cosmid identifier
        geneid = re.sub(r'Gene:CELE', 'WormBase:', geneid)
        if sp == 'CAEEL':
            if re.match(r'(Gene|ENSEMBLGenome):\w+\.\d+', geneid):
                geneid = re.sub(
                    r'(?:Gene|ENSEMBLGenome):(\w+\.\d+)',
                    r'WormBase:\1', geneid)

        if sp == 'DROME':
            if re.match(r'(ENSEMBLGenome):\w+\.\d+', geneid):
                geneid = re.sub(
                    r'(?:ENSEMBLGenome):(\w+\.\d+)', r'FlyBase:\1', geneid)

        # rewrite GeneID --> NCBIGene
        geneid = re.sub(r'GeneID', 'NCBIGene', geneid)

        # rewrite Gene:Dmel --> FlyBase
        geneid = re.sub(r'Gene:Dmel_', 'FlyBase:', geneid)
        # rewrite Gene:CG --> FlyBase:CG
        geneid = re.sub(r'Gene:CG', 'FlyBase:CG', geneid)

        # rewrite ENSEMBLGenome:FBgn --> FlyBase:FBgn
        geneid = re.sub(r'ENSEMBLGenome:FBgn', 'FlyBase:FBgn', geneid)

        # rewrite Gene:<ensembl ids> --> ENSEMBL:<id>
        geneid = re.sub(r'Gene:ENS', 'ENSEMBL:ENS', geneid)

        # rewrite Gene:<Xenbase ids> --> Xenbase:<id>
        geneid = re.sub(r'Gene:Xenbase:', 'Xenbase:', geneid)

        # TODO this would be much better done as
        # if foo not in curie_map:
        # if re.match(r'(Gene|ENSEMBLGenome):', geneid) or \
        #        re.match(r'Gene_ORFName', geneid) or \
        #        re.match(r'Gene_Name', geneid):
        #    # logger.warning(
        #    #"Found an identifier I don't know how to fix (species %s): %s",
        #    #   sp, geneid)

        pfxlcl = re.split(r':', geneid)
        pfx = pfxlcl[0]
        if pfx is None or pfx not in curie_map.get():
            logger.warning("No curie prefix for (species %s): %s", sp, geneid)
            geneid = None
        return geneid

    def getTestSuite(self):
        import unittest
        from tests.test_panther import PantherTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(
            PantherTestCase)

        return test_suite
