import tarfile
import re
import logging

from dipper.sources.Source import Source
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Model import Model

__author__ = 'nicole'

LOG = logging.getLogger(__name__)


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
    "protein" identifiers in conf.yaml.

    By default, this will produce a file with ALL orthologous relationships.
    IF YOU WANT ONLY A SUBSET, YOU NEED TO PROVIDE A FILTER UPON CALLING THIS
    WITH THE TAXON IDS

    """
    PNTHDL = 'ftp://ftp.pantherdb.org/ortholog/current_release'
    # see: ftp://ftp.pantherdb.org/ortholog/current_release/README
    files = {
        'refgenome': {
            'file': 'RefGenomeOrthologs.tar.gz',
            'url': PNTHDL+'/RefGenomeOrthologs.tar.gz',
            'columns': [  # expected
                'thing1',           # pipe-list species|DB=id|protdb=pdbid
                'thing2',           # pipe-list species|DB=id|protdb=pdbid
                'orthology_class',  # [LDO, O, P, X ,LDX]  see: localtt
                'ancestor_taxon',   # ancestor sometimes piped-pair  or 'ND'
                'panther_id'        # ortholog cluster id I assume
            ]
        },
        'hcop': {  # HUGO Comparison of Orthology Prediction
            'file': 'Orthologs_HCOP.tar.gz',
            'url': PNTHDL + '/Orthologs_HCOP.tar.gz',
            'columns': [  # expected
                'thing1',           # pipe-list species|DB=id|protdb=pdbid
                'thing2',           # pipe-list species|DB=id|protdb=pdbid
                'orthology_class',  # [LDO, O, P, X ,LDX]  see: localtt
                'ancestor_taxon',   # ancestor sometimes piped-pair  or 'ND'
                'panther_id'        # ortholog cluster id I assume
            ]
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized, data_release_version=None, tax_ids=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            name='panther',
            ingest_title='Protein ANalysis THrough Evolutionary Relationships',
            ingest_url='http://pantherdb.org',
            ingest_logo='source-panther.jpg',
            license_url=None,
            data_rights='http://www.pantherdb.org/terms/disclaimer.jsp'
            # file_handle=None
        )
        self.dataset.set_citation(
            'http://pantherdb.org/publications.jsp#HowToCitePANTHER')

        default_species = ['Homo sapiens', 'Mus musculus', 'Danio rerio']
        if tax_ids is None:
            # self.tax_ids = ['9606', '10090', '7955']
            self.tax_ids = [self.globaltt[x].split(':')[-1] for x in default_species]
        else:
            self.tax_ids = [str(x) for x in tax_ids]

        if 'protein' in self.all_test_ids:
            self.test_ids = self.all_test_ids['protein']
        else:
            LOG.warning("not configured with protein test ids.")
            self.test_ids = []

    def fetch(self, is_dl_forced=False):
        """
        :return: None

        """
        self.get_files(is_dl_forced)
        # TODO the version number is tricky to get
        # we can't get it from redirects of the url
        # TODO use the remote timestamp of the file?

    def parse(self, limit=None):
        """
        :return: None
        """

        if self.test_only:
            self.test_mode = True

        if self.tax_ids is None:
            LOG.info("No taxon filter set; Dumping all orthologous associations.")
        else:
            LOG.info("Only the following taxa will be dumped: %s", self.tax_ids)

        self._get_orthologs(limit)

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
        LOG.info("getting orthologs")

        if self.test_mode:
            graph = self.testgraph

        else:
            graph = self.graph
        model = Model(graph)
        unprocessed_gene_ids = set()  # may be faster to make a set after

        for src_key in self.files:
            src_file = '/'.join((self.rawdir, self.files[src_key]['file']))
            matchcounter = line_counter = 0
            col = self.files[src_key]['columns']
            reader = tarfile.open(src_file, 'r:gz')

            # assume that the first entry is the item
            fname = reader.getmembers()[0]
            LOG.info("Parsing %s", fname.name)

            with reader.extractfile(fname) as csvfile:
                for line in csvfile:
                    # skip comment lines
                    if re.match(r'^#', line.decode()):
                        LOG.info("Skipping header line")
                        continue
                    line_counter += 1

                    # a little feedback to the user since there's so many
                    if line_counter % 1000000 == 0:
                        LOG.info(
                            "Processed %d lines from %s",
                            line_counter, fname.name)

                    line = line.decode().strip()
                    row = line.split('\t')
                    # parse each row. ancestor_taxons is unused
                    # HUMAN|Ensembl=ENSG00000184730|UniProtKB=Q0VD83
                    #   	MOUSE|MGI=MGI=2176230|UniProtKB=Q8VBT6
                    #       	LDO	Euarchontoglires	PTHR15964
                    thing1 = row[col.index('thing1')].strip()
                    thing2 = row[col.index('thing2')].strip()
                    orthology_class = row[col.index('orthology_class')].strip()
                    # ancestor_taxons  = row[col.index('')].strip()
                    panther_id = row[col.index('panther_id')].strip()

                    (species_a, gene_a, protein_a) = thing1.split('|')
                    (species_b, gene_b, protein_b) = thing2.split('|')

                    # skip the entries that don't have homolog relationships
                    # with the test ids
                    if self.test_mode and not (
                            re.sub(r'UniProtKB=', '',
                                   protein_a) in self.test_ids or
                            re.sub(r'UniProtKB=', '', protein_b)
                            in self.test_ids):
                        continue

                    # map the taxon abbreviations to ncbi taxon id numbers
                    taxon_a = self.resolve(species_a).split(':')[1].strip()
                    taxon_b = self.resolve(species_b).split(':')[1].strip()

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
                    if self.tax_ids is not None and \
                            (taxon_a not in self.tax_ids) and \
                            (taxon_b not in self.tax_ids):
                        continue
                    else:
                        matchcounter += 1
                        if limit is not None and matchcounter > limit:
                            break

                    # ### end code block for filtering on taxon

                    # fix the gene identifiers
                    gene_a = re.sub(r'=', ':', gene_a)
                    gene_b = re.sub(r'=', ':', gene_b)

                    clean_gene = self._clean_up_gene_id(
                        gene_a, species_a, self.curie_map)
                    if clean_gene is None:
                        unprocessed_gene_ids.add(gene_a)
                    gene_a = clean_gene
                    clean_gene = self._clean_up_gene_id(
                        gene_b, species_b, self.curie_map)
                    if clean_gene is None:
                        unprocessed_gene_ids.add(gene_b)
                    gene_b = clean_gene

                    # a special case here; mostly some rat genes
                    # they use symbols instead of identifiers.  will skip
                    if gene_a is None or gene_b is None:
                        continue

                    rel = self.resolve(orthology_class)

                    evidence_id = self.globaltt['phylogenetic evidence']

                    # add the association and relevant nodes to graph
                    assoc = OrthologyAssoc(graph, self.name, gene_a, gene_b, rel)
                    assoc.add_evidence(evidence_id)

                    # add genes to graph;
                    # assume labels will be taken care of elsewhere
                    model.addClassToGraph(gene_a, None)
                    model.addClassToGraph(gene_b, None)

                    # might as well add the taxon info for completeness
                    graph.addTriple(
                        gene_a, self.globaltt['in taxon'], 'NCBITaxon:' + taxon_a)
                    graph.addTriple(
                        gene_b, self.globaltt['in taxon'], 'NCBITaxon:' + taxon_b)

                    assoc.add_association_to_graph()

                    # note this is incomplete...
                    # it won't construct the full family hierarchy,
                    # just the top-grouping
                    assoc.add_gene_family_to_graph(
                        ':'.join(('PANTHER', panther_id)))

                    if not self.test_mode \
                            and limit is not None and line_counter > limit:
                        break
                # make report on unprocessed_gene_ids

            LOG.info("finished processing %s", src_file)
            LOG.warning(
                "The following gene ids were unable to be processed: %s",
                str(unprocessed_gene_ids))

    @staticmethod
    def _clean_up_gene_id(geneid, species, curie_map):
        """
        A series of identifier rewriting to conform with
        standard gene identifiers.
        :param geneid:
        :param species:
        :return:
        """
        # special case for MGI
        geneid = re.sub(r'MGI:MGI:', 'MGI:', geneid)

        # rewrite Ensembl --> ENSEMBL
        geneid = re.sub(r'Ensembl', 'ENSEMBL', geneid)

        # rewrite Gene:CELE --> WormBase
        # these are old-school cosmid identifier
        geneid = re.sub(r'Gene:CELE', 'WormBase:', geneid)
        if species == 'CAEEL':
            if re.match(r'(Gene|ENSEMBLGenome):\w+\.\d+', geneid):
                geneid = re.sub(
                    r'(?:Gene|ENSEMBLGenome):(\w+\.\d+)',
                    r'WormBase:\1', geneid)

        if species == 'DROME':
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
        # if foo not in selfcurie_map:
        # if re.match(r'(Gene|ENSEMBLGenome):', geneid) or \
        #        re.match(r'Gene_ORFName', geneid) or \
        #        re.match(r'Gene_Name', geneid):
        #    # LOG.warning(
        #    #"Found an identifier I don't know how to fix (species %s): %s",
        #    #   species, geneid)

        pfxlcl = re.split(r':', geneid)
        pfx = pfxlcl[0]
        if pfx is None or pfx not in curie_map:
            # LOG.warning( "No curie prefix for (species %s): %s", species, geneid)
            geneid = None
        return geneid

    def getTestSuite(self):
        import unittest
        from tests.test_panther import PantherTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(
            PantherTestCase)

        return test_suite
