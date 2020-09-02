import tarfile
import re
import logging

from dipper.sources.Source import Source
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Model import Model
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv


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
    "protein" identifiers in resources/test_id.yaml.

    By default, this will produce a file with ALL orthologous relationships.
    IF YOU WANT ONLY A SUBSET, YOU NEED TO PROVIDE A FILTER UPON CALLING THIS
    WITH THE TAXON IDS

    """
    PNTHDL = 'ftp://ftp.pantherdb.org/ortholog/current_release'
    # see: ftp://ftp.pantherdb.org/ortholog/current_release/README
    panther_format = [
        'Gene',                                 # species1|DB=id1|protdb=pdbid1
        'Ortholog',                             # species2|DB=id2|protdb=pdbid2
        'Type of ortholog',                     # [LDO, O, P, X ,LDX]  see: localtt
        'Common ancestor for the orthologs',    # unused
        'Panther Ortholog ID'                   # panther_id
    ]

    files = {
        'RefGenomeOrthologs': {
            'file': 'RefGenomeOrthologs.tar.gz',
            'url': PNTHDL + '/RefGenomeOrthologs.tar.gz',
            'columns': panther_format
        },
        'Orthologs_HCOP': {  # HUGO Comparison of Orthology Prediction
            'file': 'Orthologs_HCOP.tar.gz',
            'url': PNTHDL + '/Orthologs_HCOP.tar.gz',
            'columns': panther_format
        },
        'current_release': {
            'file': 'current_release.ver',               # preprocessed via DipperCache
            'url': 'ftp://ftp.pantherdb.org/ortholog/',  # symlink in the ftp listing
            'columns': ['version']
        }
    }

    def __init__(
            self,
            graph_type,
            are_bnodes_skolemized,
            data_release_version=None,
            tax_ids=None
    ):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
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
        # TODO the version number is in 'current_release.ver'; where to put it now?

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

        for src_key in ['RefGenomeOrthologs', 'Orthologs_HCOP']:
            self._get_orthologs(src_key, limit)

    def _get_orthologs(self, src_key, limit):
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

        We presently have a filter to select only orthology relationships where
        each of the pair is found in self.tax_ids.

        Genes are also added to a grouping class defined with a PANTHER id.

        Triples:
        <gene1_id> RO:othologous <gene2_id>
        <assoc_id> :hasSubject <gene1_id>
        <assoc_id> :hasObject <gene2_id>
        <assoc_id> :hasPredicate <RO:orthologous>
        <assoc_id> dcterms:evidence ECO:phylogenetic_evidence

        <panther_id> rdf:type DATA:gene_family
        <panther_id> RO:has_member <gene1_id>
        <panther_id> RO:has_member <gene2_id>

        :param limit:
        :return:

        """
        LOG.info("reading orthologs")

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)

        unprocessed_gene_ids = []

        src_file = '/'.join((self.rawdir, self.files[src_key]['file']))
        matchcounter = line_counter = 0
        col = self.files[src_key]['columns']
        reader = tarfile.open(src_file, 'r:gz')

        LOG.info("Parsing %s", src_key)

        with reader.extractfile(src_key) as csvfile:
            # there are no comments or headers
            for line in csvfile:
                # a little feedback to the user since there's so many ... bah strace
                # if line_counter % 1000000 == 0:
                #    LOG.info("Processed %d lines from %s", line_counter, fname.name)

                # parse each row. ancestor_taxons is unused
                # HUMAN|Ensembl=ENSG00000184730|UniProtKB=Q0VD83
                #   	MOUSE|MGI=MGI=2176230|UniProtKB=Q8VBT6
                #       	LDO	Euarchontoglires	PTHR15964

                row = line.decode().split('\t')
                thing1 = row[col.index('Gene')].strip()
                thing2 = row[col.index('Ortholog')].strip()
                orthology_type = row[col.index('Type of ortholog')].strip()
                # ancestor_taxons  = row[
                #    col.index('Common ancestor for the orthologs')].strip()
                panther_id = row[
                    col.index('Panther Ortholog ID')].strip()

                (species_a, gene_a, protein_a) = thing1.split('|')
                (species_b, gene_b, protein_b) = thing2.split('|')

                # for testing skip entries without homolog relationships to test ids
                if self.test_mode and not (
                        protein_a[9:] in self.test_ids or
                        protein_b[9:] in self.test_ids):
                    continue

                # map the species abbreviations to ncbi taxon id numbers
                taxon_a = self.resolve(species_a).split(':')[1].strip()
                taxon_b = self.resolve(species_b).split(':')[1].strip()

                # ###
                # keep orthologous relationships to genes in the given tax_ids
                # using AND will get you only those associations where
                # gene1 AND gene2 are in the taxid list (most-filter)
                # using OR will get you any associations where
                # gene1 OR gene2 are in the taxid list (some-filter)
                if self.tax_ids is not None and (
                        taxon_a not in self.tax_ids) and (
                        taxon_b not in self.tax_ids):
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
                    unprocessed_gene_ids.append(gene_a)
                    continue
                gene_a = clean_gene
                clean_gene = self._clean_up_gene_id(gene_b, species_b)
                if clean_gene is None:
                    unprocessed_gene_ids.append(gene_b)
                    continue
                gene_b = clean_gene

                rel = self.resolve(orthology_type)

                evidence_id = self.globaltt['phylogenetic evidence']

                # add the association and relevant nodes to graph
                assoc = OrthologyAssoc(graph, self.name, gene_a, gene_b, rel)
                assoc.add_evidence(evidence_id)

                # add genes to graph;  assume labels will be taken care of elsewhere
                model.addType(gene_a, self.globaltt['gene'])
                model.addType(gene_b, self.globaltt['gene'])

                # might as well add the taxon info for completeness
                graph.addTriple(
                    gene_a, self.globaltt['in taxon'], 'NCBITaxon:' + taxon_a
                )
                graph.addTriple(
                    gene_b, self.globaltt['in taxon'], 'NCBITaxon:' + taxon_b
                )

                assoc.add_association_to_graph(
                    blv.terms['GeneToGeneHomologyAssociation']
                )

                # note this is incomplete...
                # it won't construct the full family hierarchy,
                # just the top-grouping
                assoc.add_gene_family_to_graph('PANTHER:' + panther_id)

                if not self.test_mode and\
                        limit is not None and line_counter > limit:
                    break

            LOG.info("finished processing %s", src_file)
            LOG.warning(
                "The following gene ids were unable to be processed: %s",
                str(set(unprocessed_gene_ids)))

    def _clean_up_gene_id(self, geneid, species):
        """
        A series of identifier rewriting to conform with
        standard gene identifiers.
        :param geneid:
        :param species:
        :return:
        """
        # no curie should have more than one colon. generalize as
        geneid = ':'.join(geneid.split(':')[-2:])
        # which keeps the penultimate & ultimite tokens (aka last two)

        if species == 'CAEEL':
            if geneid[0:14] == 'EnsemblGenome:':
                geneid = 'WormBase:' + geneid[14:]
            elif geneid[:9] == 'Gene:CELE':  # Gene:CELE --> WormBase: (old cosmid ids)
                geneid = 'WormBase:' + geneid[9:]
            elif geneid[0:5] == 'Gene:':
                geneid = 'WormBase' + geneid[5:]

        elif species == 'DROME':
            if geneid[0:14] == 'EnsemblGenome:':
                geneid = 'FlyBase:' + geneid[14:]
            elif geneid[:10] == 'Gene:Dmel_':   # Gene:Dmel_ --> FlyBase:
                geneid = 'FlyBase:' + geneid[10:]
            elif geneid[:7] == 'Gene:CG':
                geneid = 'FlyBase:' + geneid[5:]  # Gene:CG --> FlyBase:CG

        else:

            if geneid[:8] == 'Ensembl:':         # Ensembl --> ENSEMBL   (why?)
                geneid = 'ENSEMBL:' + geneid[8:]
            elif geneid[:7] == 'GeneID:':        # GeneID: --> NCBIGene:
                geneid = 'NCBIGene:' + geneid[7:]
            elif geneid[:8] == 'Gene:ENS':       # Gene:<ensembl ids> --> ENSEMBL:<id>
                geneid = 'ENSEMBL:' + geneid[5:]
            elif geneid[:12] == 'Gene:Xenbase':  # Gene:<Xenbase ids> --> Xenbase:<id>
                geneid = 'Xenbase:' + geneid[12:]

        pfx = geneid.split(':')[0]
        if pfx is None or pfx not in self.curie_map:
            # LOG.warning( "No curie prefix for (species %s): %s", species, geneid)
            geneid = None
        return geneid

    def getTestSuite(self):
        import unittest
        from tests.test_panther import PantherTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(
            PantherTestCase)

        return test_suite
