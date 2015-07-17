import csv
import logging
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.OrthologyAssoc import OrthologyAssoc
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Pathway import Pathway
from dipper import curie_map
from dipper import config

logger = logging.getLogger(__name__)


class KEGG(Source):
    files = {
        'disease': {'file': 'disease',
                    'url': 'http://rest.genome.jp/list/disease'},
        'pathway': {'file': 'pathway',
                    'url': 'http://rest.genome.jp/list/pathway'},
        'hsa_genes': {'file': 'hsa_genes',
                      'url': 'http://rest.genome.jp/list/hsa'},
        'ortholog_classes': {'file': 'ortholog_classes',
                             'url': 'http://rest.genome.jp/list/orthology'},
        'disease_gene': {'file': 'disease_gene',
                         'url': 'http://rest.kegg.jp/link/disease/hsa'},
        'omim2disease': {'file': 'omim2disease',
                         'url': 'http://rest.genome.jp/link/disease/omim'},
        'omim2gene': {'file': 'omim2gene',
                      'url': 'http://rest.genome.jp/link/omim/hsa'},
        'ncbi': {'file': 'ncbi',
                 'url': 'http://rest.genome.jp/conv/ncbi-geneid/hsa'},
        'hsa_gene2pathway': {'file': 'human_gene2pathway',
                             'url': 'http://rest.kegg.jp/link/pathway/hsa'},
        'hsa_orthologs': {'file': 'hsa_orthologs',
                          'url': 'http://rest.kegg.jp/link/orthology/hsa'},
        'mmu_orthologs': {'file': 'mmu_orthologs',
                          'url': 'http://rest.kegg.jp/link/orthology/mmu'},
        'rno_orthologs': {'file': 'rno_orthologs',
                          'url': 'http://rest.kegg.jp/link/orthology/rno'},
        'dme_orthologs': {'file': 'dme_orthologs',
                          'url': 'http://rest.kegg.jp/link/orthology/dme'},
        'dre_orthologs': {'file': 'dre_orthologs',
                          'url': 'http://rest.kegg.jp/link/orthology/dre'},
        'cel_orthologs': {'file': 'cel_orthologs',
                          'url': 'http://rest.kegg.jp/link/orthology/cel'}
    }

    # TODO http://rest.kegg.jp/link/pathway/pubmed
    # TODO http://rest.kegg.jp/link/pathway/ds  # disease pathway assoc

    test_ids = {
        "pathway": ["path:map00010", "path:map00195", "path:map00100", "path:map00340", "path:hsa05223"],
        "disease": ["ds:H00015", "ds:H00026", "ds:H00712", "ds:H00736", "ds:H00014"],
        "genes": ["hsa:100506275", "hsa:285958", "hsa:286410", "hsa:6387",
                  "hsa:1080", "hsa:11200", "hsa:1131", "hsa:1137", "hsa:126", "hsa:1277", "hsa:1278",
                  "hsa:1285", "hsa:1548", "hsa:1636", "hsa:1639", "hsa:183", "hsa:185", "hsa:1910",
                  "hsa:207", "hsa:2099", "hsa:2483", "hsa:2539", "hsa:2629", "hsa:2697", "hsa:3161",
                  "hsa:3845", "hsa:4137", "hsa:4591", "hsa:472", "hsa:4744", "hsa:4835", "hsa:4929",
                  "hsa:5002", "hsa:5080", "hsa:5245", "hsa:5290", "hsa:53630", "hsa:5630", "hsa:5663",
                  "hsa:580", "hsa:5888", "hsa:5972", "hsa:6311", "hsa:64327", "hsa:6531", "hsa:6647",
                  "hsa:672", "hsa:675", "hsa:6908", "hsa:7040", "hsa:7045", "hsa:7048", "hsa:7157",
                  "hsa:7251", "hsa:7490", "hsa:7517", "hsa:79728", "hsa:83893", "hsa:83990", "hsa:841",
                  "hsa:8438", "hsa:8493", "hsa:860", "hsa:9568", "hsa:9627", "hsa:9821", "hsa:999",
                  "hsa:3460"],
        "orthology_classes": ["ko:K00010", "ko:K00027", "ko:K00042", "ko:K00088"]
    }

    def __init__(self):
        Source.__init__(self, 'kegg')

        # update the dataset object with details about this resource
        self.dataset = Dataset('kegg', 'KEGG', 'http://www.genome.jp/kegg/', None, None)

        # source-specific warnings.  will be cleared when resolved.
        # check to see if there's any ids configured in the config; otherwise, warn
        if 'test_ids' not in config.get_config() or 'disease' not in config.get_config()['test_ids']:
            logger.warn("not configured with disease test ids.")
        else:
            self.test_ids['disease'] += config.get_config()['test_ids']['disease']

        return

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

        # TODO add versioning information from info rest call, like http://rest.kegg.jp/info/pathway

        return

    def parse(self, limit=None):
        """

        :param limit:
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %s rows fo each file", str(limit))

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True
        self.label_hash = {}
        self.omim_disease_hash = {}  # to hold the mappings of omim:kegg ids
        self.kegg_disease_hash = {}  # to hold the mappings of kegg:omim ids

        self._process_diseases(limit)
        self._process_genes(limit)
        self._process_genes_kegg2ncbi(limit)
        self._process_omim2gene(limit)
        self._process_omim2disease(limit)
        self._process_kegg_disease2gene(limit)

        self._process_pathways(limit)

        # self._process_ortholog_classes(limit)  #add in when refactoring for #141
        # for f in ['hsa_orthologs', 'mmu_orthologs', 'rno_orthologs','dme_orthologs','dre_orthologs','cel_orthologs']:
        #     file = '/'.join((self.rawdir, self.files[f]['file']))
        #     self._process_orthologs(file, limit)  # DONE #

        logger.info("Finished parsing")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return

    def _process_pathways(self, limit=None):
        """
        This method adds the KEGG pathway IDs.  These are the canonical pathways as defined in KEGG.
        We also encode the graphical depiction which maps 1:1 with the identifier.

        Triples created:
        <pathway_id> is a GO:signal_transduction
        <pathway_id> rdfs:label <pathway_name>
        <gene_id> RO:involved_in <pathway_id>
        :param limit:
        :return:
        """

        logger.info("Processing pathways")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        path = Pathway(g, self.nobnodes)
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['pathway']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pathway_id, pathway_name) = row

                if self.testMode and pathway_id not in self.test_ids['pathway']:
                    continue

                pathway_id = 'KEGG-'+pathway_id.strip()

                path.addPathway(pathway_id, pathway_name)

                # we know that the pathway images from kegg map 1:1 here.  so add those
                image_filename = re.sub('KEGG-path:','',pathway_id) + '.png'
                image_url = 'http://www.genome.jp/kegg/pathway/map/'+image_filename
                gu.addDepiction(g, pathway_id, image_url)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with pathways")
        return

    def _process_diseases(self, limit=None):
        """
        This method processes the KEGG disease IDs.

        Triples created:
        <disease_id> is a class
        <disease_id> rdfs:label <disease_name>
        :param limit:
        :return:
        """

        logger.info("Processing diseases")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (disease_id, disease_name) = row

                disease_id = 'KEGG-'+disease_id.strip()
                if disease_id not in self.label_hash:
                    self.label_hash[disease_id] = disease_name

                if self.testMode and disease_id not in self.test_ids['disease']:
                    continue

                # Add the disease as a class....we don't get all of these from MONDO yet
                # see https://github.com/monarch-initiative/human-disease-ontology/issues/3
                gu.addClassToGraph(g, disease_id, disease_name)
                # not typing the diseases as DOID:4 yet because I don't want to bulk up the graph unnecessarily

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with diseases")
        return

    def _process_genes(self, limit=None):
        """
        This method processes the KEGG gene IDs.
        The label for the gene is pulled as the first symbol in the list of gene symbols; the rest
        are added as synonyms.  The long-form of the gene name is added as a definition.
        This is hardcoded to just processes human genes.

        Triples created:
        <gene_id> is a SO:gene
        <gene_id> rdfs:label <gene_name>

        :param limit:
        :return:
        """

        logger.info("Processing genes")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['hsa_genes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, gene_name) = row

                gene_id = 'KEGG-'+gene_id.strip()

                # the gene listing has a bunch of labels that are delimited, like:
                # DST, BP240, BPA, BPAG1, CATX-15, CATX15, D6S1101, DMH, DT, EBSB2, HSAN6, MACF2; dystonin; K10382 dystonin
                # it looks like the list is semicolon delimited (symbol, name, gene_class)
                # where the symbol is a comma-delimited list

                # here, we split them up.  we will take the first abbreviation and make it the symbol
                # then take the rest as synonyms

                gene_stuff = re.split(';', gene_name)
                symbollist = re.split(',', gene_stuff[0])
                first_symbol = symbollist[0].strip()

                if gene_id not in self.label_hash:
                    self.label_hash[gene_id] = first_symbol

                if self.testMode and gene_id not in self.test_ids['genes']:
                    continue

                # Add the gene as a class.
                geno.addGene(gene_id, first_symbol)

                # add the long name as the description
                if len(gene_stuff) > 1:
                    description = gene_stuff[1].strip()
                    gu.addDefinition(g, gene_id, description)

                # add the rest of the symbols as synonyms
                for i in enumerate(symbollist, start=1):
                    gu.addSynonym(g, gene_id, i[1].strip())

                # TODO add the KO here?

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with genes")
        return

    def _process_ortholog_classes(self, limit=None):
        """
        This method add the KEGG orthology classes to the graph.

        Triples created:
        <orthology_class_id> is a class
        <orthology_class_id> has label <orthology_symbols>
        <orthology_class_id> has description <orthology_description>
        :param limit:
        :return:
        """

        logger.info("Processing ortholog classes")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['ortholog_classes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (orthology_class_id, orthology_class_name) = row

                if self.testMode and orthology_class_id not in self.test_ids['ortholog_classes']:
                    continue

                # FIXME: What's the proper route for this?
                # The orthology class is essentially a KEGG gene ID that is species agnostic.
                # Add the ID and label as a class. Would it be considered a gene as well?

                other_labels = re.split(';', orthology_class_name)
                orthology_label = other_labels[0]  # the first one is the label we'll use

                orthology_class_id = 'KEGG-'+orthology_class_id.strip()

                orthology_type = OrthologyAssoc.terms['gene_family']
                gu.addClassToGraph(g, orthology_class_id, orthology_label, orthology_type)
                if len(other_labels) > 1:
                    # add the rest as synonyms
                    # todo skip the first
                    for s in other_labels:
                        gu.addSynonym(g, orthology_class_id, s)

                    # add the last one as the description
                    gu.addDescription(g, orthology_class_id, other_labels[len(other_labels)-1])

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with ortholog classes")
        return

    def _process_orthologs(self, raw, limit=None):
        """
        This method maps orthologs for a species to the KEGG orthology classes.

        Triples created:
        <gene_id> is a class
        <orthology_class_id> is a class

        <assoc_id> has subject <gene_id>
        <assoc_id> has object <orthology_class_id>
        :param limit:
        :return:
        """

        logger.info("Processing orthologs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, orthology_class_id) = row

                orthology_class_id = 'KEGG:'+orthology_class_id.strip()
                gene_id = 'KEGG:'+gene_id.strip()

                # note that the panther_id references a group of orthologs, and is not 1:1 with the rest
                assoc_id = self.make_id(''.join((gene_id, orthology_class_id)))

                rel = OrthologyAssoc.ortho_rel['orthologous']
                # add the association and relevant nodes to graph
                assoc = OrthologyAssoc(assoc_id, gene_id, orthology_class_id, None, None)
                assoc.setRelationship(rel)
                assoc.loadAllProperties(g)

                # add gene and orthology class to graph; assume labels will be taken care of elsewhere
                gu.addClassToGraph(g, gene_id, None)
                gu.addClassToGraph(g, orthology_class_id, None)
                assoc.addAssociationToGraph(g)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with orthologs")
        return

    def _process_kegg_disease2gene(self, limit=None):
        """
        This method creates an association between diseases and their associated genes.
        We are being conservative here, and only processing those diseases for which there
        is no mapping to OMIM.

        Triples created:
        <alternate_locus> is an Individual
        <alternate_locus> has type <variant_locus>
        <alternate_locus> is an allele of  <gene_id>

        <assoc_id> has subject <disease_id>
        <assoc_id> has object <gene_id>
        :param limit:
        :return:
        """

        logger.info("Processing KEGG disease to gene")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())

        noomimset = set()
        raw = '/'.join((self.rawdir, self.files['disease_gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, disease_id) = row

                if self.testMode and gene_id not in self.test_ids['genes']:
                    continue

                gene_id = 'KEGG-'+gene_id.strip()
                disease_id = 'KEGG-'+disease_id.strip()

                # Make an association ID.
                rel = gu.object_properties['is_marker_for']
                assoc_id = self.make_association_id(self.name, gene_id, rel,
                                                    disease_id, None, None)

                # only add diseases for which there is no omim id and not a grouping class
                if disease_id not in self.kegg_disease_hash:
                    # add as a class
                    disease_label = None
                    if disease_id in self.label_hash:
                        disease_label = self.label_hash[disease_id]
                    if re.search('includ', str(disease_label)):
                        # they use 'including' when it's a grouping class
                        logger.info("Skipping this association because it's a grouping class: %s", disease_label)
                        continue
                    gu.addClassToGraph(g, disease_id, disease_label, 'DOID:4')  # type this disease_id as a disease
                    noomimset.add(disease_id)
                    alt_locus_id = self._make_variant_locus_id(gene_id, disease_id)
                    alt_label = self.label_hash[alt_locus_id]
                    gu.addIndividualToGraph(g, alt_locus_id, alt_label, geno.genoparts['variant_locus'])
                    geno.addAlleleOfGene(alt_locus_id, gene_id)
                    # Add the disease to gene relationship.
                    assoc = G2PAssoc(assoc_id, alt_locus_id, disease_id, None, None)
                    assoc.setRelationship(rel)
                    assoc.loadAllProperties(g)
                    assoc.addAssociationToGraph(g)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with KEGG disease to gene")
        logger.info("Found %d diseases with no omim id", len(noomimset))

        return

    def _process_omim2gene(self, limit=None):
        """
        This method maps the OMIM IDs and KEGG gene ID. Currently split based on the link_type field.
        Equivalent link types are mapped as gene XRefs.
        Reverse link types are mapped as disease to gene associations.
        Original link types are currently skipped.

        Triples created:
        <kegg_gene_id> is a Gene
        <omim_gene_id> is a Gene
        <kegg_gene_id>> hasXref <omim_gene_id>

        <assoc_id> has subject <omim_disease_id>
        <assoc_id> has object <kegg_gene_id>
        :param limit:
        :return:
        """

        logger.info("Processing OMIM to KEGG gene")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['omim2gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (kegg_gene_id, omim_id, link_type) = row

                if self.testMode and kegg_gene_id not in self.test_ids['genes']:
                    continue

                kegg_gene_id = 'KEGG-'+kegg_gene_id.strip()
                omim_id = re.sub('omim', 'OMIM', omim_id)
                if link_type == 'equivalent':
                    # these are genes!  so add them as a class then make equivalence
                    gu.addClassToGraph(g, omim_id, None)
                    geno.addGene(kegg_gene_id, None)
                    gu.addEquivalentClass(g, kegg_gene_id, omim_id)
                elif link_type == 'reverse':
                    # make an association between an OMIM ID and the KEGG gene ID
                    # we do this with omim ids because they are more atomic than KEGG ids

                    alt_locus_id = self._make_variant_locus_id(kegg_gene_id, omim_id)
                    alt_label = self.label_hash[alt_locus_id]
                    gu.addIndividualToGraph(g, alt_locus_id, alt_label, geno.genoparts['variant_locus'])
                    geno.addAlleleOfGene(alt_locus_id, kegg_gene_id)

                    # Add the disease to gene relationship.
                    rel = gu.object_properties['is_marker_for']
                    assoc_id = self.make_association_id(self.name, alt_locus_id, rel, omim_id, None, None)
                    assoc = G2PAssoc(assoc_id, alt_locus_id, omim_id, None, None)
                    assoc.setRelationship(rel)
                    assoc.addAssociationToGraph(g)
                    assoc.loadAllProperties(g)
                elif link_type == 'original':
                    # these are sometimes a gene, and sometimes a disease
                    logger.info('Unable to handle original link for %s-%s', kegg_gene_id, omim_id)
                else:
                    # don't know what these are
                    logger.warn('Unhandled link type for %s-%s: %s', kegg_gene_id, omim_id, link_type)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with OMIM to KEGG gene")
        return

    def _process_omim2disease(self, limit=None):
        """
        This method maps the KEGG disease IDs to the corresponding OMIM disease IDs.
        Currently this only maps KEGG diseases and OMIM diseases that have a 1:1 mapping.

        Triples created:
        <kegg_disease_id> is a class
        <omim_disease_id> is a class
        <kegg_disease_id> hasXref <omim_disease_id>
        :param limit:
        :return:
        """

        logger.info("Processing 1:1 KEGG disease to OMIM disease mappings")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['omim2disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (omim_disease_id, kegg_disease_id, link_type) = row

                kegg_disease_id = 'KEGG-'+kegg_disease_id.strip()
                omim_disease_id = re.sub('omim', 'OMIM', omim_disease_id)

                # Create hash for the links from OMIM ID -> KEGG ID
                if omim_disease_id not in self.omim_disease_hash:
                    self.omim_disease_hash[omim_disease_id] = [kegg_disease_id]
                else:
                    self.omim_disease_hash[omim_disease_id].append(kegg_disease_id)

                # Create hash for the links from KEGG ID -> OMIM ID
                if kegg_disease_id not in self.kegg_disease_hash:
                    self.kegg_disease_hash[kegg_disease_id] = [omim_disease_id]
                else:
                    self.kegg_disease_hash[kegg_disease_id].append(omim_disease_id)

        # Now process the disease hashes and only process 1:1 omim disease:KEGG disease entries.
        for omim_disease_id in self.omim_disease_hash:
            if self.testMode and omim_disease_id not in self.test_ids['disease']:
                continue

            if (not self.testMode) and (limit is not None and line_counter > limit):
                break
            line_counter += 1

            if len(self.omim_disease_hash[omim_disease_id]) == 1:
                kegg_disease_id = ''.join(self.omim_disease_hash.get(omim_disease_id))
                if len(self.kegg_disease_hash[kegg_disease_id]) == 1:
                    # add ids, and deal with the labels separately
                    gu.addClassToGraph(g, kegg_disease_id, None)
                    gu.addClassToGraph(g, omim_disease_id, None)
                    gu.addEquivalentClass(g, kegg_disease_id, omim_disease_id)  # safe?
                    # gu.addXref(g, kegg_disease_id, omim_disease_id)

        logger.info("Done with KEGG disease to OMIM disease mappings.")
        return

    def _process_genes_kegg2ncbi(self, limit=None):
        """
        This method maps the KEGG human gene IDs to the corresponding NCBI Gene IDs.

        Triples created:
        <kegg_gene_id> is a class
        <ncbi_gene_id> is a class
        <kegg_gene_id> equivalentClass <ncbi_gene_id>
        :param limit:
        :return:
        """

        logger.info("Processing KEGG gene IDs to NCBI gene IDs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0

        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['ncbi']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (kegg_gene_id, ncbi_gene_id, link_type) = row

                if self.testMode and kegg_gene_id not in self.test_ids['genes']:
                    continue

                # Adjust the NCBI gene ID prefix.
                ncbi_gene_id = re.sub('ncbi-geneid', 'NCBIGene', ncbi_gene_id)
                kegg_gene_id = 'KEGG-'+kegg_gene_id

                # Adding the KEGG gene ID to the graph here is redundant, unless there happens to be
                # additional gene IDs in this table not present in the genes table.
                gu.addClassToGraph(g, kegg_gene_id, None)
                gu.addClassToGraph(g, ncbi_gene_id, None)
                gu.addEquivalentClass(g, kegg_gene_id, ncbi_gene_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with KEGG gene IDs to NCBI gene IDs")
        return

    def _make_variant_locus_id(self, gene_id, disease_id):
        """
        we actually want the association between the gene and the disease to be via an alternate locus
        not the "wildtype" gene itself.
        so we make an anonymous alternate locus, and put that in the association.
        we also make the label for the anonymous class, and add it to the label hash

        :param gene_id:
        :param disease_id:
        :return:
        """
        alt_locus_id = '_'+gene_id+'-'+disease_id+'VL'
        if self.nobnodes:
            alt_locus_id = ':'+alt_locus_id
        alt_label = self.label_hash.get(gene_id)
        disease_label = self.label_hash.get(disease_id)
        if alt_label is not None and alt_label != '':
            alt_label = 'some variant of ' + str(alt_label)
            if disease_label is not None and disease_label != '':
                alt_label += ' that is associated with ' + str(disease_label)
        else:
            alt_label = None

        self.label_hash[alt_locus_id] = alt_label

        return alt_locus_id

    def getTestSuite(self):
        import unittest
        from tests.test_kegg import KEGGTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(KEGGTestCase)

        return test_suite