import csv
import logging
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Family import Family
from dipper.models.Reference import Reference
from dipper.models.Pathway import Pathway
from dipper.models.Model import Model
from dipper import config
from dipper.utils.DipperUtil import DipperUtil


logger = logging.getLogger(__name__)


class KEGG(Source):
    files = {
        'disease': {
            'file': 'disease',
            'url': 'http://rest.genome.jp/list/disease'},
        'pathway': {
            'file': 'pathway',
            'url': 'http://rest.genome.jp/list/pathway'},
        'hsa_genes': {
            'file': 'hsa_genes',
            'url': 'http://rest.genome.jp/list/hsa'},
        'ortholog_classes': {
            'file': 'ortholog_classes',
            'url': 'http://rest.genome.jp/list/orthology'},
        'disease_gene': {
            'file': 'disease_gene',
            'url': 'http://rest.kegg.jp/link/disease/hsa'},
        'omim2disease': {
            'file': 'omim2disease',
            'url': 'http://rest.genome.jp/link/disease/omim'},
        'omim2gene': {
            'file': 'omim2gene',
            'url': 'http://rest.genome.jp/link/omim/hsa'},
        'ncbi': {
            'file': 'ncbi',
            'url': 'http://rest.genome.jp/conv/ncbi-geneid/hsa'},
        'hsa_gene2pathway': {
            'file': 'human_gene2pathway',
            'url': 'http://rest.kegg.jp/link/pathway/hsa'},
        'hsa_orthologs': {
            'file': 'hsa_orthologs',
            'url': 'http://rest.kegg.jp/link/orthology/hsa'},
        'mmu_orthologs': {
            'file': 'mmu_orthologs',
            'url': 'http://rest.kegg.jp/link/orthology/mmu'},
        'rno_orthologs': {
            'file': 'rno_orthologs',
            'url': 'http://rest.kegg.jp/link/orthology/rno'},
        'dme_orthologs': {
            'file': 'dme_orthologs',
            'url': 'http://rest.kegg.jp/link/orthology/dme'},
        'dre_orthologs': {
            'file': 'dre_orthologs',
            'url': 'http://rest.kegg.jp/link/orthology/dre'},
        'cel_orthologs': {
            'file': 'cel_orthologs',
            'url': 'http://rest.kegg.jp/link/orthology/cel'},
        'pathway_pubmed': {
            'file': 'pathway_pubmed',
            'url': 'http://rest.kegg.jp/link/pathway/pubmed'},
        'pathway_disease': {
            'file': 'pathway_disease',
            'url': 'http://rest.kegg.jp/link/pathway/ds'},
        # 'pathway_pathway': { # TEC 2016Mar03 does not exist here.
        #    'file': 'pathway_eq',
        #    'url': 'http://rest.kegg.jp/link/pathway/pathway'},
        'pathway_ko': {
            'file': 'pathway_ko',
            'url': 'http://rest.kegg.jp/link/pathway/ko'},
    }

    test_ids = {
        "pathway": [
            "path:map00010", "path:map00195", "path:map00100", "path:map00340",
            "path:hsa05223"],
        "disease": [
            "ds:H00015", "ds:H00026", "ds:H00712", "ds:H00736", "ds:H00014"],
        "genes": [
            "hsa:100506275", "hsa:285958", "hsa:286410", "hsa:6387",
            "hsa:1080", "hsa:11200", "hsa:1131", "hsa:1137", "hsa:126",
            "hsa:1277", "hsa:1278", "hsa:1285", "hsa:1548", "hsa:1636",
            "hsa:1639", "hsa:183", "hsa:185", "hsa:1910", "hsa:207",
            "hsa:2099", "hsa:2483", "hsa:2539", "hsa:2629", "hsa:2697",
            "hsa:3161", "hsa:3845", "hsa:4137", "hsa:4591", "hsa:472",
            "hsa:4744", "hsa:4835", "hsa:4929", "hsa:5002", "hsa:5080",
            "hsa:5245", "hsa:5290", "hsa:53630", "hsa:5630", "hsa:5663",
            "hsa:580", "hsa:5888", "hsa:5972", "hsa:6311", "hsa:64327",
            "hsa:6531", "hsa:6647", "hsa:672", "hsa:675", "hsa:6908",
            "hsa:7040", "hsa:7045", "hsa:7048", "hsa:7157", "hsa:7251",
            "hsa:7490", "hsa:7517", "hsa:79728", "hsa:83893", "hsa:83990",
            "hsa:841", "hsa:8438", "hsa:8493", "hsa:860", "hsa:9568",
            "hsa:9627", "hsa:9821", "hsa:999", "hsa:3460"],
        "orthology_classes": [
            "ko:K00010", "ko:K00027", "ko:K00042", "ko:K00088"]
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'kegg')

        # update the dataset object with details about this resource
        self.dataset = Dataset('kegg', 'KEGG', 'http://www.genome.jp/kegg/',
                               None, None,
                               'http://www.kegg.jp/kegg/legal.html')

        # check to see if there are any ids configured in the config;
        # otherwise, warn
        if 'test_ids' not in config.get_config() or\
                'disease' not in config.get_config()['test_ids']:
            logger.warning("not configured with disease test ids.")
        else:
            self.test_ids['disease'] += \
                config.get_config()['test_ids']['disease']

        self.label_hash = {}
        self.omim_disease_hash = {}  # to hold the mappings of omim:kegg ids
        self.kegg_disease_hash = {}  # to hold the mappings of kegg:omim ids

        return

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

        # TODO add versioning information from info rest call,
        # like http://rest.kegg.jp/info/pathway

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

        self._process_diseases(limit)
        self._process_genes(limit)
        self._process_genes_kegg2ncbi(limit)
        self._process_omim2gene(limit)
        self._process_omim2disease(limit)
        self._process_kegg_disease2gene(limit)

        self._process_pathways(limit)
        self._process_pathway_pubmed(limit)
        # self._process_pathway_pathway(limit)
        self._process_pathway_disease(limit)
        self._process_pathway_ko(limit)

        self._process_ortholog_classes(limit)
        # TODO add in when refactoring for #141
        # for f in ['hsa_orthologs', 'mmu_orthologs', 'rno_orthologs',
        #           'dme_orthologs','dre_orthologs','cel_orthologs']:
        #     file = '/'.join((self.rawdir, self.files[f]['file']))
        #     self._process_orthologs(file, limit)  # DONE #

        logger.info("Finished parsing")

        return

    def _process_pathways(self, limit=None):
        """
        This method adds the KEGG pathway IDs.
        These are the canonical pathways as defined in KEGG.
        We also encode the graphical depiction
        which maps 1:1 with the identifier.

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
        model = Model(g)
        line_counter = 0
        path = Pathway(g)
        raw = '/'.join((self.rawdir, self.files['pathway']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pathway_id, pathway_name) = row

                if self.testMode and \
                        pathway_id not in self.test_ids['pathway']:
                    continue

                pathway_id = 'KEGG-'+pathway_id.strip()
                path.addPathway(pathway_id, pathway_name)

                # we know that the pathway images from kegg map 1:1 here.
                # so add those
                image_filename = re.sub(r'KEGG-path:', '', pathway_id) + '.png'
                image_url = \
                    'http://www.genome.jp/kegg/pathway/map/'+image_filename
                model.addDepiction(pathway_id, image_url)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
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
        model = Model(g)
        raw = '/'.join((self.rawdir, self.files['disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (disease_id, disease_name) = row

                disease_id = 'KEGG-'+disease_id.strip()
                if disease_id not in self.label_hash:
                    self.label_hash[disease_id] = disease_name

                if self.testMode and\
                        disease_id not in self.test_ids['disease']:
                    continue

                # Add the disease as a class.
                # we don't get all of these from MONDO yet see:
                # https://github.com/monarch-initiative/human-disease-ontology/issues/3
                model.addClassToGraph(disease_id, disease_name)
                # not typing the diseases as DOID:4 yet because
                # I don't want to bulk up the graph unnecessarily

                if (not self.testMode) and (
                        limit is not None and line_counter > limit):
                    break

        logger.info("Done with diseases")
        return

    def _process_genes(self, limit=None):
        """
        This method processes the KEGG gene IDs.
        The label for the gene is pulled as
        the first symbol in the list of gene symbols;
        the rest are added as synonyms.
        The long-form of the gene name is added as a definition.
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
        model = Model(g)
        line_counter = 0
        family = Family(g)
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['hsa_genes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, gene_name) = row

                gene_id = 'KEGG-'+gene_id.strip()

                # the gene listing has a bunch of labels
                # that are delimited, as:
                # DST, BP240, BPA, BPAG1, CATX-15, CATX15, D6S1101, DMH, DT,
                # EBSB2, HSAN6, MACF2; dystonin; K10382 dystonin
                # it looks like the list is semicolon delimited
                # (symbol, name, gene_class)
                # where the symbol is a comma-delimited list

                # here, we split them up.
                # we will take the first abbreviation and make it the symbol
                # then take the rest as synonyms

                gene_stuff = re.split('r;', gene_name)
                symbollist = re.split(r',', gene_stuff[0])
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
                    model.addDefinition(gene_id, description)

                # add the rest of the symbols as synonyms
                for i in enumerate(symbollist, start=1):
                    model.addSynonym(gene_id, i[1].strip())

                if len(gene_stuff) > 2:
                    ko_part = gene_stuff[2]
                    ko_match = re.search(r'K\d+', ko_part)
                    if ko_match is not None and len(ko_match.groups()) == 1:
                        ko = 'KEGG-ko:'+ko_match.group(1)
                        family.addMemberOf(gene_id, ko)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        logger.info("Done with genes")
        return

    def _process_ortholog_classes(self, limit=None):
        """
        This method add the KEGG orthology classes to the graph.

        If there's an embedded enzyme commission number,
        that is added as an xref.

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
        model = Model(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, self.files['ortholog_classes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (orthology_class_id, orthology_class_name) = row

                if self.testMode and \
                        orthology_class_id not in \
                        self.test_ids['orthology_classes']:
                    continue

                # The orthology class is essentially a KEGG gene ID
                # that is species agnostic.
                # Add the ID and label as a gene family class

                other_labels = re.split(r'[;,]', orthology_class_name)
                # the first one is the label we'll use
                orthology_label = other_labels[0]

                orthology_class_id = 'KEGG-'+orthology_class_id.strip()

                orthology_type = OrthologyAssoc.terms['gene_family']
                model.addClassToGraph(orthology_class_id, orthology_label,
                                      orthology_type)
                if len(other_labels) > 1:
                    # add the rest as synonyms
                    # todo skip the first
                    for s in other_labels:
                        model.addSynonym(orthology_class_id, s.strip())

                    # add the last one as the description
                    d = other_labels[len(other_labels)-1]
                    model.addDescription(orthology_class_id, d)

                    # add the enzyme commission number (EC:1.2.99.5)as an xref
                    # sometimes there's two, like [EC:1.3.5.1 1.3.5.4]
                    # can also have a dash, like EC:1.10.3.-
                    ec_matches = re.findall(r'((?:\d+|\.|-){5,7})', d)
                    if ec_matches is not None:
                        for ecm in ec_matches:
                            model.addXref(orthology_class_id, 'EC:'+ecm)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
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
        model = Model(g)
        line_counter = 0
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, orthology_class_id) = row

                orthology_class_id = 'KEGG:'+orthology_class_id.strip()
                gene_id = 'KEGG:'+gene_id.strip()

                # note that the panther_id references a group of orthologs,
                # and is not 1:1 with the rest

                # add the KO id as a gene-family grouping class
                OrthologyAssoc(g, self.name, gene_id, None)\
                    .add_gene_family_to_graph(orthology_class_id)

                # add gene and orthology class to graph;
                # assume labels will be taken care of elsewhere
                model.addClassToGraph(gene_id, None)
                model.addClassToGraph(orthology_class_id, None)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        logger.info("Done with orthologs")
        return

    def _process_kegg_disease2gene(self, limit=None):
        """
        This method creates an association between diseases and
        their associated genes. We are being conservative here, and only
        processing those diseases for which there is no mapping to OMIM.

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
        model = Model(g)
        line_counter = 0
        geno = Genotype(g)
        rel = model.object_properties['is_marker_for']
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

                # only add diseases for which
                # there is no omim id and not a grouping class
                if disease_id not in self.kegg_disease_hash:
                    # add as a class
                    disease_label = None
                    if disease_id in self.label_hash:
                        disease_label = self.label_hash[disease_id]
                    if re.search(r'includ', str(disease_label)):
                        # they use 'including' when it's a grouping class
                        logger.info(
                            "Skipping this association because " +
                            "it's a grouping class: %s",
                            disease_label)
                        continue
                    # type this disease_id as a disease
                    model.addClassToGraph(disease_id, disease_label, 'DOID:4')
                    noomimset.add(disease_id)
                    alt_locus_id = self._make_variant_locus_id(gene_id,
                                                               disease_id)
                    alt_label = self.label_hash[alt_locus_id]
                    model.addIndividualToGraph(alt_locus_id, alt_label,
                                               geno.genoparts['variant_locus'])
                    geno.addAffectedLocus(alt_locus_id, gene_id)
                    model.addBlankNodeAnnotation(alt_locus_id)
                    # Add the disease to gene relationship.
                    assoc = G2PAssoc(g, self.name, alt_locus_id, disease_id, rel)
                    assoc.add_association_to_graph()

                if (not self.testMode) and (
                        limit is not None and line_counter > limit):
                    break

        logger.info("Done with KEGG disease to gene")
        logger.info("Found %d diseases with no omim id", len(noomimset))

        return

    def _process_omim2gene(self, limit=None):
        """
        This method maps the OMIM IDs and KEGG gene ID.
        Currently split based on the link_type field.
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
        model = Model(g)
        line_counter = 0
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['omim2gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (kegg_gene_id, omim_id, link_type) = row

                if self.testMode and \
                        kegg_gene_id not in self.test_ids['genes']:
                    continue

                kegg_gene_id = 'KEGG-'+kegg_gene_id.strip()
                omim_id = re.sub(r'omim', 'OMIM', omim_id)
                if link_type == 'equivalent':
                    # these are genes!
                    # so add them as a class then make equivalence
                    model.addClassToGraph(omim_id, None)
                    geno.addGene(kegg_gene_id, None)
                    if not DipperUtil.is_omim_disease(omim_id):
                        model.addEquivalentClass(kegg_gene_id, omim_id)
                elif link_type == 'reverse':
                    # make an association between an OMIM ID & the KEGG gene ID
                    # we do this with omim ids because
                    # they are more atomic than KEGG ids

                    alt_locus_id = self._make_variant_locus_id(kegg_gene_id,
                                                               omim_id)
                    alt_label = self.label_hash[alt_locus_id]
                    model.addIndividualToGraph(alt_locus_id, alt_label,
                                               geno.genoparts['variant_locus'])
                    geno.addAffectedLocus(alt_locus_id, kegg_gene_id)
                    model.addBlankNodeAnnotation(alt_locus_id)

                    # Add the disease to gene relationship.
                    rel = model.object_properties['is_marker_for']
                    assoc = G2PAssoc(g, self.name, alt_locus_id, omim_id, rel)
                    assoc.add_association_to_graph()

                elif link_type == 'original':
                    # these are sometimes a gene, and sometimes a disease
                    logger.info('Unable to handle original link for %s-%s',
                                kegg_gene_id, omim_id)
                else:
                    # don't know what these are
                    logger.warning('Unhandled link type for %s-%s: %s',
                                   kegg_gene_id, omim_id, link_type)

                if (not self.testMode) and (
                        limit is not None and line_counter > limit):
                    break

        logger.info("Done with OMIM to KEGG gene")

        return

    def _process_omim2disease(self, limit=None):
        """
        This method maps the KEGG disease IDs to
        the corresponding OMIM disease IDs.
        Currently this only maps KEGG diseases and OMIM diseases that are 1:1.

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
        model = Model(g)
        raw = '/'.join((self.rawdir, self.files['omim2disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (omim_disease_id, kegg_disease_id, link_type) = row

                kegg_disease_id = 'KEGG-'+kegg_disease_id.strip()
                omim_disease_id = re.sub(r'omim', 'OMIM', omim_disease_id)

                # Create hash for the links from OMIM ID -> KEGG ID
                if omim_disease_id not in self.omim_disease_hash:
                    self.omim_disease_hash[omim_disease_id] = [kegg_disease_id]
                else:
                    self.omim_disease_hash[
                        omim_disease_id].append(kegg_disease_id)

                # Create hash for the links from KEGG ID -> OMIM ID
                if kegg_disease_id not in self.kegg_disease_hash:
                    self.kegg_disease_hash[kegg_disease_id] = [omim_disease_id]
                else:
                    self.kegg_disease_hash[
                        kegg_disease_id].append(omim_disease_id)

        # Now process the disease hashes
        # and only pass 1:1 omim disease:KEGG disease entries.
        for omim_disease_id in self.omim_disease_hash:
            if self.testMode and \
                    omim_disease_id not in self.test_ids['disease']:
                continue

            if (not self.testMode) and (
                    limit is not None and line_counter > limit):
                break
            line_counter += 1

            if len(self.omim_disease_hash[omim_disease_id]) == 1:
                kegg_disease_id = \
                    ''.join(self.omim_disease_hash.get(omim_disease_id))
                if len(self.kegg_disease_hash[kegg_disease_id]) == 1:
                    # add ids, and deal with the labels separately
                    model.addClassToGraph(kegg_disease_id, None)
                    model.addClassToGraph(omim_disease_id, None)
                    # TODO is this safe?
                    model.addEquivalentClass(kegg_disease_id, omim_disease_id)
            else:
                pass
                # gu.addXref(g, omim_disease_id, kegg_disease_id)
                # TODO add xrefs if >1:1 mapping?

        logger.info("Done with KEGG disease to OMIM disease mappings.")
        return

    def _process_genes_kegg2ncbi(self, limit=None):
        """
        This method maps the KEGG human gene IDs
            to the corresponding NCBI Gene IDs.

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
        model = Model(g)
        line_counter = 0

        raw = '/'.join((self.rawdir, self.files['ncbi']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (kegg_gene_id, ncbi_gene_id, link_type) = row

                if self.testMode and \
                        kegg_gene_id not in self.test_ids['genes']:
                    continue

                # Adjust the NCBI gene ID prefix.
                ncbi_gene_id = re.sub(r'ncbi-geneid', 'NCBIGene', ncbi_gene_id)
                kegg_gene_id = 'KEGG-'+kegg_gene_id

                # Adding the KEGG gene ID to the graph here is redundant,
                # unless there happens to be additional gene IDs in this table
                # not present in the genes table.
                model.addClassToGraph(kegg_gene_id, None)
                model.addClassToGraph(ncbi_gene_id, None)
                model.addEquivalentClass(kegg_gene_id, ncbi_gene_id)

                if (not self.testMode) and (
                        limit is not None and line_counter > limit):
                    break

        logger.info("Done with KEGG gene IDs to NCBI gene IDs")
        return

    def _process_pathway_pubmed(self, limit):
        """
        Indicate that a pathway is annotated directly to a paper (is about)
            via it's pubmed id.
        :param limit:
        :return:
        """
        logger.info("Processing KEGG pathways to pubmed ids")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, self.files['pathway_pubmed']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pubmed_id, kegg_pathway_num) = row

                if self.testMode and \
                        kegg_pathway_num not in self.test_ids['pathway']:
                    continue

                pubmed_id = pubmed_id.upper()
                # will look like KEGG-path:map04130
                kegg_id = 'KEGG-'+kegg_pathway_num

                r = Reference(
                    g, pubmed_id, Reference.ref_types['journal_article'])
                r.addRefToGraph()
                g.addTriple(pubmed_id,
                            model.object_properties['is_about'], kegg_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_pathway_disease(self, limit):
        """
        We make a link between the pathway identifiers,
        and any diseases associated with them.
        Since we model diseases as processes, we make a triple saying that
        the pathway may be causally upstream of or within the disease process.

        :param limit:
        :return:

        """
        logger.info("Processing KEGG pathways to disease ids")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0

        model = Model(g)
        raw = '/'.join((self.rawdir, self.files['pathway_disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (disease_id, kegg_pathway_num) = row

                if self.testMode and \
                        kegg_pathway_num not in self.test_ids['pathway']:
                    continue

                disease_id = 'KEGG-'+disease_id
                # will look like KEGG-path:map04130 or KEGG-path:hsa04130
                pathway_id = 'KEGG-'+kegg_pathway_num

                g.addTriple(
                    pathway_id,
                    model.object_properties[
                        'causally_upstream_of_or_within'],
                    disease_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_pathway_pathway(self, limit):
        """
        There are "map" and "ko" identifiers for pathways.
        This makes equivalence mapping between them, where they exist.
        :param limit:
        :return:

        """
        logger.info("Processing KEGG pathways to other ids")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0

        model = Model(g)
        raw = '/'.join((self.rawdir, self.files['pathway_pathway']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pathway_id_1, pathway_id_2) = row

                if self.testMode and \
                        pathway_id_1 not in self.test_ids['pathway']:
                    continue

                pathway_id_1 = 'KEGG-'+pathway_id_1
                # will look like KEGG-path:map04130 or KEGG-path:ko04130
                pathway_id_2 = 'KEGG-'+pathway_id_2

                if pathway_id_1 != pathway_id_2:
                    model.addEquivalentClass(pathway_id_1, pathway_id_2)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_pathway_ko(self, limit):
        """
        This adds the kegg orthologous group (gene) to the canonical pathway.
        :param limit:

        :return:
        """
        logger.info("Processing KEGG pathways to kegg ortholog classes")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0

        raw = '/'.join((self.rawdir, self.files['pathway_ko']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (ko_id, pathway_id) = row

                if self.testMode and \
                        pathway_id not in self.test_ids['pathway']:
                    continue

                pathway_id = 'KEGG-'+pathway_id
                ko_id = 'KEGG-'+ko_id

                p = Pathway(g)
                p.addGeneToPathway(ko_id, pathway_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _make_variant_locus_id(self, gene_id, disease_id):
        """
        We actually want the association between the gene and the disease
        to be via an alternate locus not the "wildtype" gene itself.
        so we make an anonymous alternate locus,
        and put that in the association
        We also make the label for the anonymous class,
        and add it to the label hash

        :param gene_id:
        :param disease_id:
        :return:

        """
        alt_locus_id = '_:'+re.sub(r':', '', gene_id) +'-'+re.sub(r':', '', disease_id)+'VL'
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
