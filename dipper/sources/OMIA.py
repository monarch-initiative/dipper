import logging
import xml.etree.ElementTree as ET
import re
import gzip
import io
import shutil

from dipper.sources.Source import Source
from dipper.sources.OMIM import OMIM, filter_keep_phenotype_entry_ids
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.assoc.DispositionAssoc import DispositionAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.sources.NCBIGene import NCBIGene
from dipper.utils.DipperUtil import DipperUtil
from dipper.models.Model import Model

logger = logging.getLogger(__name__)


class OMIA(Source):
    """
    This is the parser for the
    [Online Mendelian Inheritance in Animals
    (OMIA)](http://www.http://omia.angis.org.au),
    from which we process inherited disorders, other (single-locus) traits,
    and genes in >200 animal species (other than human and mouse and rats).

    We generate the omia graph to include the following information:
    * genes
    * animal taxonomy, and breeds as instances of those taxa
        (breeds are akin to "strains" in other taxa)
    * animal diseases, along with species-specific subtypes of those diseases
    * publications (and their mapping to PMIDs, if available)
    * gene-to-phenotype associations (via an anonymous variant-locus
    * breed-to-phenotype associations

    We make links between OMIA and OMIM in two ways:
    1.  mappings between OMIA and OMIM are created as OMIA --> hasdbXref OMIM
    2.  mappings between a breed and OMIA disease are created
        to be a model for the mapped OMIM disease,
        IF AND ONLY IF it is a 1:1 mapping.
        there are some 1:many mappings,
        and these often happen if the OMIM item is a gene.

    Because many of these species are not covered in
    the PANTHER orthology datafiles, we also pull any orthology
    relationships from the gene_group files from NCBI.

    """

    files = {
        'data': {
            'file': 'omia.xml.gz',
            'url': 'http://omia.angis.org.au/dumps/omia.xml.gz'},
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'omia')

        self.dataset = Dataset(
            'omia', 'Online Mendelian Inheritance in Animals',
            'http://omia.angis.org.au', None, None,
            'http://sydney.edu.au/disclaimer.shtml')

        self.id_hash = {
            'article': {},
            'phene': {},
            'breed': {},
            'taxon': {},
            'gene': {}
        }
        self.label_hash = {}
        # used to store the omia to omim phene mappings
        self.omia_omim_map = {}
        # used to store the unique genes that have phenes
        # (for fetching orthology)
        self.annotated_genes = set()

        self.test_ids = {
            'disease': [
                'OMIA:001702', 'OMIA:001867', 'OMIA:000478', 'OMIA:000201',
                'OMIA:000810', 'OMIA:001400'],
            'gene': [
                492297, 434, 492296, 3430235, 200685834, 394659996, 200685845,
                28713538, 291822383],
            'taxon': [9691, 9685, 9606, 9615, 9913, 93934, 37029, 9627, 9825],
            # to be filled in during parsing of breed table
            # for lookup by breed-associations
            'breed': []
        }
        # to store a map of omia ids and any molecular info
        # to write a report for curation
        self.stored_omia_mol_gen = {}
        self.g = self.graph
        return

    def fetch(self, is_dl_forced=False):
        """
        :param is_dl_forced:
        :return:
        """
        self.get_files(is_dl_forced)

        ncbi = NCBIGene(self.graph_type, self.are_bnodes_skized)
        # ncbi.fetch()
        gene_group = ncbi.files['gene_group']
        self.fetch_from_url(
            gene_group['url'], '/'.join((ncbi.rawdir, gene_group['file'])),
            False)

        return

    def parse(self, limit=None):
        # names of tables to iterate - probably don't need all these:
        # Article_Breed, Article_Keyword, Article_Gene, Article_Keyword,
        # Article_People, Article_Phene, Articles, Breed, Breed_Phene,
        # Genes_gb, Group_Categories, Group_MPO, Inherit_Type, Keywords,
        # Landmark, Lida_Links, OMIA_Group, OMIA_author, Omim_Xref, People,
        # Phene, Phene_Gene, Publishers, Resources, Species_gb, Synonyms

        self.scrub()

        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        if self.testMode:
            self.g = self.testgraph
        else:
            self.g = self.graph

        # we do three passes through the file
        # first process species (two others reference this one)
        self.process_species(limit)

        # then, process the breeds, genes, articles, and other static stuff
        self.process_classes(limit)

        # next process the association data
        self.process_associations(limit)

        # process the vertebrate orthology for genes
        # that are annotated with phenotypes
        ncbi = NCBIGene(self.graph_type, self.are_bnodes_skized)
        ncbi.add_orthologs_by_gene_group(self.g, self.annotated_genes)

        logger.info("Done parsing.")

        self.write_molgen_report()

        return

    def scrub(self):
        """
        The XML file seems to have mixed-encoding;
        we scrub out the control characters
        from the file for processing.

        i.e.?i
        omia.xml:1555328.28: PCDATA invalid Char value 2
        <field name="journal">Bulletin et Memoires de la Societe Centrale de Medic

        :return:

        """

        logger.info(
            "Scrubbing out the nasty characters that break our parser.")

        myfile = '/'.join((self.rawdir, self.files['data']['file']))
        tmpfile = '/'.join((self.rawdir, self.files['data']['file']+'.tmp.gz'))
        t = gzip.open(tmpfile, 'wb')
        du = DipperUtil()
        with gzip.open(myfile, 'rb') as f:
            filereader = io.TextIOWrapper(f, newline="")
            for l in filereader:
                l = du.remove_control_characters(l) + '\n'
                t.write(l.encode('utf-8'))
        t.close()
        # TEC I do not like this at all. original data must be preserved as is.
        # also may be heavy handed as chars which do not break the parser
        # are stripped as well (i.e. tabs and newlines)
        
        # move the temp file
        logger.info("Replacing the original data with the scrubbed file.")
        shutil.move(tmpfile, myfile)
        return

    # ###################### XML LOOPING FUNCTIONS ##################

    def process_species(self, limit):
        """
        Loop through the xml file and process the species.
        We add elements to the graph, and store the
        id-to-label in the label_hash dict.
        :param limit:
        :return:
        """

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = gzip.open(myfile, 'rb')
        filereader = io.TextIOWrapper(f, newline="")

        filereader.readline()  # remove the xml declaration line

        for event, elem in ET.iterparse(filereader):
            # Species ids are == genbank species ids!
            self.process_xml_table(
                elem, 'Species_gb', self._process_species_table_row, limit)

        f.close()

        return

    def process_classes(self, limit):
        """
        Loop through the xml file and process the articles,
        breed, genes, phenes, and phenotype-grouping classes.
        We add elements to the graph,
        and store the id-to-label in the label_hash dict,
        along with the internal key-to-external id in the id_hash dict.
        The latter are referenced in the association processing functions.

        :param limit:
        :return:

        """

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = gzip.open(myfile, 'rb')
        filereader = io.TextIOWrapper(f, newline="")

        filereader.readline()  # remove the xml declaration line

        parser = ET.XMLParser(encoding='utf-8')

        for event, elem in ET.iterparse(filereader, parser=parser):
            self.process_xml_table(
                elem, 'Articles', self._process_article_row, limit)
            self.process_xml_table(
                elem, 'Breed', self._process_breed_row, limit)
            self.process_xml_table(
                elem, 'Genes_gb', self._process_gene_row, limit)
            self.process_xml_table(
                elem, 'OMIA_Group', self._process_omia_group_row, limit)
            self.process_xml_table(
                elem, 'Phene', self._process_phene_row, limit)
            self.process_xml_table(
                elem, 'Omim_Xref', self._process_omia_omim_map, limit)

        f.close()

        # post-process the omia-omim associations to filter out the genes
        # (keep only phenotypes/diseases)
        self.clean_up_omim_genes()

        return

    def process_associations(self, limit):
        """
        Loop through the xml file and process the article-breed, article-phene,
        breed-phene, phene-gene associations, and the external links to LIDA.

        :param limit:
        :return:

        """

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = gzip.open(myfile, 'rb')
        filereader = io.TextIOWrapper(f, newline="")

        filereader.readline()  # remove the xml declaration line

        for event, elem in ET.iterparse(filereader):
            self.process_xml_table(
                elem, 'Article_Breed', self._process_article_breed_row, limit)
            self.process_xml_table(
                elem, 'Article_Phene', self._process_article_phene_row, limit)
            self.process_xml_table(
                elem, 'Breed_Phene', self._process_breed_phene_row, limit)
            self.process_xml_table(
                elem, 'Lida_Links', self._process_lida_links_row, limit)
            self.process_xml_table(
                elem, 'Phene_Gene', self._process_phene_gene_row, limit)
            self.process_xml_table(
                elem, 'Group_MPO', self._process_group_mpo_row, limit)

        f.close()

        return

    # ############ INDIVIDUAL TABLE-LEVEL PROCESSING FUNCTIONS ################

    def _process_species_table_row(self, row):
        # gb_species_id, sci_name, com_name, added_by, date_modified
        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
        sci_name = row['sci_name']
        com_name = row['com_name']
        model = Model(self.g)
        if self.testMode and \
                (int(row['gb_species_id']) not in self.test_ids['taxon']):
            return

        model.addClassToGraph(tax_id)
        if com_name != '':
            model.addSynonym(tax_id, com_name)
            self.label_hash[tax_id] = com_name  # for lookup later
        else:
            self.label_hash[tax_id] = sci_name

        return

    def _process_breed_row(self, row):
        model = Model(self.g)
        # in test mode, keep all breeds of our test species
        if self.testMode and \
                (int(row['gb_species_id']) not in self.test_ids['taxon']):
            return

        # save the breed keys in the test_ids for later processing
        self.test_ids['breed'] += [int(row['breed_id'])]

        breed_id = self.make_breed_id(row['breed_id'])

        self.id_hash['breed'][row['breed_id']] = breed_id
        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
        breed_label = row['breed_name']
        species_label = self.label_hash.get(tax_id)
        if species_label is not None:
            breed_label = breed_label + ' ('+species_label+')'

        model.addIndividualToGraph(breed_id, breed_label, tax_id)
        self.label_hash[breed_id] = breed_label

        return

    def _process_phene_row(self, row):
        model = Model(self.g)
        phenotype_id = None
        sp_phene_label = row['phene_name']
        if sp_phene_label == '':
            sp_phene_label = None
        if 'omia_id' not in row:
            logger.info("omia_id not present for %s", row['phene_id'])
            omia_id = self._make_internal_id('phene', phenotype_id)
        else:
            omia_id = 'OMIA:'+str(row['omia_id'])

        if self.testMode and not\
                (int(row['gb_species_id']) in self.test_ids['taxon'] and
                 omia_id in self.test_ids['disease']):
            return
        # add to internal hash store for later lookup
        self.id_hash['phene'][row['phene_id']] = omia_id

        descr = row['summary']
        if descr == '':
            descr = None

        # omia label
        omia_label = self.label_hash.get(omia_id)

        # add the species-specific subclass (TODO please review this choice)
        gb_species_id = row['gb_species_id']

        if gb_species_id != '':
            sp_phene_id = '-'.join((omia_id, gb_species_id))
        else:
            logger.error(
                "No species supplied in species-specific phene table for %s",
                omia_id)
            return

        species_id = 'NCBITaxon:'+str(gb_species_id)
        # use this instead
        species_label = self.label_hash.get('NCBITaxon:'+gb_species_id)
        if sp_phene_label is None and \
                omia_label is not None and species_label is not None:
            sp_phene_label = ' '.join((omia_label, 'in', species_label))
        model.addClassToGraph(
            sp_phene_id, sp_phene_label, omia_id, descr)
        # add to internal hash store for later lookup
        self.id_hash['phene'][row['phene_id']] = sp_phene_id
        self.label_hash[sp_phene_id] = sp_phene_label
        # add each of the following descriptions,
        # if they are populated, with a tag at the end.
        for item in [
                'clin_feat', 'history', 'pathology', 'mol_gen', 'control']:
            if row[item] is not None and row[item] != '':
                model.addDescription(
                    sp_phene_id, row[item] + ' ['+item+']')
        # if row['symbol'] is not None:  # species-specific
        # CHECK ME - sometimes spaces or gene labels
        #     gu.addSynonym(g, sp_phene, row['symbol'])

        model.addOWLPropertyClassRestriction(
            sp_phene_id, model.object_properties['in_taxon'],
            species_id)

        # add inheritance as an association
        inheritance_id = self._map_inheritance_term_id(row['inherit'])
        if inheritance_id is not None:
            assoc = DispositionAssoc(self.g, self.name, sp_phene_id, inheritance_id)
            assoc.add_association_to_graph()

        if row['characterised'] == 'Yes':
            self.stored_omia_mol_gen[omia_id] = {
                'mol_gen': row['mol_gen'],
                'map_info': row['map_info'],
                'species': row['gb_species_id']}

        return

    def write_molgen_report(self):
        import csv
        logger.info("Writing G2P report for OMIA")
        f = '/'.join((self.outdir, 'omia_molgen_report.txt'))

        with open(f, 'w', newline='\n') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            # write header
            h = ['omia_id', 'molecular_description', 'mapping_info', 'species']
            writer.writerow(h)
            for phene in self.stored_omia_mol_gen:
                writer.writerow((str(phene),
                                 self.stored_omia_mol_gen[phene]['mol_gen'],
                                 self.stored_omia_mol_gen[phene]['map_info'],
                                 self.stored_omia_mol_gen[phene]['species']))

        logger.info(
            "Wrote %d potential G2P descriptions for curation to %s",
            len(self.stored_omia_mol_gen), f)

        return

    def _process_article_row(self, row):
        model = Model(self.g)
        # don't bother in test mode
        if self.testMode:
            return

        iarticle_id = self._make_internal_id('article', row['article_id'])
        self.id_hash['article'][row['article_id']] = iarticle_id
        rtype = None
        if row['journal'] != '':
            rtype = Reference.ref_types['journal_article']
        reference = Reference(self.g, iarticle_id, rtype)

        if row['title'] is not None:
            reference.setTitle(row['title'].strip())
        if row['year'] is not None:
            reference.setYear(row['year'])
        reference.addRefToGraph()

        if row['pubmed_id'] is not None:
            pmid = 'PMID:'+str(row['pubmed_id'])
            self.id_hash['article'][row['article_id']] = pmid
            model.addSameIndividual(iarticle_id, pmid)
            model.addComment(pmid, iarticle_id.replace("_:", ''))

        return

    def _process_omia_group_row(self, row):
        model = Model(self.g)
        omia_id = 'OMIA:'+row['omia_id']

        if self.testMode and omia_id not in self.test_ids['disease']:
            return

        group_name = row['group_name']
        group_summary = row['group_summary']

        disease_id = None
        group_category = row.get('group_category')
        disease_id = \
            self.map_omia_group_category_to_ontology_id(group_category)
        if disease_id is not None:
            model.addClassToGraph(disease_id, None)
            if disease_id == 'MP:0008762':  # embryonic lethal
                # add this as a phenotype association
                # add embryonic onset
                assoc = D2PAssoc(self.g, self.name, omia_id, disease_id)
                assoc.add_association_to_graph()
                disease_id = None
        else:
            logger.info(
                "No disease superclass defined for %s:  %s",
                omia_id, group_name)
            # default to general disease  FIXME this may not be desired
            disease_id = 'DOID:4'

        if group_summary == '':
            group_summary = None
        if group_name == '':
            group_name = None

        model.addClassToGraph(omia_id, group_name, disease_id, group_summary)

        self.label_hash[omia_id] = group_name

        return

    def _process_gene_row(self, row):
        model = Model(self.g)
        geno = Genotype(self.g)
        if self.testMode and row['gene_id'] not in self.test_ids['gene']:
            return
        gene_id = 'NCBIGene:'+str(row['gene_id'])
        self.id_hash['gene'][row['gene_id']] = gene_id
        gene_label = row['symbol']
        self.label_hash[gene_id] = gene_label
        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
        gene_type_id = NCBIGene.map_type_of_gene(row['gene_type'])
        model.addClassToGraph(gene_id, gene_label, gene_type_id)
        geno.addTaxon(tax_id, gene_id)

        return

    def _process_article_breed_row(self, row):
        model = Model(self.g)
        # article_id, breed_id, added_by
        # don't bother putting these into the test... too many!

        # and int(row['breed_id']) not in self.test_ids['breed']:
        if self.testMode:
            return

        article_id = self.id_hash['article'].get(row['article_id'])
        breed_id = self.id_hash['breed'].get(row['breed_id'])

        # there's some missing data (article=6038).  in that case skip
        if article_id is not None:
            self.g.addTriple(
                article_id, model.object_properties['is_about'], breed_id)
        else:
            logger.warning("Missing article key %s", str(row['article_id']))

        return

    def _process_article_phene_row(self, row):
        """
        Linking articles to species-specific phenes.

        :param row:
        :return:
        """
        model = Model(self.g)
        # article_id, phene_id, added_by
        # look up the article in the hashmap
        phenotype_id = self.id_hash['phene'].get(row['phene_id'])
        article_id = self.id_hash['article'].get(row['article_id'])

        omia_id = self._get_omia_id_from_phene_id(phenotype_id)
        if self.testMode and omia_id not in self.test_ids['disease'] \
                or phenotype_id is None or article_id is None:
            return

        # make a triple, where the article is about the phenotype
        self.g.addTriple(
            article_id,
            model.object_properties['is_about'], phenotype_id)

        return

    def _process_breed_phene_row(self, row):
        model = Model(self.g)
        # Linking disorders/characteristic to breeds
        # breed_id, phene_id, added_by
        breed_id = self.id_hash['breed'].get(row['breed_id'])
        phene_id = self.id_hash['phene'].get(row['phene_id'])

        # get the omia id
        omia_id = self._get_omia_id_from_phene_id(phene_id)

        if (self.testMode and not (
                omia_id in self.test_ids['disease'] and
                int(row['breed_id']) in self.test_ids['breed']) or
                breed_id is None or phene_id is None):
            return

        # FIXME we want a different relationship here
        assoc = G2PAssoc(
            self.g, self.name, breed_id, phene_id,
            model.object_properties['has_phenotype'])
        assoc.add_association_to_graph()

        # add that the breed is a model of the human disease
        # use the omia-omim mappings for this
        # we assume that we have already scrubbed out the genes
        # from the omim list, so we can make the model associations here

        omim_ids = self.omia_omim_map.get(omia_id)
        eco_id = "ECO:0000214"   # biological aspect of descendant evidence
        if omim_ids is not None and len(omim_ids) > 0:
            if len(omim_ids) > 1:
                logger.info(
                    "There's 1:many omia:omim mapping: %s, %s",
                    omia_id, str(omim_ids))
            for i in omim_ids:
                assoc = G2PAssoc(
                    self.g, self.name, breed_id, i,
                    model.object_properties['model_of'])
                assoc.add_evidence(eco_id)
                assoc.add_association_to_graph()
                aid = assoc.get_association_id()

                breed_label = self.label_hash.get(breed_id)
                if breed_label is None:
                    breed_label = "this breed"

                m = re.search(r'\((.*)\)', breed_label)
                if m:
                    sp_label = m.group(1)
                else:
                    sp_label = ''

                phene_label = self.label_hash.get(phene_id)
                if phene_label is None:
                    phene_label = "phenotype"
                elif phene_label.endswith(sp_label):
                    # some of the labels we made already include the species;
                    # remove it to make a cleaner desc
                    phene_label = re.sub(r' in '+sp_label, '', phene_label)
                desc = ' '.join(
                    ("High incidence of", phene_label, "in", breed_label,
                     "suggests it to be a model of disease", i + "."))
                model.addDescription(aid, desc)
        return

    def _process_lida_links_row(self, row):
        model = Model(self.g)
        # lidaurl, omia_id, added_by
        omia_id = 'OMIA:'+row['omia_id']
        lidaurl = row['lidaurl']

        if self.testMode and omia_id not in self.test_ids['disease']:
            return

        model.addXref(omia_id, lidaurl, True)

        return

    def _process_phene_gene_row(self, row):
        geno = Genotype(self.g)
        model = Model(self.g)
        gene_id = self.id_hash['gene'].get(row['gene_id'])
        phene_id = self.id_hash['phene'].get(row['phene_id'])

        omia_id = self._get_omia_id_from_phene_id(phene_id)

        if self.testMode and not (
                omia_id in self.test_ids['disease'] and
                row['gene_id'] in self.test_ids['gene']) or\
                gene_id is None or phene_id is None:
            return

        # occasionally some phenes are missing!  (ex: 406)
        if phene_id is None:
            logger.warning("Phene id %s is missing", str(row['phene_id']))
            return

        gene_label = self.label_hash[gene_id]
        # some variant of gene_id has phenotype d
        vl = '_:'+re.sub(r'NCBIGene:', '', str(gene_id)) + 'VL'
        geno.addAllele(vl, 'some variant of ' + gene_label)
        geno.addAlleleOfGene(vl, gene_id)
        geno.addAffectedLocus(vl, gene_id)
        model.addBlankNodeAnnotation(vl)
        assoc = G2PAssoc(self.g, self.name, vl, phene_id)
        assoc.add_association_to_graph()

        # add the gene id to the set of annotated genes
        # for later lookup by orthology
        self.annotated_genes.add(gene_id)

        return

    def _process_omia_omim_map(self, row):
        """
        Links OMIA groups to OMIM equivalents.
        :param row:
        :return:
        """
        # omia_id, omim_id, added_by
        model = Model(self.g)
        omia_id = 'OMIA:'+row['omia_id']
        omim_id = 'OMIM:'+row['omim_id']

        # also store this for use when we say that a given animal is
        # a model of a disease
        if omia_id not in self.omia_omim_map:
            self.omia_omim_map[omia_id] = set()
        self.omia_omim_map[omia_id].add(omim_id)

        if self.testMode and omia_id not in self.test_ids['disease']:
            return

        model.addXref(omia_id, omim_id)

        return

    def map_omia_group_category_to_ontology_id(self, category_num):
        """
        Using the category number in the OMIA_groups table,
        map them to a disease id.
        This may be superceeded by other MONDO methods.

        Platelet disorders will be more specific once
        https://github.com/obophenotype/human-disease-ontology/issues/46
        is fulfilled.

        :param category_num:
        :return:

        """

        category_map = {
            1: 'DOID:0014667',      # Inborn error of metabolism
            2: 'MESH:D004392',      # Dwarfism
            3: 'DOID:1682',         # congenital heart disease
            4: 'DOID:74',           # blood system disease
            5: 'DOID:3211',         # lysosomal storage disease
            6: 'DOID:16',           # integumentary system disease
            # --> retinal degeneration ==> OMIA:000830
            7: 'DOID:8466',         # progressive retinal atrophy
            8: 'DOID:0050572',      # Cone–rod dystrophy
            9: 'MESH:C536122',      # stationary night blindness
            10: 'Orphanet:98553',   # developmental retinal disorder
            11: 'DOID:5679',        # retinal disorder
            12: 'Orphanet:90771',   # Disorder of Sex Development
            #  - what to do about this one?
            13: 'MP:0008762',       # embryonic lethal
            # - not sure what to do with this
            14: None,               # blood group
            # FIXME make me more specific
            15: 'DOID:2218',        # intrinsic platelet disorder
            # FIXME make me more specific
            16: 'DOID:2218',        # extrinsic platelet disorder
            17: None  # transgenic ???
        }

        disease_id = None
        if category_num is not None and int(category_num) in category_map:
            disease_id = category_map.get(int(category_num))
            logger.info(
                "Found %s for category %s", str(disease_id), str(category_num))
        else:
            logger.info(
                "There's a group category I don't know anything about: %s",
                str(category_num))

        return disease_id

    def _process_group_mpo_row(self, row):
        """
        Make OMIA to MP associations
        :param row:
        :return:
        """
        omia_id = 'OMIA:'+row['omia_id']
        mpo_num = int(row['MPO_no'])
        mpo_id = 'MP:'+str(mpo_num).zfill(7)

        assoc = D2PAssoc(self.g, self.name, omia_id, mpo_id)
        assoc.add_association_to_graph()

        return

    def clean_up_omim_genes(self):
        omim = OMIM(self.graph_type, self.are_bnodes_skized)
        # get all the omim ids
        allomimids = set()
        for omia in self.omia_omim_map:
            allomimids.update(self.omia_omim_map[omia])

        entries_that_are_phenotypes = omim.process_entries(
            list(allomimids), filter_keep_phenotype_entry_ids, None, None)
        logger.info(
            "Filtered out %d/%d entries that are genes or features",
            len(allomimids)-len(entries_that_are_phenotypes), len(allomimids))

        # now iterate again and remove those non-phenotype ids
        removed_count = 0
        for omia in self.omia_omim_map:
            ids = self.omia_omim_map[omia]
            cleanids = set()
            for i in ids:
                if i in entries_that_are_phenotypes:
                    cleanids.add(i)
                else:
                    removed_count += 1  # keep track of how many we've removed
            self.omia_omim_map[omia] = cleanids

        logger.info(
            "Removed %d omim ids from the omia-to-omim map", removed_count)

        return

    def _make_internal_id(self, prefix, key):
        iid = '_:'+''.join(('omia', prefix, 'key', str(key)))
        return iid

    def make_breed_id(self, key):
        breed_id = 'OMIA-breed:'+str(key)

        return breed_id

    @staticmethod
    def _get_omia_id_from_phene_id(phene_id):
        omia_id = None
        if phene_id is not None:
            m = re.match(r'OMIA:\d+', str(phene_id))
            if m:
                omia_id = m.group(0)

        return omia_id

    @staticmethod
    def _map_inheritance_term_id(inheritance_symbol):

        inherit_map = {
            'A':  None,  # Autosomal
            'ACD': 'GENO:0000143',  # Autosomal co-dominant
            'ADV': None,  # autosomal dominant with variable expressivity
            'AID': 'GENO:0000259',  # autosomal incompletely dominant
            'ASD': 'GENO:0000145',  # autosomal semi-dominant
            # autosomal recessive, semi-lethal
            # using generic autosomal recessive
            'ASL': 'GENO:0000150',
            'D': 'GENO:0000147',  # autosomal dominant
            'M': None,  # multifactorial
            'MAT': None,  # Maternal
            # probably autosomal recessive
            # using generic autosomal recessive
            'PR':  'GENO:0000150',
            'R': 'GENO:0000150',  # Autosomal Recessive
            # Recessive Embryonic Lethal
            # using plain recessive
            'REL': 'GENO:0000148',
            # Autosomal Recessive Lethal
            # using plain autosomal recessive
            'RL': 'GENO:0000150',
            'S': 'GENO:0000146',  # Sex-linked   <--using allosomal dominant
            'SLi': None,  # Sex-limited
            'UD': 'GENO:0000144',  # Dominant
            'X': None,  # x-linked    # HP:0001417 ?
            # X-linked Dominant     <-- temp using allosomal dominant  FIXME
            'XLD': 'GENO:0000146',
            # X-linked Recessive    <-- temp using allosomal recessive  FIXME
            'XLR': 'GENO:0000149',
            'Y': None,  # Y-linked
            'Z': None,  # Z-linked
            # Z-linked recessive    <-- temp using allosomal recessive  FIXME
            'ZR': 'GENO:0000149',
            '999': None,  # Z-linked incompletely dominant
        }

        inheritance_id = inherit_map.get(inheritance_symbol)
        if inheritance_id is None and inheritance_symbol is not None:
            logger.warning(
                "No inheritance id is mapped for %s", inheritance_symbol)

        return inheritance_id

    def getTestSuite(self):
        import unittest
        from tests.test_omia import OMIATestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(OMIATestCase)

        return test_suite
