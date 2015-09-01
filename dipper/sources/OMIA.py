import logging
import xml.etree.ElementTree as ET

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.utils.GraphUtils import GraphUtils
from dipper.sources.NCBIGene import NCBIGene
from dipper import config
from dipper import curie_map


logger = logging.getLogger(__name__)


class OMIA(Source):
    """
    """

    files = {
        'data': {
            'file': 'omia.xml',
            'url': 'http://omia.angis.org.au/dumps/omia.xml.gz'}  # TODO change to gz
    }

    def __init__(self):
        Source.__init__(self, 'omia')

        self.load_bindings()

        self.dataset = Dataset('omia', 'Online Mendelian Inheritance in Animals', 'http://omia.angis.org.au',
                                None,
                                None,
                                'http://sydney.edu.au/disclaimer.shtml')

        # check to see if there's any ids configured in the config; otherwise, warn
        if 'test_ids' not in config.get_config() or 'disease' not in config.get_config()['test_ids']:
            logger.warn("not configured with disease test ids.")

        self.id_hash = {
            'article': {},
            'phene': {},
            'breed': {},
            'taxon': {},
            'gene': {}
        }
        self.label_hash = {}

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

        # we do two passes through the file
        # forst process species (two others reference this one)
        self.process_species(limit)

        # then, process the breeds, genes, articles, and other static stuff
        self.process_classes(limit)
        # next process the association data

        self._process_associations(limit)

        self.load_core_bindings()
        self.load_bindings()

        logger.info("Done parsing.")

        return

    def process_species(self, limit):
        """
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())

        # names of tables to iterate - probably don't need all these:
        # Article_Breed, Article_Keyword, Article_Gene, Article_Keyword, Article_People, Article_Phene,
        # Articles, Breed, Breed_Phene, Genes_gb, Group_Categories, Group_MPO, Inherit_Type, Keywords,
        # Landmark, Lida_Links, OMIA_Group, OMIA_author, Omim_Xref, People, Phene, Phene_Gene,
        # Publishers, Resources, Species_gb, Synonyms

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = open(myfile, 'r', errors='replace', encoding='utf-8')
        f.readline()  # remove the xml declaration line

        parser = ET.XMLParser(encoding='utf-8')
        # if i use utf-8 encoding, i get to line 1510906 - but there's some "control-B" chars to clean up

        for event, elem in ET.iterparse(myfile, parser=parser):
            if elem.tag == 'table_data':
                # get the element name and id
                # id = elem.get('id') # some internal identifier
                species = elem.find("[@name='Species_gb']")

                if species is not None:
                    # gb_species_id, sci_name, com_name, added_by, date_modified
                    # Gene info from Genbank
                    # Species ids are == genbank species ids!
                    logger.info("Processing Species")

                    row = {}
                    for r in species.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
                        sci_name = row['sci_name']
                        com_name = row['com_name']

                        gu.addClassToGraph(g, tax_id, sci_name)
                        if com_name != '':
                            gu.addSynonym(g, tax_id, com_name)
                            self.label_hash[tax_id] = com_name  # for lookup later
                        else:
                            self.label_hash[tax_id] = sci_name

                    ########################################
                    elem.clear()  # discard the element  (DO THIS LAST)

            if self.testMode and limit is not None and line_counter > limit:
                return

        print(self.label_hash)
        return

    def process_classes(self, limit):
        """
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())

        # names of tables to iterate - probably don't need all these:
        # Article_Breed, Article_Keyword, Article_Gene, Article_Keyword, Article_People, Article_Phene,
        # Articles, Breed, Breed_Phene, Genes_gb, Group_Categories, Group_MPO, Inherit_Type, Keywords,
        # Landmark, Lida_Links, OMIA_Group, OMIA_author, Omim_Xref, People, Phene, Phene_Gene,
        # Publishers, Resources, Species_gb, Synonyms

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = open(myfile, 'r', errors='replace', encoding='utf-8')
        f.readline()  # remove the xml declaration line

        parser = ET.XMLParser(encoding='utf-8')
        # if i use utf-8 encoding, i get to line 1510906 - but there's some "control-B" chars to clean up

        for event, elem in ET.iterparse(myfile, parser=parser):
            if elem.tag == 'table_data':
                # get the element name and id
                # id = elem.get('id') # some internal identifier
                articles = elem.find("[@name='Articles']")
                breed_data = elem.find("[@name='Breed']")
                genes = elem.find("[@name='Genes_gb']")
                omia_group = elem.find("[@name='OMIA_Group']")
                phene_data = elem.find("[@name='Phene']")


                if breed_data is not None:
                    logger.info("Processing Breeds")
                    row = {}
                    # breed_id, breed_name, gb_species_id, added_by, date_modified
                    for r in breed_data.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        breed_id = self._make_internal_id('breed', row['breed_id'])
                        self.id_hash['breed'][row['breed_id']] = breed_id
                        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
                        breed_label = row['breed_name']
                        species_label = self.label_hash.get(tax_id)
                        if species_label is not None:
                            breed_label = breed_label + ' ('+species_label+')'

                        gu.addIndividualToGraph(g, breed_id, breed_label, tax_id)   # should these be classes?

                elif phene_data is not None:
                    logger.info("Processing Phenes")
                    # These are the species-specific descriptions of a disease.
                    # Create the taxon-specific diseases as subclass restrictions on the disease id,
                    # which itself is a subclass of the main omia disease class
                    row = {}
                    # the following are the keys to this table:
                    # phene_id, omia_id, gb_species_id, phene_name, summary, symbol, marker, clin_feat,
                    # gen_test, inherit, inherit_text, mol_gen, map_info, history,
                    # control, pathology, prevalence, defect, singlelocus, characterised, added_by, date_modified

                    for r in phene_data.findall('row'):
                        for f in r.findall('field'):  # get the elements of the row
                            ats = f.attrib
                            row[ats['name']] = f.text

                        phenotype_id = omia_id = None
                        sp_phene_label = row['phene_name']
                        if sp_phene_label == '':
                            sp_phene_label = None
                        if 'omia_id' not in row:
                            logger.info("omia_id not present for %s", row['phene_id'])
                            omia_id = self._make_internal_id('phene', phenotype_id)
                        else:
                            omia_id = 'OMIA:'+str(row['omia_id'])
                        self.id_hash['phene'][row['phene_id']] = omia_id  # add to internal hash store for later lookup

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
                            logger.error("No species supplied in species-specific phene table for %s", omia_id)
                            continue

                        species_id = 'NCBITaxon:'+str(gb_species_id)
                        species_label = self.label_hash.get('NCBITaxon:'+gb_species_id)  # use this instead
                        if sp_phene_label is None and omia_label is not None and species_label is not None:
                            sp_phene_label = ' '.join((omia_label, 'in', species_label))
                        gu.addClassToGraph(g, sp_phene_id, sp_phene_label, omia_id, descr)
                        # add each of the following descriptions, if they are populated, with a tag at the end.
                        for item in ['clin_feat', 'history', 'pathology', 'mol_gen']:
                            if row[item] is not None and row[item] != '':
                                gu.addDescription(g, sp_phene_id, row[item] + ' ['+item+']')
                        # if row['symbol'] is not None:  # species-specific
                        #     gu.addSynonym(g, sp_phene, row['symbol'])  # CHECK ME - sometimes spaces or gene labels

                        # TODO add inheritance method!
                        gu.addOWLPropertyClassRestriction(g, sp_phene_id, gu.object_properties['in_taxon'], species_id)

                elif articles is not None:
                    logger.info("Processing Articles")
                    row = {
                        'article_id': None,
                        'title': None,
                        'journal': None,
                        'volume': None,
                        'pages': None,
                        'year': None,
                        'locus': None,
                        'abstract': None,
                        'publisher': None,
                        'pubmed_id': None,
                        'library': None,
                        'added_by': None,
                        'date_modified': None,
                    }
                    for r in articles.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        iarticle_id = self._make_internal_id('article', row['article_id'])
                        self.id_hash['article'][row['article_id']] = iarticle_id
                        rtype = None
                        if row['journal'] != '':
                            rtype = Reference.ref_types['journal_article']
                        r = Reference(iarticle_id, rtype)

                        if row['title'] is not None:
                            r.setTitle(row['title'])
                        if row['year'] is not None:
                            r.setYear(row['year'])
                        r.addRefToGraph(g)

                        if row['pubmed_id'] is not None:
                            pmid = 'PMID:'+str(row['pubmed_id'])
                            self.id_hash['article'][row['article_id']] = pmid
                            gu.addSameIndividual(g, iarticle_id, pmid)


                elif omia_group is not None:
                    logger.info("Processing OMIA groups")
                    # OMIA groups are the actual omia classes, and their general classification
                    # this the group name and summary are the class label and definition.
                    # omia_id, group_name, group_summary, group_category, added_by, date_modified

                    row = {}
                    for r in omia_group.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        omia_id = 'OMIA:'+row['omia_id']

                        group_name = row['group_name']
                        group_summary = row['group_summary']
                        group_category = row['group_category']
                        # TODO something to group category

                        if group_summary == '':
                            group_summary = None
                        if group_name == '':
                            group_name = None
                        gu.addClassToGraph(g, omia_id, group_name, None, group_summary)
                        self.label_hash[omia_id] = group_name

                elif genes is not None:
                    logger.info("Processing Genes")
                    # gene_id, gb_species_id, pubmed_id, symbol, gene_desc, gene_type, added_by, date_modified
                    # Gene info from Genbank
                    # Gene ids are == genbank gene ids!
                    row = {}
                    for r in genes.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text
                        gene_id = 'NCBIGene:'+str(row['gene_id'])
                        self.id_hash['gene'][row['gene_id']] = gene_id
                        gene_label =row['symbol']
                        self.label_hash[gene_id] = gene_label
                        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
                        gene_type_id = NCBIGene._map_type_of_gene(row['gene_type'])
                        gu.addClassToGraph(g, gene_id, gene_label, gene_type_id)
                        geno.addTaxon(tax_id, gene_id)


                    ########################################
                    elem.clear()  # discard the element  (DO THIS LAST)

            if self.testMode and limit is not None and line_counter > limit:
                return

        return

    def _process_associations(self, limit):
        """
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())

        # names of tables to iterate - probably don't need all these:
        # Article_Breed, Article_Keyword, Article_Gene, Article_Keyword, Article_People, Article_Phene,
        # Articles, Breed, Breed_Phene, Genes_gb, Group_Categories, Group_MPO, Inherit_Type, Keywords,
        # Landmark, Lida_Links, OMIA_Group, OMIA_author, Omim_Xref, People, Phene, Phene_Gene,
        # Publishers, Resources, Species_gb, Synonyms

        myfile = '/'.join((self.rawdir, self.files['data']['file']))

        f = open(myfile, 'r', errors='replace', encoding='utf-8')
        f.readline()  # remove the xml declaration line

        parser = ET.XMLParser(encoding='utf-8')
        # if i use utf-8 encoding, i get to line 1510906 - but there's some "control-B" chars to clean up

        for event, elem in ET.iterparse(myfile, parser=parser):
            if elem.tag == 'table_data':
                # get the element name and id
                # id = elem.get('id') # some internal identifier
                article_breed_data = elem.find("[@name='Article_Breed']")
                article_phene = elem.find("[@name='Article_Phene']")
                breed_phene = elem.find("[@name='Breed_Phene']")
                lida_links = elem.find("[@name='Lida_Links']")
                phene_gene = elem.find("[@name='Phene_Gene']")

                # process each section of the xml, each of which are tables.
                if article_breed_data is not None:
                    logger.info("Processing article-breed associations")
                    row = {}
                    # article_id, breed_id, added_by
                    for r in article_breed_data.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text
                            # make a mentions between the breed and the article

                        # TODO look up the keys to get the actual identifiers
                        article_id = self.id_hash['article'].get(row['article_id'])
                        breed_id = self._make_internal_id('breed', row['breed_id'])
                        # there's some missing data (article=6038).  in that case skip
                        if article_id is not None:
                            gu.addTriple(g, article_id, gu.object_properties['is_about'], breed_id)
                        else:
                            logger.warn("Missing article key %s", str(row['article_id']))


                elif article_phene is not None:
                    logger.info("Processing article-phene associations")
                    # The description for this is: Linking articles to species-specific phenes.
                    # article_id, phene_id, added_by
                    row = {}
                    for r in article_phene.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                    # look up the article in the hashmap
                    phenotype_id = self.id_hash['phene'][row['phene_id']]
                    article_id = self.id_hash['article'][row['article_id']]

                    # make a triple, where the article is about the phenotype
                    gu.addTriple(g, article_id, gu.object_properties['is_about'], phenotype_id)



                elif breed_phene is not None:
                    logger.info("Processing breed-phene associations")
                    # Linking disorders/characteristic to breeds
                    # breed_id, phene_id, added_by
                    row = {}
                    for r in breed_phene.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        breed_id = self.id_hash['breed'][row['breed_id']]
                        phene_id = self.id_hash['phene'][row['phene_id']]

                        # FIXME we want a different relationship here
                        assoc = G2PAssoc(self.name, breed_id, phene_id, gu.object_properties['is_marker_for'])
                        assoc.add_association_to_graph(g)


                elif lida_links is not None:
                    logger.info("Processing external links")
                    # External links to LIDA
                    # lidaurl, omia_id, added_by
                    row = {}
                    for r in lida_links.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        omia_id = 'OMIA:'+row['omia_id']
                        lidaurl = row['lidaurl']

                        gu.addXref(g, omia_id, lidaurl, True)

                elif phene_gene is not None:
                    logger.info("Processing phene-gene associations")
                    # gene_id, phene_id, added_by

                    row = {}
                    for r in phene_gene.findall('row'):
                        for f in r.findall('field'):
                            ats = f.attrib
                            row[ats['name']] = f.text

                        gene_id = self.id_hash['gene'][row['gene_id']]
                        phene_id = self.id_hash['phene'].get(row['phene_id'])

                        # occasionally some phenes are missing!  (ex: 406)
                        if phene_id is None:
                            logger.warn("Phene id %s is missing", str(row['phene_id']))
                            continue

                        gene_label = self.label_hash[gene_id]
                        # some variant of gene_id has phenotype d
                        vl = gene_id + 'VL'
                        geno.addAllele(vl, 'some variant of ' + gene_label)
                        assoc = G2PAssoc(self.name, vl, phene_id)
                        assoc.add_association_to_graph(g)



                    ########################################
                    elem.clear()  # discard the element  (DO THIS LAST)

            if self.testMode and limit is not None and line_counter > limit:
                return

        return

    def _make_internal_id(self, prefix, key):

        iid = '_'+''.join(('omia', prefix, 'key', str(key)))
        if self.nobnodes:
            iid = ':'+iid

        return iid


    # def getTestSuite(self):
    #     import unittest
    #     from tests.test_orphanet import OMIATestCase
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(OMIATestCase)
    #
    #     return test_suite
