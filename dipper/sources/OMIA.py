import logging
import xml.etree.ElementTree as ET
import re
import gzip
import io
import shutil
import csv

from dipper.sources.OMIMSource import OMIMSource
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.assoc.D2PAssoc import D2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.sources.NCBIGene import NCBIGene
from dipper.utils.DipperUtil import DipperUtil
from dipper.models.Model import Model

LOG = logging.getLogger(__name__)


class OMIA(OMIMSource):
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
        to be a 'is model of' the mapped OMIM disease,
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
            # CNAME broken? urllib not following redirects??
            # 'url': 'http://omia.angis.org.au/dumps/omia.xml.gz'
            'url': 'http://compldb.angis.org.au/dumps/omia.xml.gz',
            # see dipper/resources/omia/omia_xml.*  for xml xpaths and more
        },
        'causal_mutations':  {  # not used yet
            'file':  'causal_mutations.tab',
            'columns': [  # expected
                'gene_symbol',
                'ncbi_gene_id',
                'OMIA_id',
                'ncbi_tax_id',
                'OMIA_url',
                'phene_name'],
            'url': 'http://omia.org/curate/causal_mutations/?format=gene_table',
        },
    }

    def __init__(self, graph_type, are_bnodes_skolemized, data_release_version=None):

        super().__init__(
            graph_type=graph_type,
            are_bnodes_skolemized=are_bnodes_skolemized,
            name='omia',
            ingest_title='Online Mendelian Inheritance in Animals',
            ingest_url='https://omia.org',
            ingest_logo='source-omia.png',
            # ingest_desc=None,
            license_url=None,
            data_rights='http://sydney.edu.au/disclaimer.shtml',
            # file_handle=None
        )

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
                '492297', '434', '492296', '3430235', '200685834', '394659996',
                '200685845', '28713538', '291822383'],
            'taxon': [
                '9691', '9685', '9606', '9615', '9913', '93934', '37029', '9627',
                '9825'],
            # to be filled in during parsing of breed table
            # for lookup by breed-associations
            'breed': []
        }
        # to store a map of omia ids and any molecular info
        # to write a report for curation
        self.stored_omia_mol_gen = {}
        self.graph = self.graph

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
            gene_group['url'], '/'.join((ncbi.rawdir, gene_group['file'])), False)

    def parse(self, limit=None):
        # names of tables to iterate - probably don't need all these:
        # Article_Breed, Article_Keyword, Article_Gene, Article_Keyword,
        # Article_People, Article_Phene, Articles, Breed, Breed_Phene,
        # Genes_gb, Group_Categories, Group_MPO, Inherit_Type, Keywords,
        # Landmark, Lida_Links, OMIA_Group, OMIA_author, Omim_Xref, People,
        # Phene, Phene_Gene, Publishers, Resources, Species_gb, Synonyms

        self.scrub()

        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        if self.test_mode:
            self.graph = self.testgraph
        else:
            self.graph = self.graph

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
        ncbi.add_orthologs_by_gene_group(self.graph, self.annotated_genes)

        LOG.info("Done parsing.")

        self.write_molgen_report()

    def scrub(self):
        """
        The XML file seems to have mixed-encoding;
        we scrub out the control characters
        from the file for processing.

        i.e.?
        omia.xml:1555328.28: PCDATA invalid Char value 2
        <field name="journal">Bulletin et Memoires de la Societe Centrale de Medic

        :return:

        """

        LOG.info("Scrubbing out the nasty characters that break our parser.")

        myfile = '/'.join((self.rawdir, self.files['data']['file']))
        tmpfile = '/'.join((self.rawdir, self.files['data']['file']+'.tmp.gz'))
        tmp = gzip.open(tmpfile, 'wb')
        du = DipperUtil()
        with gzip.open(myfile, 'rb') as readbin:
            filereader = io.TextIOWrapper(readbin, newline="")
            for line in filereader:
                line = du.remove_control_characters(line) + '\n'
                tmp.write(line.encode('utf-8'))
        tmp.close()
        # TEC I do not like this at all. original data must be preserved as is.
        # also may be heavy handed as chars which do not break the parser
        # are stripped as well (i.e. tabs and newlines)
        # move the temp file
        LOG.info("Replacing the original data with the scrubbed file.")
        shutil.move(tmpfile, myfile)

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
        with gzip.open(myfile, 'rb') as readbin:
            filereader = io.TextIOWrapper(readbin, newline="")
            filereader.readline()  # remove the xml declaration line
            for event, elem in ET.iterparse(filereader):  # iterparse is not deprecated
                # Species ids are == NCBITaxon ids
                self.process_xml_table(
                    elem, 'Species_gb', self._process_species_table_row, limit)

    def process_classes(self, limit):
        """
        After all species have been processed .
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

        with gzip.open(myfile, 'rb') as readbin:
            filereader = io.TextIOWrapper(readbin, newline="")
            filereader.readline()  # remove the xml declaration line

            # iterparse is not deprecated
            for event, elem in ET.iterparse(filereader):
                self.process_xml_table(
                    elem, 'Articles', self._process_article_row, limit)
                self.process_xml_table(elem, 'Breed', self._process_breed_row, limit)
                self.process_xml_table(elem, 'Genes_gb', self._process_gene_row, limit)
                self.process_xml_table(
                    elem, 'OMIA_Group', self._process_omia_group_row, limit)
                self.process_xml_table(elem, 'Phene', self._process_phene_row, limit)
                self.process_xml_table(
                    elem, 'Omim_Xref', self._process_omia_omim_map, limit)

        # post-process the omia-omim associations to filter out the genes
        # (keep only phenotypes/diseases)
        self.clean_up_omim_genes()

    def process_associations(self, limit):
        """
        Loop through the xml file and process the article-breed, article-phene,
        breed-phene, phene-gene associations, and the external links to LIDA.

        :param limit:
        :return:

        """

        myfile = '/'.join((self.rawdir, self.files['data']['file']))
        with gzip.open(myfile, 'rb') as readbin:
            filereader = io.TextIOWrapper(readbin, newline="")
            filereader.readline()  # remove the xml declaration line
            for event, elem in ET.iterparse(filereader):  # iterparse is not deprecated
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

    # ############ INDIVIDUAL TABLE-LEVEL PROCESSING FUNCTIONS ################

    def _process_species_table_row(self, row):  # row is expected as a dict
        # gb_species_id, sci_name, com_name, added_by, date_modified
        tax_id = 'NCBITaxon:' + str(row['gb_species_id'])
        sci_name = row['sci_name']
        com_name = row['com_name']
        model = Model(self.graph)
        if self.test_mode and row['gb_species_id'] not in self.test_ids['taxon']:
            return

        model.addClassToGraph(tax_id)
        if com_name != '':
            model.addSynonym(tax_id, com_name)
            self.label_hash[tax_id] = com_name  # for lookup later
        else:
            self.label_hash[tax_id] = sci_name

    def _process_breed_row(self, row):
        model = Model(self.graph)
        # in test mode, keep all breeds of our test species
        if self.test_mode and row['gb_species_id'] not in self.test_ids['taxon']:
            return

        # save the breed keys in the test_ids for later processing
        self.test_ids['breed'] += [row['breed_id']]

        breed_id = 'OMIA-breed:' + str(row['breed_id'])

        self.id_hash['breed'][row['breed_id']] = breed_id
        tax_id = 'NCBITaxon:' + str(row['gb_species_id'])
        breed_label = row['breed_name']
        species_label = self.label_hash.get(tax_id)
        if species_label is not None:
            breed_label = breed_label + ' ('+species_label+')'

        model.addIndividualToGraph(breed_id, breed_label, tax_id)
        self.label_hash[breed_id] = breed_label

    def _process_phene_row(self, row):
        model = Model(self.graph)
        phenotype_id = None
        sp_phene_label = row['phene_name']
        if sp_phene_label == '':
            sp_phene_label = None
        if 'omia_id' not in row:
            LOG.info("omia_id not present for %s", row['phene_id'])
            omia_id = self._make_internal_id('phene', phenotype_id)
        else:
            omia_id = 'OMIA:' + str(row['omia_id'])

        if self.test_mode and not(  # demorgan this
                row['gb_species_id'] in self.test_ids['taxon'] and omia_id
                in self.test_ids['disease']):
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
            LOG.error(
                "No species supplied in species-specific phene table for %s", omia_id)
            return

        species_id = 'NCBITaxon:' + str(gb_species_id)
        # use this instead
        species_label = self.label_hash.get('NCBITaxon:'+gb_species_id)
        if sp_phene_label is None and omia_label is not None \
                and species_label is not None:
            sp_phene_label = ' '.join((omia_label, 'in', species_label))
        model.addClassToGraph(sp_phene_id, sp_phene_label, omia_id, descr)
        # add to internal hash store for later lookup
        self.id_hash['phene'][row['phene_id']] = sp_phene_id
        self.label_hash[sp_phene_id] = sp_phene_label
        # add each of the following descriptions,
        # if they are populated, with a tag at the end.
        for item in ['clin_feat', 'history', 'pathology', 'mol_gen', 'control']:
            if row[item] is not None and row[item] != '':
                model.addDescription(sp_phene_id, row[item] + ' ['+item+']')
        # if row['symbol'] is not None:  # species-specific
        # CHECK ME - sometimes spaces or gene labels
        #     gu.addSynonym(g, sp_phene, row['symbol'])

        model.addOWLPropertyClassRestriction(
            sp_phene_id, self.globaltt['in taxon'],
            species_id)

        # add inheritance as an association
        inheritance_id = None
        if row['inherit'] is not None and row['inherit'] in self.localtt:
            inheritance_id = self.resolve(row['inherit'])
        elif row['inherit'] is not None and row['inherit'] != '':
            LOG.info('Unhandled inheritance type:\t%s', row['inherit'])

        if inheritance_id is not None:  # observable related to genetic disposition
            assoc = D2PAssoc(
                self.graph, self.name, sp_phene_id, inheritance_id,
                rel=self.globaltt['has disposition'])
            assoc.add_association_to_graph()

        if row['characterised'] == 'Yes':
            self.stored_omia_mol_gen[omia_id] = {
                'mol_gen': row['mol_gen'],
                'map_info': row['map_info'],
                'species': row['gb_species_id']}

    def write_molgen_report(self):
        LOG.info("Writing G2P report for OMIA")
        filename = '/'.join((self.outdir, 'omia_molgen_report.txt'))

        with open(filename, 'w', newline='\n') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(  # write header
                ['omia_id', 'molecular_description', 'mapping_info', 'species'])
            for phene in self.stored_omia_mol_gen:
                writer.writerow((
                    str(phene),
                    self.stored_omia_mol_gen[phene]['mol_gen'],
                    self.stored_omia_mol_gen[phene]['map_info'],
                    self.stored_omia_mol_gen[phene]['species']))

        LOG.info(
            "Wrote %d potential G2P descriptions for curation to %s",
            len(self.stored_omia_mol_gen), filename)

    def _process_article_row(self, row):
        model = Model(self.graph)
        # don't bother in test mode
        if self.test_mode:
            return

        iarticle_id = self._make_internal_id('article', row['article_id'])
        self.id_hash['article'][row['article_id']] = iarticle_id
        rtype = None
        if row['journal'] != '':
            rtype = self.globaltt['journal article']
        reference = Reference(self.graph, iarticle_id, rtype)

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

    def _process_omia_group_row(self, row):
        model = Model(self.graph)
        omia_id = 'OMIA:' + row['omia_id']

        if self.test_mode and omia_id not in self.test_ids['disease']:
            return

        group_name = row['group_name']
        group_summary = row['group_summary']
        # default to general disease seems the only reasonable choice
        disease_id = self.globaltt['disease or disorder']
        group_category = 'group_category:' + str(row['group_category'])
        disease_id = self.resolve(group_category, False)

        if disease_id == 'group_category:None':
            disease_id = self.globaltt['disease']
        elif disease_id == group_category:
            LOG.info(
                "No disease superclass defined for %s:  %s  with parent %s",
                omia_id, group_name, group_category)
            disease_id = self.globaltt['disease']
        else:
            if disease_id == self.globaltt['embryonic lethality']:
                # add this as a phenotype association
                # add embryonic onset
                assoc = D2PAssoc(self.graph, self.name, omia_id, disease_id)
                assoc.add_association_to_graph()
                # disease_id = None
        model.addClassToGraph(disease_id, None)

        if group_summary == '':
            group_summary = None
        if group_name == '':
            group_name = None

        model.addClassToGraph(
            omia_id, group_name, description=group_summary, class_type=disease_id)

        self.label_hash[omia_id] = group_name

    def _process_gene_row(self, row):
        model = Model(self.graph)
        geno = Genotype(self.graph)
        if self.test_mode and row['gene_id'] not in self.test_ids['gene']:
            return
        gene_id = 'NCBIGene:' + str(row['gene_id'])
        self.id_hash['gene'][row['gene_id']] = gene_id
        gene_label = row['symbol']
        self.label_hash[gene_id] = gene_label
        tax_id = 'NCBITaxon:'+str(row['gb_species_id'])
        if row['gene_type'] is not None:
            gene_type_id = self.resolve(row['gene_type'])
            model.addClassToGraph(gene_id, gene_label, gene_type_id)
        geno.addTaxon(tax_id, gene_id)

    def _process_article_breed_row(self, row):

        # article_id, breed_id, added_by
        # don't bother putting these into the test... too many!

        # and row['breed_id'] not in self.test_ids['breed']:
        if self.test_mode:
            return

        article_id = self.id_hash['article'].get(row['article_id'])
        breed_id = self.id_hash['breed'].get(row['breed_id'])

        # there's some missing data (article=6038).  in that case skip
        if article_id is not None:
            self.graph.addTriple(article_id, self.globaltt['is_about'], breed_id)
        else:
            LOG.warning("Missing article key %s", str(row['article_id']))

    def _process_article_phene_row(self, row):
        """
        Linking articles to species-specific phenes.

        :param row:
        :return:
        """
        # article_id, phene_id, added_by
        # look up the article in the hashmap
        phenotype_id = self.id_hash['phene'].get(row['phene_id'])
        article_id = self.id_hash['article'].get(row['article_id'])

        omia_id = self._get_omia_id_from_phene_id(phenotype_id)
        if self.test_mode or omia_id not in self.test_ids['disease'] \
                or phenotype_id is None or article_id is None:
            return

        # make a triple, where the article is about the phenotype
        self.graph.addTriple(article_id, self.globaltt['is_about'], phenotype_id)

    def _process_breed_phene_row(self, row):
        model = Model(self.graph)
        # Linking disorders/characteristic to breeds
        # breed_id, phene_id, added_by
        breed_id = self.id_hash['breed'].get(row['breed_id'])
        phene_id = self.id_hash['phene'].get(row['phene_id'])

        # get the omia id
        omia_id = self._get_omia_id_from_phene_id(phene_id)

        if breed_id is None or phene_id is None or (
                self.test_mode and (
                    omia_id not in self.test_ids['disease'] or
                    row['breed_id'] not in self.test_ids['breed'])):
            return

        # FIXME we want a different relationship here
        assoc = G2PAssoc(
            self.graph, self.name, breed_id, phene_id, self.globaltt['has phenotype'])
        assoc.add_association_to_graph()

        # add that the breed is a model of the human disease
        # use the omia-omim mappings for this
        # we assume that we have already scrubbed out the genes
        # from the omim list, so we can make the model associations here

        omim_ids = self.omia_omim_map.get(omia_id)
        eco_id = self.globaltt['biological aspect of descendant evidence']
        if omim_ids is not None and omim_ids:
            # if len(omim_ids) > 1:
            #    LOG.info(
            #        "There's 1:many omia:omim mapping: %s, %s", omia_id, str(omim_ids))
            # else:
            #    oid = list(omim_ids)[0]
            #    LOG.info("OMIA %s is mapped to OMIM %s", omia_id, oid)

            for oid in omim_ids:
                assoc = G2PAssoc(
                    self.graph, self.name, breed_id, oid, self.globaltt['is model of'])
                assoc.add_evidence(eco_id)
                assoc.add_association_to_graph()
                aid = assoc.get_association_id()

                breed_label = self.label_hash.get(breed_id)
                if breed_label is None:  # get taxon label?
                    breed_label = "this breed"

                mch = re.search(r'\((.*)\)', breed_label)
                if mch:
                    sp_label = mch.group(1)
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
                     "suggests it to be a model of disease", oid + "."))
                model.addDescription(aid, desc)
        else:
            LOG.warning("No OMIM Disease associated with %s", omia_id)

    def _process_lida_links_row(self, row):
        model = Model(self.graph)
        # lidaurl, omia_id, added_by
        omia_id = 'OMIA:' + row['omia_id']
        lidaurl = row['lidaurl']

        if self.test_mode and omia_id not in self.test_ids['disease']:
            return

        model.addXref(omia_id, lidaurl, True)

    def _process_phene_gene_row(self, row):
        geno = Genotype(self.graph)
        model = Model(self.graph)
        gene_id = self.id_hash['gene'].get(row['gene_id'])
        phene_id = self.id_hash['phene'].get(row['phene_id'])

        omia_id = self._get_omia_id_from_phene_id(phene_id)

        if self.test_mode and not (
                omia_id in self.test_ids['disease'] and
                row['gene_id'] in self.test_ids['gene']
                ) or gene_id is None or phene_id is None:
            return

        # occasionally some phenes are missing!  (ex: 406)
        if phene_id is None:
            LOG.warning("Phene id %s is missing", str(row['phene_id']))
            return

        gene_label = self.label_hash[gene_id]
        # some variant of gene_id has phenotype d
        var = '_:' + gene_id.split(':')[-1] + 'VL'
        geno.addAllele(var, 'some variant of ' + gene_label)
        geno.addAlleleOfGene(var, gene_id)
        geno.addAffectedLocus(var, gene_id)
        model.addBlankNodeAnnotation(var)
        assoc = G2PAssoc(self.graph, self.name, var, phene_id)
        assoc.add_association_to_graph()

        # add the gene id to the set of annotated genes
        # for later lookup by orthology
        self.annotated_genes.add(gene_id)

    def _process_omia_omim_map(self, row):
        """
        Links OMIA groups to OMIM equivalents.
        :param row:
        :return:
        """
        # omia_id, omim_id, added_by
        model = Model(self.graph)
        omia_id = 'OMIA:' + row['omia_id']
        omim_id = 'OMIM:' + row['omim_id']

        # also store this for use when we say that a given animal is
        # a model of a disease
        if omia_id not in self.omia_omim_map:
            self.omia_omim_map[omia_id] = set()
        self.omia_omim_map[omia_id].add(omim_id)

        if self.test_mode and omia_id not in self.test_ids['disease']:
            return

        model.addXref(omia_id, omim_id)

    def _process_group_mpo_row(self, row):
        """
        Make OMIA to MP associations
        :param row:
        :return:
        """
        omia_id = 'OMIA:' + row['omia_id']
        mpo_num = row['MPO_no']
        mpo_id = 'MP:' + str(mpo_num).zfill(7)

        assoc = D2PAssoc(self.graph, self.name, omia_id, mpo_id)
        assoc.add_association_to_graph()

    def clean_up_omim_genes(self):
        '''
            Attempt to limit omim links to diseases and not genes/locus
        '''
        # get all the omim ids
        allomim_curie = set()
        for omia in self.omia_omim_map:
            allomim_curie.update(self.omia_omim_map[omia])
        # strip the curie prefix
        allomimids = set([o.split(':')[-1] for o in allomim_curie])

        LOG.info("Have %i omim_ids before filtering", len(allomimids))
        LOG.info("Exists %i omim_ids replaceable", len(self.omim_replaced))
        if self.omim_replaced:
            LOG.info(
                "Sample of each (all & replace) look like: %s , %s",
                list(allomimids)[0],
                list(self.omim_replaced.keys())[0])

        # deal with replaced identifiers
        replaced = allomimids & self.omim_replaced.keys()

        if replaced is not None and replaced:
            LOG.warning("These OMIM ID's are past their pull date: %s", str(replaced))
            for oid in replaced:
                allomimids.remove(oid)
                replacements = self.omim_replaced[oid]
                for rep in replacements:
                    allomimids.update(rep)

        # guard against omim identifiers which have been removed
        obsolete = [
            o for o in self.omim_type
            if self.omim_type[o] == self.globaltt['obsolete']]
        removed = allomimids & set(obsolete)
        if removed is not None and removed:
            LOG.warning("These OMIM ID's are gone: %s", str(removed))
            for oid in removed:
                allomimids.remove(oid)

        # get a list of omim ids which we consider to be for disease / phenotype
        omim_phenotypes = set([
            omim for omim in self.omim_type if self.omim_type[omim] in (
                self.globaltt['phenotype'],
                self.globaltt['has_affected_feature'],
                self.globaltt['heritable_phenotypic_marker'])])
        LOG.info(
            "Have %i omim_ids globally typed as phenotypes from OMIM",
            len(omim_phenotypes))

        entries_that_are_phenotypes = allomimids & omim_phenotypes

        LOG.info(
            "Filtered out %d/%d entries that are genes or features",
            len(allomimids - entries_that_are_phenotypes), len(allomimids))

        # now iterate again and remove those non-phenotype ids
        # this could be redone with set operations
        removed_count = 0
        for omia in self.omia_omim_map:
            cleanids = set()
            for dirty_curie in self.omia_omim_map[omia]:
                dirty_num = dirty_curie.split(':')[-1]
                if dirty_num in entries_that_are_phenotypes:
                    cleanids.add(dirty_curie)
                else:
                    removed_count += 1  # keep track of how many we've removed
            self.omia_omim_map[omia] = cleanids

        LOG.info("Removed %d omim ids from the omia-to-omim map", removed_count)

    @staticmethod
    def _make_internal_id(prefix, key):
        ''' more blank nodes '''
        return '_:' + ''.join(('omia', prefix, 'key', str(key)))

    @staticmethod
    def _get_omia_id_from_phene_id(phene_id):
        omia_id = None
        if phene_id is not None:
            mch = re.match(r'OMIA:\d+', str(phene_id))
            if mch:
                omia_id = mch.group(0)
        return omia_id

    def getTestSuite(self):
        import unittest
        from tests.test_omia import OMIATestCase
        test_suite = unittest.TestLoader().loadTestsFromTestCase(OMIATestCase)
        return test_suite
