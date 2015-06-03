import logging
import csv
import re
import unicodedata
import pysftp
from datetime import datetime
import stat
import os

from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Dataset import Dataset
from dipper import config
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from dipper.models.Genotype import Genotype
from dipper.models.G2PAssoc import G2PAssoc

logger = logging.getLogger(__name__)


class Coriell(Source):
    """
    The Coriell Catalog provided to Monarch includes metadata and descriptions of NIGMS, NINDS, NHGRI, and
    NIA cell lines.  These lines are made available for research purposes.
    Here, we create annotations for the cell lines as models of the diseases from which
    they originate.

    We create a handle for a patient from which the given cell line is derived (since there may be multiple
    cell lines created from a given patient).  A genotype is assembled for a patient, which includes
    a karyotype (if specified) or a collection of variants.  Both the genotype (has_genotype) and disease are linked
    to the patient (has_phenotype), and the cell line is listed as derived from the patient.
    The cell line is classified by it's [CLO cell type](http://www.ontobee.org/browser/index.php?o=clo),
    which itself is linked to a tissue of origin.

    Unfortunately, the omim numbers listed in this file are both for genes and diseases; we have no way of knowing
    a priori if a designated omim number is a gene or disease; so we presently link the patient to any omim id
    via the has_phenotype relationship.

    Notice: The Coriell catalog is delivered to Monarch in a specific format, and requires ssh rsa fingerprint
    identification.  Other groups wishing to get this data in it's raw form will need to contact Coriell
    for credentials.  This needs to be placed into your configuration file for it to work.
    """

    terms = {
        'cell_line_repository': 'CLO:0000008',
        'race': 'SIO:001015',
        'ethnic_group': 'EFO:0001799',
        'age': 'EFO:0000246',
        'sampling_time': 'EFO:0000689',
    }

    files = {
        'ninds': {'file': 'NINDS_latest.csv',
                  'id': 'NINDS',
                  'label': 'NINDS Human Genetics DNA and Cell line Repository',
                  'page': 'https://catalog.coriell.org/1/NINDS'},
        'nigms': {'file': 'NIGMS_latest.csv',
                  'id': 'NIGMS',
                  'label': 'NIGMS Human Genetic Cell Repository',
                  'page': 'https://catalog.coriell.org/1/NIGMS'},
        'nia': {'file': 'NIA_latest.csv',
                'id': 'NIA',
                'label': 'NIA Aging Cell Repository',
                'page': 'https://catalog.coriell.org/1/NIA'},
        'nhgri': {'file': 'NHGRI_latest.csv',
                  'id': 'NHGRI',
                  'label': 'NHGRI Sample Repository for Human Genetic Research',
                  'page': 'https://catalog.coriell.org/1/NHGRI'}
    }

    # the following will house the specific cell lines to use for test output
    test_lines = ['ND02380', 'ND02381', 'ND02383', 'ND02384', 'GM17897', 'GM17898', 'GM17896', 'GM17944', 'GM17945',
                  'ND00055', 'ND00094', 'ND00136', 'GM17940', 'GM17939', 'GM20567', 'AG02506', 'AG04407', 'AG07602',
                  'AG07601', 'GM19700', 'GM19701', 'GM19702', 'GM00324', 'GM00325', 'GM00142']

    def __init__(self):
        Source.__init__(self, 'coriell')

        self.load_bindings()

        self.dataset = Dataset('coriell', 'Coriell', 'http://ccr.coriell.org/', None)

        # data-source specific warnings (will be removed when issues are cleared)

        logger.warn('We assume that if a species is not provided, that it is a Human-derived cell line')
        logger.warn('We map all omim ids as a disease/phenotype entity, but should be fixed in the future')

        # check if config exists; if it doesn't, error out and let user know
        if 'dbauth' not in config.get_config() or 'coriell' not in config.get_config()['dbauth']:
            logger.error("not configured with FTP user/password.")

        return

    def fetch(self, is_dl_forced):
        """
        Here we connect to the coriell sftp server using private connection details.  They dump bi-weekly files
        with a timestamp in the filename.  For each catalog, we poll the remote site and pull the most-recently
        updated file, renaming it to our local *_latest.csv.

        Be sure to have pg user/password connection details in your conf.json file, like:
        dbauth : {
        "coriell" : {"user" : "<username>", "password" : "<password>", "host" : <host>, "private_key"=path/to/rsa_key}
        }

        :param is_dl_forced:
        :return:
        """
        host = config.get_config()['dbauth']['coriell']['host']
        user = config.get_config()['dbauth']['coriell']['user']
        passwd = config.get_config()['dbauth']['coriell']['password']
        key = config.get_config()['dbauth']['coriell']['private_key']

        with pysftp.Connection(host, username=user, password=passwd, private_key=key) as sftp:
            files_by_repo = {'NIGMS': [], 'NIA': [], 'NHGRI': [], 'NINDS': []}
            most_recent_files = {}
            for attr in sftp.listdir_attr():
                # for each catalog, get the most-recent filename
                m = re.match('(NIGMS|NIA|NHGRI|NINDS)', attr.filename)
                if m is not None and len(m.groups()) > 0:
                    files_by_repo[m.group(1)] += [attr]
            # sort each array in hash, and get the name and time of the most-recent file for each catalog
            for r in files_by_repo:
                files_by_repo[r].sort(key=lambda x: x.st_mtime, reverse=True)
                if len(files_by_repo[r]) < 1:
                    # pass
                    logger.error("There were no files to download for the %s catalog", r)
                most_recent_file = files_by_repo[r][0]
                most_recent_files[r] = most_recent_file
            for r in most_recent_files:
                f = most_recent_files[r]
                target_name = r+"_latest.csv"
                target_name = '/'.join((self.rawdir, target_name))
                # check if the local file is out of date, if so, download.  otherwise, skip.
                # we rename (for simplicity) the original file
                st = None
                if os.path.exists(target_name):
                    st = os.stat(target_name)
                if st is None or f.st_mtime > st[stat.ST_CTIME]:
                    if st is None:
                        logger.info("File does not exist locally; downloading...")
                    else:
                        logger.info("There's a new version of %s catalog available; downloading...", r)
                    sftp.get(f.filename)
                    logger.info("Fetched %s", f.filename)
                    # move to the "latest" name
                    os.rename(f.filename, target_name)
                    # in our file hash, set the filename
                else:
                    logger.info("File %s exists; using local copy", f.filename)
                t = datetime.utcfromtimestamp(f.st_mtime).strftime("%Y-%m-%d")
                logger.info("Local file date: %s", datetime.utcfromtimestamp(st[stat.ST_CTIME]))
                self.dataset.setFileAccessUrl(f.filename)
                filedate = datetime.utcfromtimestamp(f.st_mtime).strftime("%Y-%m-%d")
                self.dataset.setVersion(filedate)
        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        for f in ['ninds', 'nigms', 'nia', 'nhgri']:
            file = '/'.join((self.rawdir, self.files[f]['file']))
            self._process_repository(self.files[f]['id'], self.files[f]['label'], self.files[f]['page'])
            self._process_data(file, limit)

        logger.info("Finished parsing.")

        self.load_bindings()

        return

    def _process_data(self, raw, limit=None):
        """
        This function will process the data files from Coriell.  We make the assumption that any alleles
        listed are variants (alternates to w.t.)

        Triples: (examples)

            :NIGMSrepository a CLO_0000008 #repository
                label : NIGMS Human Genetic Cell Repository
                foaf:page https://catalog.coriell.org/0/sections/collections/NIGMS/?SsId=8

            line_id a CL_0000057,  #fibroblast line
                derives_from patient_id
                part_of :NIGMSrepository
                RO:model_of OMIM:disease_id

            patient id a foaf:person,
                label: "fibroblast from patient 12345 with disease X"
                member_of family_id  #what is the right thing here?
                SIO:race EFO:caucasian  #subclass of EFO:0001799
                in_taxon NCBITaxon:9606
                dc:description Literal(remark)
                RO:has_phenotype OMIM:disease_id
                GENO:has_genotype genotype_id

            family_id a owl:NamedIndividual
                foaf:page "https://catalog.coriell.org/0/Sections/BrowseCatalog/FamilyTypeSubDetail.aspx?PgId=402&fam=2104&coll=GM"

            genotype_id a intrinsic_genotype
                GENO:has_alternate_part allelic_variant_id
                #we don't necessarily know much about the genotype, other than the allelic variant.
                #also there's the sex here)

            pub_id mentions cell_line_id

        :param raw:
        :param limit:
        :return:
        """
        logger.info("Processing Data from %s", raw)
        gu = GraphUtils(curie_map.get())

        if self.testMode:      # set the graph to build
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        geno = Genotype(g)

        gu.loadProperties(g, geno.object_properties, gu.OBJPROP)
        gu.loadAllProperties(g)

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                if not row:
                    pass
                else:
                    line_counter += 1

                    (catalog_id, description, omim_number, sample_type, cell_line_available,
                     dna_in_stock, dna_ref, gender, age, race, ethnicity, affected, karyotype,
                     relprob, mutation, gene, family_id, collection, url, cat_remark, pubmed_ids,
                     family_member, variant_id, dbsnp_id, species) = row

                    # example:
                    # GM00003,HURLER SYNDROME,607014,Fibroblast,Yes,No,,Female,26 YR,Caucasian,,,,
                    # parent,,,39,NIGMS Human Genetic Cell Repository,
                    # http://ccr.coriell.org/Sections/Search/Sample_Detail.aspx?Ref=GM00003,
                    # 46;XX; clinically normal mother of a child with Hurler syndrome; proband not in Repository,,
                    # 2,,18343,Homo sapiens

                    if self.testMode and catalog_id not in self.test_lines:
                        continue    # skip rows not in our test lines, when in test mode

                    # #############    BUILD REQUIRED VARIABLES    #############

                    # Make the cell line ID
                    cell_line_id = 'Coriell:'+catalog_id.strip()

                    # Map the cell/sample type
                    cell_type = self._map_cell_type(sample_type)

                    # Make a cell line label
                    line_label = collection.partition(' ')[0]+'-'+catalog_id.strip()

                    # Map the repository/collection
                    repository = self._map_collection(collection)

                    # patients are uniquely identified by one of:
                    # dbsnp id (which is == an individual haplotype)
                    # family id + family member (if present) OR   #probands are usually family member zero
                    # cell line id
                    # since some patients have >1 cell line derived from them, we must make sure that the genotype
                    # is attached to the patient, and can be inferred to the cell line
                    # examples of repeated patients are:  famid=1159, member=1; fam=152,member=1

                    # Make the patient ID

                    if family_id != '':
                        patient_id = 'MONARCH:Coriell-'+family_id+'-'+family_member  # FIXME this won't resolve like this
                    else:
                        patient_id = cell_line_id   # otherwise, just default to the cell line as the patient id

                    # properties of the individual patients:  sex, family id, member/relproband, description
                    # descriptions are really long and ugly SCREAMING text, so need to clean up
                    # the control cases are so odd with this labeling scheme; but we'll deal with it as-is for now.
                    short_desc = (description.split(';')[0]).capitalize()
                    if affected == 'Yes':
                        affected = 'affected'
                    else:
                        affected = 'unaffected'
                    gender = gender.lower()
                    patient_label = ' '.join((affected, gender, relprob))
                    if relprob == 'proband':
                        patient_label = ' '.join((patient_label, 'with', short_desc))
                    else:
                        patient_label = ' '.join((patient_label, 'of proband with', short_desc))

                    # #############    BUILD THE CELL LINE    #############

                    # Adding the cell line as a typed individual.
                    gu.addIndividualToGraph(g, cell_line_id, line_label, cell_type)

                    # Cell line derives from patient
                    geno.addDerivesFrom(cell_line_id, patient_id)

                    # Cell line a member of repository
                    gu.addMember(g, repository, cell_line_id)

                    if cat_remark != '':
                        gu.addDescription(g, cell_line_id, cat_remark)

                    # Cell age_at_sampling
                    # TODO add in the age nodes once these are modeled properly in #78
                    # if (age != ''):
                        # this would give a BNode that is an instance of Age.  but i don't know how to connect
                        # the age node to the cell line?  we need to ask @mbrush.
                        # age_id = '_'+re.sub('\s+','_',age)
                        # gu.addIndividualToGraph(g,age_id,age,self.terms['age'])
                        # gu.addTriple(g,age_id,self.properties['has_measurement'],age,True)

                    # that the cell line is a model for a disease is dealt with below

                    # #############    BUILD THE PATIENT    #############

                    # Add the patient ID as an individual.
                    gu.addPerson(g, patient_id, patient_label)
                    # TODO map relationship to proband as a class (what ontology?)

                    # Add race of patient
                    # FIXME: Adjust for subcategories based on ethnicity field
                    # EDIT: There are 743 different entries for ethnicity... Too many to map?
                    # Perhaps add ethnicity as a literal in addition to the mapped race?
                    # Need to adjust the ethnicity txt (if using) to initial capitalization to remove ALLCAPS

                    # TODO race should go into the individual's background and abstracted out to the Genotype class
                    # punting for now.
                    # if race != '':
                    #    mapped_race = self._map_race(race)
                    #    if mapped_race is not None:
                    #        gu.addTriple(g,patient_id,self.terms['race'],mapped_race)
                    #        gu.addSubclass(g,self.terms['ethnic_group'],mapped_race)

                    # #############    BUILD THE FAMILY    #############

                    # Add triples for family_id, if present.
                    # a family is just a small population
                    if family_id != '':
                        family_comp_id = 'CoriellFamily:'+family_id

                        # Add the family ID as a named individual
                        gu.addIndividualToGraph(g, family_comp_id, None, geno.genoparts['population'])

                        # Add the patient as a member of the family
                        gu.addMemberOf(g, patient_id, family_comp_id)

                    # #############    BUILD THE GENOTYPE   #############

                    # the important things to pay attention to here are:
                    # karyotype (used when there isn't a typical mutation), esp for chr rearrangements
                    # mutation = protein-level mutation as a label, often from omim
                    # gene = gene symbol - TODO get id
                    # variant_id = omim variant ids (; delimited)
                    # dbsnp_id = snp individual ids, is == patient id == haplotype?

                    # note GM00633 is a good example of chromosomal variation - do we have enough to capture this?
                    # GM00325 has both abnormal karyotype and variation

                    # make an assumption that if the taxon is blank, that it is human!
                    if species is None or species == '':
                        species = 'Homo sapiens'
                    taxon = self._map_species(species)

                    # if there's a dbSNP id, this is actually the individual's genotype
                    genotype_id = None
                    genotype_label = None
                    if dbsnp_id != '':
                        genotype_id = 'dbSNPIndividual:'+dbsnp_id.strip()

                    omim_map = {}
                    gvc_id = None

                    # some of the karyotypes are encoded with terrible hidden codes. remove them here
                    # i've seen a <98> character
                    karyotype = self.remove_control_characters(karyotype)
                    # TODO break down the karyotype into parts and map into GENO. depends on #77

                    if gene != '':
                        vl = gene+'['+mutation+']'
                    if karyotype.strip() != '' and karyotype != '46;XX' and karyotype != '46;XY':
                        mutation = mutation.strip()
                        gvc_id = self.make_id(variant_id.strip()+karyotype)
                        if mutation.strip() != '':
                            gvc_label = ';'.join((vl, karyotype))
                        else:
                            gvc_label = karyotype
                    elif variant_id.strip() != '':
                        gvc_id = self.make_id(variant_id)
                        gvc_label = vl

                    if variant_id.strip() != '':
                        # split the variants and add them as part of the genotype
                        # we don't necessarily know their zygosity, just that they are part of the genotype
                        # variant ids are from OMIM, so prefix as such
                        # we assume that the variants will be defined in OMIM rather than here.
                        # TODO sort the variant_id list, if the omim prefix is the same, then assume it's the locus
                        # make a hashmap of the omim id to variant id list; then build the genotype
                        # hashmap is also useful for removing the "genes" from the list of phenotypes
                        omim_map = {}

                        locus_num = None
                        for v in variant_id.split(';'):
                            # handle omim-style and odd var ids like 610661.p.R401X
                            m = re.match('(\d+)\.+(.*)', v.strip())

                            if m is not None and len(m.groups()) == 2:
                                (locus_num, var_num) = m.groups()

                            if locus_num is not None and locus_num not in omim_map:
                                omim_map[locus_num] = [var_num]
                            else:
                                omim_map[locus_num] = omim_map[locus_num]+[var_num]

                        for o in omim_map:

                            allele1_id = 'OMIM:'+o+'.'+omim_map.get(o)[0]
                            gu.addIndividualToGraph(g, allele1_id, None, geno.genoparts['variant_locus'])
                            if len(omim_map.get(o)) > 1:
                                # build a vslc
                                vslc_id = self.make_id(o + ' '.join(omim_map.get(o)))
                                vslc_label = vl
                                allele2_id = 'OMIM:'+o+'.'+omim_map.get(o)[1]
                                gu.addIndividualToGraph(g, allele2_id, None, geno.genoparts['variant_locus'])

                                # making the assumption that the alleles are variants!
                                geno.addPartsToVSLC(vslc_id, allele1_id, allele2_id,
                                                    geno.object_properties['has_alternate_part'])
                                geno.addVSLCtoParent(vslc_id, gvc_id)
                                gu.addIndividualToGraph(g, vslc_id, vslc_label,
                                                        geno.genoparts['variant_single_locus_complement'])
                                # note, we do not add the gene id to the variation because we don't have it.
                                # this will have to come from OMIM or ClinVar.

                            else:
                                geno.addParts(allele1_id, gvc_id, geno.properties['has_alternate_part'])
                    #else:
                        #genotype_id = None
                        #if affected == 'affected':
                        #    logger.info('affected individual without recorded variants: %s, %s',cell_line_id,description)

                    if affected == 'unaffected':
                        # let's just say that this person is wildtype
                        gu.addType(g, patient_id, geno.genoparts['wildtype'])
                    elif genotype_id is None:
                        # make an anonymous genotype id
                        genotype_id = '_geno'+cell_line_id

                    # add the gvc
                    if gvc_id is not None:
                        gu.addIndividualToGraph(g, gvc_id, gvc_label, geno.genoparts['genomic_variation_complement'])

                        # add the gvc to the genotype
                        if genotype_id is not None:
                            if affected == 'unaffected':
                                rel = geno.object_properties['has_reference_part']
                            else:
                                rel = geno.object_properties['has_alternate_part']
                            geno.addParts(gvc_id, genotype_id, rel)
                        genotype_label = gvc_label

                    if genotype_id is not None:
                        geno.addGenotype(genotype_id, genotype_label, geno.genoparts['intrinsic_genotype'])
                        geno.addTaxon(taxon, genotype_id)
                        # add that the patient has the genotype
                        gu.addTriple(g, patient_id, geno.properties['has_genotype'], genotype_id)
                    else:
                        geno.addTaxon(taxon, patient_id)

                    # TODO: Add sex/gender  (as part of the karyotype?)

                    ##############    DEAL WITH THE DISEASES   #############

                    # we associate the disease to the patient
                    if affected == 'affected':
                        if omim_number != '':
                            for d in omim_number.split(';'):
                                if d is not None and d != '':
                                    # if the omim number is in omim_map, then it is a gene not a pheno
                                    if d not in omim_map:
                                        disease_id = 'OMIM:'+d.strip()
                                        gu.addClassToGraph(g, disease_id, None)   # assume the label is taken care of

                                        # add the association: the patient has the disease
                                        assoc_id = self.make_id(patient_id+disease_id)
                                        assoc = G2PAssoc(assoc_id, patient_id, disease_id, None, None)
                                        assoc.addAssociationToGraph(g)

                                        # also, this line is a model of this disease
                                        # TODO abstract out model into it's own association class?
                                        gu.addTriple(g, cell_line_id, gu.properties['model_of'], disease_id)
                                    else:
                                        logger.info('removing %s from disease list since it is a gene', d)

                                        # BEWARE
                                        # Coriell currently seems to sometimes annotate variants to phenotype ids rather
                                        # than gene ids!  YIKES
                                        # we are checking with them now to find out why.
                                        # ex: GM16467 is annotated to 278720.0006, but it should be 613208.0006
                                        # so this current solution is aggressive, and might mean missing data

                    # #############    ADD PUBLICATIONS   #############

                    if pubmed_ids != '':
                        for s in pubmed_ids.split(';'):
                            pubmed_id = 'PMID:'+s.strip()
                            gu.addTriple(g, pubmed_id, gu.properties['mentions'], cell_line_id)

                    if not self.testMode and (limit is not None and line_counter > limit):
                        break

            Assoc().loadAllProperties(g)

        return

    def _process_repository(self, id, label, page):
        """
        This function will process the data supplied internally about the repository from Coriell.

        Triples:
            Repository a CLO_0000008 #repository
            rdf:label Literal (label)
            foaf:page Literal (page)


        :param raw:
        :param limit:
        :return:
        """
        # #############    BUILD THE CELL LINE REPOSITORY    #############
        for g in [self.graph, self.testgraph]:
            # FIXME: How to devise a label for each repository?
            gu = GraphUtils(curie_map.get())
            repo_id = 'CoriellCollection:'+id
            repo_label = label
            repo_page = page
            #gu.addClassToGraph(g,repo_id,repo_label)
            gu.addIndividualToGraph(g, repo_id, repo_label, self.terms['cell_line_repository'])
            gu.addPage(g, repo_id, repo_page)

        return

    def _map_cell_type(self, sample_type):
        ctype = None
        type_map = {
            'Adipose stromal cell': 'CL:0002570',  # FIXME: mesenchymal stem cell of adipose
            'Amniotic fluid-derived cell line': 'CL:0002323',  # FIXME: amniocyte?
            'B-Lymphocyte': 'CL:0000236',  # B cell
            'Chorionic villus-derived cell line': 'CL:0000000',  # FIXME: No Match
            'Endothelial': 'CL:0000115',  # endothelial cell
            'Epithelial': 'CL:0000066',  # epithelial cell
            'Erythroleukemic cell line': 'CL:0000000',  # FIXME: No Match. "Abnormal precursor (virally transformed) of mouse erythrocytes that can be grown in culture and induced to differentiate by treatment with, for example, DMSO."
            'Fibroblast': 'CL:0000057',  # fibroblast
            'Keratinocyte': 'CL:0000312',  # keratinocyte
            'Melanocyte': 'CL:0000148',  # melanocyte
            'Mesothelial': 'CL:0000077',
            'Microcell hybrid': 'CL:0000000',  # FIXME: No Match
            'Myoblast': 'CL:0000056',  # myoblast
            'Smooth muscle': 'CL:0000192',  # smooth muscle cell
            'Stem cell': 'CL:0000034',  # stem cell
            'T-Lymphocyte': 'CL:0000084',  # T cell
            'Tumor-derived cell line': 'CL:0002198'  # FIXME: No Match. "Cells isolated from a mass of neoplastic cells, i.e., a growth formed by abnormal cellular proliferation." Oncocyte? CL:0002198
        }
        if sample_type.strip() in type_map:
            ctype = type_map.get(sample_type)
        else:
            logger.error("Cell type not mapped: %s", sample_type)

        return ctype

    def _map_race(self, race):
        rtype = None
        type_map = {
            'African American': 'EFO:0003150',
            #'American Indian': 'EFO',
            'Asian': 'EFO:0003152',
            'Asian; Other': 'EFO:0003152',  # FIXME: Asian?
            'Asiatic Indian': 'EFO:0003153',  # Asian Indian
            'Black': 'EFO:0003150',  # FIXME: African American? There is also African.
            'Caucasian': 'EFO:0003156',
            'Chinese': 'EFO:0003157',
            'East Indian': 'EFO:0003158',  # Eastern Indian
            'Filipino': 'EFO:0003160',
            'Hispanic/Latino': 'EFO:0003169',  # Hispanic: EFO:0003169, Latino: EFO:0003166
            'Japanese': 'EFO:0003164',
            'Korean': 'EFO:0003165',
            #'More than one race': 'EFO',
            #'Not Reported': 'EFO',
            #'Other': 'EFO',
            'Pacific Islander': 'EFO:0003154',  # Asian/Pacific Islander
            'Polynesian': 'EFO:0003154',  # Asian/Pacific Islander
            #'Unknown': 'EFO',
            'Vietnamese': 'EFO:0003152',  # Asian
        }
        if race.strip() in type_map:
            rtype = type_map.get(race)
        else:
            logger.warn("Race type not mapped: %s", race)

        return rtype

    def _map_species(self, species):
        tax = None
        type_map = {
            'Mus musculus': 'NCBITaxon:10090',
            'Peromyscus peromyscus californicus': 'NCBITaxon:42520',
            'Peromyscus peromyscus maniculatus': 'NCBITaxon:10042',
            'Peromyscus peromyscus leucopus': 'NCBITaxon:10041',
            'Peromyscus peromyscus polionotus': 'NCBITaxon:42413',
            'Macaca fascicularis': 'NCBITaxon:9541',
            'Rattus norvegicus': 'NCBITaxon:10116',
            'Papio anubis': 'NCBITaxon:9555',
            'Cricetulus griseus': 'NCBITaxon:10029',
            'Geochelone elephantopus': 'NCBITaxon:66189',
            'Muntiacus muntjak': 'NCBITaxon:9888',
            'Ailurus fulgens': 'NCBITaxon:9649',
            'Sus scrofa': 'NCBITaxon:9823',
            'Bos taurus': 'NCBITaxon:9913',
            'Oryctolagus cuniculus': 'NCBITaxon:9986',
            'Macaca nemestrina': 'NCBITaxon:9545',
            'Canis familiaris': 'NCBITaxon:9615',
            'Equus caballus': 'NCBITaxon:9796',
            'Macaca mulatta': 'NCBITaxon:9544',
            'Mesocricetus auratus': 'NCBITaxon:10036',
            'Macaca nigra': 'NCBITaxon:54600',
            'Erythrocebus patas': 'NCBITaxon:9538',
            'Pongo pygmaeus': 'NCBITaxon:9600',
            'Callicebus moloch': 'NCBITaxon:9523',
            'Lagothrix lagotricha': 'NCBITaxon:9519',
            'Saguinus fuscicollis': 'NCBITaxon:9487',
            'Saimiri sciureus': 'NCBITaxon:9521',
            'Saguinus labiatus': 'NCBITaxon:78454',
            'Pan paniscus': 'NCBITaxon:9597',
            'Ovis aries': 'NCBITaxon:9940',
            'Felis catus': 'NCBITaxon:9685',
            'Homo sapiens': 'NCBITaxon:9606'
        }
        if species.strip() in type_map:
            tax = type_map.get(species)
        else:
            logger.warn("Species type not mapped: %s", species)

        return tax

    def _map_collection(self, collection):
        ctype = None
        type_map = {
            'NINDS Repository': 'CoriellCollection:NINDS',
            'NIGMS Human Genetic Cell Repository': 'CoriellCollection:NIGMS',
            'NIA Aging Cell Culture Repository': 'CoriellCollection:NIA',
            'NHGRI Sample Repository for Human Genetic Research': 'CoriellCollection:NHGRI'
        }
        if collection.strip() in type_map:
            ctype = type_map.get(collection)
        else:
            logger.warn("ERROR: Collection type not mapped: %s", collection)

        return ctype

    def remove_control_characters(self, s):
        return "".join(ch for ch in s if unicodedata.category(ch)[0] != "C")

    def getTestSuite(self):
        import unittest
        from tests.test_coriell import CoriellTestCase
        #TODO add G2PAssoc, Genotype tests

        test_suite = unittest.TestLoader().loadTestsFromTestCase(CoriellTestCase)

        return test_suite
