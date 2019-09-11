import logging
import csv
import re
from datetime import datetime
import stat
import os
import pysftp

from dipper.sources.Source import Source
from dipper import config
from dipper.models.Model import Model
from dipper.models.Genotype import Genotype
from dipper.models.Family import Family
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Reference import Reference
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.utils.DipperUtil import DipperUtil


LOG = logging.getLogger(__name__)


class Coriell(Source):
    """
    The Coriell Catalog provided to Monarch includes metadata and descriptions
    of NIGMS, NINDS, NHGRI, and NIA cell lines.  These lines are made available
    for research purposes. Here, we create annotations for the cell lines as
    models of the diseases from which they originate.

    We create a handle for a patient from which the given cell line is derived
    (since there may be multiple cell lines created from a given patient).
    A genotype is assembled for a patient, which includes a karyotype
    (if specified) and/or a collection of variants.
    Both the genotype (has_genotype) and disease are linked to the patient
    (has_phenotype), and the cell line is listed as derived from the patient.
    The cell line is classified by it's
    [CLO cell type](http://www.ontobee.org/browser/index.php?o=clo),
    which itself is linked to a tissue of origin.

    Unfortunately, the omim numbers listed in this file are both for genes
    & diseases; we have no way of knowing a priori if a designated omim number
    is a gene or disease; so we presently link the patient to any omim id via
    the has_phenotype relationship.

    Notice: The Coriell catalog is delivered to Monarch in a specific format,
    and requires ssh rsa fingerprint identification.  Other groups wishing to
    get this data in it's raw form will need to contact Coriell for credential
    This needs to be placed into your configuration file for it to work.

    """

    # head -1 NHGRI.csv | tr "[A-Z]," "[a-z]\n" | \
    #    awk  -v OFS='"' '{print "", $0, ",\t# " NR - 1}'
    column_labels = [
        "catalog_id",   # 0
        "description",  # 1
        "omim_num",     # 2
        "sample_type",  # 3
        "cell_line_available",  # 4
        "dna_instock",  # 5
        "dna_ref",      # 6
        "gender",       # 7
        "age",          # 8
        "race",         # 9
        "ethnicity",    # 10
        "affected",     # 11
        "karyotype",    # 12
        "relprob",      # 13
        "mutation",     # 14
        "gene",         # 15
        "fam",          # 16
        "collection",   # 17
        "url",          # 18
        "cat_remark",   # 19
        "pubmed_ids",   # 20
        "fammember",    # 21
        "variant_id",   # 22
        "dbsnp_id",     # 23
        "species",      # 24
    ]

    files = {
        'NINDS': {
            'file': 'NINDS.csv',
            'id': 'NINDS',
            'label': 'NINDS Human Genetics DNA and Cell line Repository',
            'page': 'https://catalog.coriell.org/1/NINDS',
            'columns': column_labels},
        'NIGMS': {
            'file': 'NIGMS.csv',
            'id': 'NIGMS',
            'label': 'NIGMS Human Genetic Cell Repository',
            'page': 'https://catalog.coriell.org/1/NIGMS',
            'columns': column_labels},
        'NIA': {
            'file': 'NIA.csv',
            'id': 'NIA',
            'label': 'NIA Aging Cell Repository',
            'page': 'https://catalog.coriell.org/1/NIA',
            'columns': column_labels},
        'NHGRI': {
            'file': 'NHGRI.csv',
            'id': 'NHGRI',
            'label': 'NHGRI Sample Repository for Human Genetic Research',
            'page': 'https://catalog.coriell.org/1/NHGRI',
            'columns': column_labels}
    }

    # the following will house the specific cell lines to use for test output
    test_lines = [
        'ND02380', 'ND02381', 'ND02383', 'ND02384', 'GM17897', 'GM17898',
        'GM17896', 'GM17944', 'GM17945', 'ND00055', 'ND00094', 'ND00136',
        'GM17940', 'GM17939', 'GM20567', 'AG02506', 'AG04407', 'AG07602'
        'AG07601', 'GM19700', 'GM19701', 'GM19702', 'GM00324', 'GM00325',
        'GM00142', 'NA17944', 'AG02505', 'GM01602', 'GM02455', 'AG00364',
        'GM13707', 'AG00780']

    def __init__(self, graph_type, are_bnodes_skolemized, skip_stats=False):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            name='coriell',
            ingest_title='Coriell Institute for Medical Research',
            ingest_url='https://ccr.coriell.org/',
            ingest_logo='https://github.com/monarch-initiative/monarch-ui/blob/master/public/img/sources/source-coriell.png',
            # website disclaimer 'https://www.coriell.org/1/About-Us/Legal-Notice'
            # wet material https://www.coriell.org/1/NINDS/About/Shared-Usage-Guidelines
            # license_url=None,
            # data_rights=None,
            # file_handle=None
        )

        # data-source specific warnings
        # (will be removed when issues are cleared)

        LOG.warning(
            'We assume that if a species is not provided, '
            'that it is a Human-derived cell line')
        LOG.warning(
            'We map all omim ids as a disease/phenotype entity, '
            'but should be fixed in the future')  # TODO

        # check if config exists; if it doesn't, error out and let user know
        if 'dbauth' not in config.get_config() or \
                'coriell' not in config.get_config()['dbauth']:
            LOG.error("not configured with FTP user/password.")
            # raise error

        return

    def fetch(self, is_dl_forced=False):
        """
        Here we connect to the coriell sftp server using private connection
        details.  They dump bi-weekly files with a timestamp in the filename.
        For each catalog, we ping the remote site and pull the most-recently
        updated file, renaming it to our local  latest.csv.

        Be sure to have pg user/password connection details in your conf.yaml
        file, like:
        dbauth : {"coriell" : {
        "user" : "<username>", "password" : "<password>",
        "host" : <host>, "private_key"=path/to/rsa_key}
        }

        :param is_dl_forced:
        :return:

        """

        host = config.get_config()['dbauth']['coriell']['host']
        key = config.get_config()['dbauth']['coriell']['private_key']
        user = config.get_config()['user']['coriell']
        passwd = config.get_config()['keys'][user]

        with pysftp.Connection(
                host, username=user, password=passwd, private_key=key) as sftp:
            # check to make sure each file is in there
            # get the remote files
            remote_files = sftp.listdir_attr()
            files_by_repo = {}
            for attr in remote_files:
                # for each catalog, get the most-recent filename
                mch = re.match('(NIGMS|NIA|NHGRI|NINDS)', attr.filename)
                if mch is not None and len(mch.groups()) > 0:
                    # there should just be one now
                    files_by_repo[mch.group(1)] = attr
            # sort each array in hash,
            # & get the name and time of the most-recent file for each catalog
            for rmt in self.files:
                LOG.info("Checking on %s catalog file", rmt)
                fname = self.files[rmt]['file']
                remotef = files_by_repo[rmt]
                target_name = '/'.join((self.rawdir, fname))
                # check if the local file is out of date, if so, download.
                # otherwise, skip.
                # we rename (for simplicity) the original file
                fstat = None
                if os.path.exists(target_name):
                    fstat = os.stat(target_name)
                    LOG.info(
                        "Local file date: %s",
                        datetime.utcfromtimestamp(fstat[stat.ST_CTIME]))
                if fstat is None or remotef.st_mtime > fstat[stat.ST_CTIME]:
                    if fstat is None:
                        LOG.info("File does not exist locally; downloading...")
                    else:
                        LOG.info(
                            "New version of %s catalog available; downloading...", rmt)
                    sftp.get(remotef.filename, target_name)
                    LOG.info(
                        "Fetched remote %s -> %s", remotef.filename, target_name)
                    fstat = os.stat(target_name)
                    filedate = datetime.utcfromtimestamp(
                        remotef.st_mtime).strftime("%Y-%m-%d")
                    LOG.info(
                        "New file date: %s",
                        datetime.utcfromtimestamp(fstat[stat.ST_CTIME]))
                else:
                    LOG.info("File %s exists; using local copy", fname)
                    filedate = datetime.utcfromtimestamp(
                        fstat[stat.ST_CTIME]).strftime("%Y-%m-%d")

                self.dataset.set_ingest_source(remotef.filename, True)
                self.dataset.set_ingest_source_file_version_date(remotef, filedate)

        return

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)
        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        for src_key in self.files:
            self._process_collection(
                self.files[src_key]['id'],
                self.files[src_key]['label'],
                self.files[src_key]['page'])
            self._process_data(src_key, limit)

        LOG.info("Finished parsing.")
        return

    def _process_data(self, src_key, limit=None):
        """
        This function will process the data files from Coriell.
        We make the assumption that any alleles listed are variants
        (alternates to w.t.)

        Triples: (examples)

        :NIGMSrepository a CLO_0000008 #repository
        label : NIGMS Human Genetic Cell Repository
        foaf:page
         https://catalog.coriell.org/0/sections/collections/NIGMS/?SsId=8

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
            foaf:page
             "https://catalog.coriell.org/0/Sections/BrowseCatalog/FamilyTypeSubDetail.aspx?PgId=402&fam=2104&coll=GM"

        genotype_id a intrinsic_genotype
            GENO:has_alternate_part allelic_variant_id
            we don't necessarily know much about the genotype,
            other than the allelic variant. also there's the sex here

        pub_id mentions cell_line_id

        :param raw:
        :param limit:
        :return:

        """

        raw = '/'.join((self.rawdir, self.files[src_key]['file']))

        LOG.info("Processing Data from %s", raw)

        if self.test_mode:      # set the graph to build
            graph = self.testgraph
        else:
            graph = self.graph

        family = Family(graph)
        model = Model(graph)

        line_counter = 1
        geno = Genotype(graph)
        diputil = DipperUtil()
        col = self.files[src_key]['columns']
        # affords access with
        # x = row[col.index('x')].strip()

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar=r'"')
            # we can keep a close watch on changing file formats
            fileheader = next(filereader, None)
            fileheader = [c.lower() for c in fileheader]
            if col != fileheader:  # assert
                LOG.error('Expected  %s to have columns: %s', raw, col)
                LOG.error('But Found %s to have columns: %s', raw, fileheader)
                raise AssertionError('Incomming data headers have changed.')

            for row in filereader:
                line_counter += 1
                if len(row) != len(col):
                    LOG.warning(
                        'Expected %i values but find %i in  row %i',
                        len(col), len(row), line_counter)
                    continue

                # (catalog_id, description, omim_number, sample_type,
                # cell_line_available, dna_in_stock, dna_ref, gender, age,
                # race, ethnicity, affected, karyotype, relprob, mutation,
                # gene, family_id, collection, url, cat_remark, pubmed_ids,
                # family_member, variant_id, dbsnp_id, species) = row

                # example:
                # GM00003,HURLER SYNDROME,607014,Fibroblast,Yes,No,
                #       ,Female,26 YR,Caucasian,,,,
                # parent,,,39,NIGMS Human Genetic Cell Repository,
                # http://ccr.coriell.org/Sections/Search/Sample_Detail.aspx?Ref=GM00003,
                # 46;XX; clinically normal mother of a child with Hurler syndrome;
                #       proband not in Repository,,
                # 2,,18343,Homo sapiens

                catalog_id = row[col.index('catalog_id')].strip()

                if self.test_mode and catalog_id not in self.test_lines:
                    # skip rows not in our test lines, when in test mode
                    continue

                # ###########    BUILD REQUIRED VARIABLES    ###########

                # Make the cell line ID
                cell_line_id = 'Coriell:' + catalog_id
                # Map the cell/sample type
                cell_type = self.resolve(row[col.index('sample_type')].strip())
                # on fail cell_type = self.globaltt['cell'] ?

                # Make a cell line label
                collection = row[col.index('collection')].strip()
                line_label = collection.partition(' ')[0] + '-' + catalog_id

                # Map the repository/collection
                repository = self.localtt[collection]

                # patients are uniquely identified by one of:
                # dbsnp id (which is == an individual haplotype)
                # family id + family member (if present) OR
                # probands are usually family member zero
                # cell line id
                # since some patients have >1 cell line derived from them,
                # we must make sure that the genotype is attached to
                # the patient, and can be inferred to the cell line
                # examples of repeated patients are:
                #   famid=1159, member=1; fam=152,member=1

                # Make the patient ID

                # make an anonymous patient
                patient_id = '_:person'
                fam_id = row[col.index('fam')].strip()
                fammember = row[col.index('fammember')].strip()
                if fam_id != '':
                    patient_id = '-'.join((patient_id, fam_id, fammember))
                else:
                    # make an anonymous patient
                    patient_id = '-'.join((patient_id, catalog_id))

                # properties of the individual patients:  sex, family id,
                # member/relproband, description descriptions are
                # really long and ugly SCREAMING text, so need to clean up
                # the control cases are so odd with this labeling scheme;
                # but we'll deal with it as-is for now.
                description = row[col.index('description')].strip()
                short_desc = (description.split(';')[0]).capitalize()

                gender = row[col.index('gender')].strip().lower()
                affected = row[col.index('affected')].strip()
                relprob = row[col.index('relprob')].strip()

                if affected == '':
                    affected = 'unspecified'
                elif affected in self.localtt:
                    affected = self.localtt[affected]
                else:
                    LOG.warning(
                        'Novel Affected status  %s at row: %i of %s',
                        affected, line_counter, raw)
                patient_label = ' '.join((affected, gender, relprob))
                if relprob == 'proband':
                    patient_label = ' '.join((
                        patient_label.strip(), 'with', short_desc))
                else:
                    patient_label = ' '.join((
                        patient_label.strip(), 'of proband with', short_desc))

                # #############    BUILD THE CELL LINE    #############

                # Adding the cell line as a typed individual.
                cell_line_reagent_id = self.globaltt['cell line']

                model.addIndividualToGraph(
                    cell_line_id, line_label, cell_line_reagent_id)

                # add the equivalent id == dna_ref
                dna_ref = row[col.index('dna_ref')].strip()
                if dna_ref != '' and dna_ref != catalog_id:
                    equiv_cell_line = 'Coriell:' + dna_ref
                    # some of the equivalent ids are not defined
                    # in the source data; so add them
                    model.addIndividualToGraph(
                        equiv_cell_line, None, cell_line_reagent_id)
                    model.addSameIndividual(cell_line_id, equiv_cell_line)

                # Cell line derives from patient
                geno.addDerivesFrom(cell_line_id, patient_id)
                geno.addDerivesFrom(cell_line_id, cell_type)

                # Cell line a member of repository
                family.addMember(repository, cell_line_id)

                cat_remark = row[col.index('cat_remark')].strip()

                if cat_remark != '':
                    model.addDescription(cell_line_id, cat_remark)

                # Cell age_at_sampling
                # TODO add the age nodes when modeled properly in #78
                # if (age != ''):
                    # this would give a BNode that is an instance of Age.
                    # but i don't know how to connect
                    # the age node to the cell line? we need to ask @mbrush
                    # age_id = '_'+re.sub('\s+','_',age)
                    # gu.addIndividualToGraph(
                    #   graph,age_id,age,self.globaltt['age'])
                    # gu.addTriple(
                    #   graph,age_id,self.globaltt['has measurement value'],age,
                    #   True)

                # #############    BUILD THE PATIENT    #############

                # Add the patient ID as an individual.
                model.addPerson(patient_id, patient_label)
                # TODO map relationship to proband as a class
                # (what ontology?)

                # Add race of patient
                # FIXME: Adjust for subcategories based on ethnicity field
                # EDIT: There are 743 different entries for ethnicity...
                # Too many to map?
                # Add ethnicity as literal in addition to the mapped race?
                # Adjust the ethnicity txt (if using)
                # to initial capitalization to remove ALLCAPS

                # TODO race should go into the individual's background
                # and abstracted out to the Genotype class punting for now.
                # if race != '':
                #    mapped_race = self.resolve(race)
                #    if mapped_race is not None:
                #        gu.addTriple(
                #           g,patient_id,self.globaltt['race'], mapped_race)
                #        model.addSubClass(
                #           mapped_race,self.globaltt['ethnic_group'])

                # #############    BUILD THE FAMILY    #############

                # Add triples for family_id, if present.
                if fam_id != '':
                    family_comp_id = 'CoriellFamily:' + fam_id

                    family_label = ' '.join(('Family of proband with', short_desc))

                    # Add the family ID as a named individual
                    model.addIndividualToGraph(
                        family_comp_id, family_label, self.globaltt['family'])

                    # Add the patient as a member of the family
                    family.addMemberOf(patient_id, family_comp_id)

                # #############    BUILD THE GENOTYPE   #############

                # the important things to pay attention to here are:
                # karyotype = chr rearrangements  (somatic?)
                # mutation = protein-level mutation as a label,
                # often from omim
                # gene = gene symbol - TODO get id
                # variant_id = omim variant ids (; delimited)
                # dbsnp_id = snp individual ids = full genotype?

                # note GM00633 is a good example of chromosomal variation
                # - do we have enough to capture this?
                # GM00325 has both abnormal karyotype and variation

                # make an assumption that if the taxon is blank,
                # that it is human!
                species = row[col.index('species')].strip()
                if species is None or species == '':
                    species = 'Homo sapiens'
                taxon = self.resolve(species)

                # if there's a dbSNP id,
                # this is actually the individual's genotype
                genotype_id = None
                genotype_label = None

                dbsnp_id = row[col.index('dbsnp_id')].strip()
                if dbsnp_id != '':
                    genotype_id = 'dbSNPIndividual:' + dbsnp_id

                omim_map = {}
                gvc_id = None

                # some of the karyotypes are encoded
                # with terrible hidden codes. remove them here
                # i've seen a <98> character
                karyotype = row[col.index('karyotype')].strip()
                karyotype = diputil.remove_control_characters(karyotype)
                karyotype_id = None
                if karyotype.strip() != '':
                    karyotype_id = '_:'+re.sub(
                        'MONARCH:', '', self.make_id(karyotype))
                    # add karyotype as karyotype_variation_complement
                    model.addIndividualToGraph(
                        karyotype_id, karyotype,
                        self.globaltt['karyotype_variation_complement'])
                    # TODO break down the karyotype into parts
                    # and map into GENO. depends on #77

                    # place the karyotype in a location(s).
                    karyo_chrs = self._get_affected_chromosomes_from_karyotype(
                        karyotype)
                    for chrom in karyo_chrs:
                        chr_id = makeChromID(chrom, taxon, 'CHR')
                        # add an anonymous sequence feature,
                        # each located on chr
                        karyotype_feature_id = '-'.join((karyotype_id, chrom))
                        karyotype_feature_label = \
                            'some karyotype alteration on chr' + str(chrom)
                        feat = Feature(
                            graph, karyotype_feature_id, karyotype_feature_label,
                            self.globaltt['sequence_alteration'])
                        feat.addFeatureStartLocation(None, chr_id)
                        feat.addFeatureToGraph()
                        geno.addParts(
                            karyotype_feature_id, karyotype_id,
                            self.globaltt['has_variant_part'])

                gene = row[col.index('gene')].strip()
                mutation = row[col.index('mutation')].strip()
                if gene != '':
                    varl = gene + '(' + mutation + ')'

                # fix the variant_id so it's always in the same order
                variant_id = row[col.index('variant_id')].strip()
                vids = variant_id.split(';')
                variant_id = ';'.join(sorted(list(set(vids))))

                if karyotype.strip() != '' and not self._is_normal_karyotype(
                        karyotype):

                    gvc_id = karyotype_id
                    if variant_id != '':
                        gvc_id = '_:' + variant_id.replace(';', '-') + '-' \
                            + re.sub(r'\w*:', '', karyotype_id)
                    if mutation.strip() != '':
                        gvc_label = '; '.join((varl, karyotype))
                    else:
                        gvc_label = karyotype
                elif variant_id.strip() != '':
                    gvc_id = '_:' + variant_id.replace(';', '-')
                    gvc_label = varl
                else:
                    # wildtype?
                    pass

                # add the karyotype to the gvc.
                # use reference if normal karyotype
                karyo_rel = self.globaltt['has_variant_part']
                if self._is_normal_karyotype(karyotype):
                    karyo_rel = self.globaltt['has_reference_part']
                if karyotype_id is not None \
                        and not self._is_normal_karyotype(karyotype) \
                        and gvc_id is not None and karyotype_id != gvc_id:
                    geno.addParts(karyotype_id, gvc_id, karyo_rel)

                if variant_id.strip() != '':
                    # split the variants & add them as part of the genotype
                    # we don't necessarily know their zygosity,
                    # just that they are part of the genotype variant ids
                    # are from OMIM, so prefix as such we assume that the
                    # sequence alts will be defined in OMIM not here
                    # TODO sort the variant_id list, if the omim prefix is
                    # the same, then assume it's the locus make a hashmap
                    # of the omim id to variant id list;
                    # then build the genotype hashmap is also useful for
                    # removing the "genes" from the list of "phenotypes"

                    # will hold gene/locus id to variant list
                    omim_map = {}

                    locus_num = None
                    for var in variant_id.split(';'):
                        # handle omim-style and odd var ids
                        # like 610661.p.R401X
                        mch = re.match(r'(\d+)\.+(.*)', var.strip())
                        if mch is not None and len(mch.groups()) == 2:
                            (locus_num, var_num) = mch.groups()

                        if locus_num is not None and locus_num not in omim_map:
                            omim_map[locus_num] = [var_num]
                        else:
                            omim_map[locus_num] += [var_num]

                    for omim in omim_map:
                        # gene_id = 'OMIM:' + omim  # TODO unused
                        vslc_id = '_:' + '-'.join(
                            [omim + '.' + a for a in omim_map.get(omim)])
                        vslc_label = varl
                        # we don't really know the zygosity of
                        # the alleles at all.
                        # so the vslcs are just a pot of them
                        model.addIndividualToGraph(
                            vslc_id, vslc_label,
                            self.globaltt['variant single locus complement'])
                        for var in omim_map.get(omim):
                            # this is actually a sequence alt
                            allele1_id = 'OMIM:' + omim + '.' + var
                            geno.addSequenceAlteration(allele1_id, None)

                            # assume that the sa -> var_loc -> gene
                            # is taken care of in OMIM
                            geno.addPartsToVSLC(
                                vslc_id, allele1_id, None,
                                self.globaltt['indeterminate'],
                                self.globaltt['has_variant_part'])

                        if vslc_id != gvc_id:
                            geno.addVSLCtoParent(vslc_id, gvc_id)

                if affected == 'unaffected':
                    # let's just say that this person is wildtype
                    model.addType(patient_id, self.globaltt['wildtype'])
                elif genotype_id is None:
                    # make an anonymous genotype id (aka blank node)
                    genotype_id = '_:geno' + catalog_id.strip()

                # add the gvc
                if gvc_id is not None:
                    model.addIndividualToGraph(
                        gvc_id, gvc_label,
                        self.globaltt['genomic_variation_complement'])

                    # add the gvc to the genotype
                    if genotype_id is not None:
                        if affected == 'unaffected':
                            rel = self.globaltt['has_reference_part']
                        else:
                            rel = self.globaltt['has_variant_part']
                        geno.addParts(gvc_id, genotype_id, rel)

                    if karyotype_id is not None \
                            and self._is_normal_karyotype(karyotype):
                        if gvc_label is not None and gvc_label != '':
                            genotype_label = '; '.join((gvc_label, karyotype))
                        elif karyotype is not None:
                            genotype_label = karyotype
                        if genotype_id is None:
                            genotype_id = karyotype_id
                        else:
                            geno.addParts(
                                karyotype_id, genotype_id,
                                self.globaltt['has_reference_part'])
                    else:
                        genotype_label = gvc_label
                        # use the catalog id as the background
                    genotype_label += ' ['+catalog_id.strip()+']'

                if genotype_id is not None and gvc_id is not None:
                    # only add the genotype if it has some parts
                    geno.addGenotype(
                        genotype_id, genotype_label,
                        self.globaltt['intrinsic_genotype'])
                    geno.addTaxon(taxon, genotype_id)
                    # add that the patient has the genotype
                    # TODO check if the genotype belongs to
                    # the cell line or to the patient
                    graph.addTriple(
                        patient_id, self.globaltt['has_genotype'], genotype_id)
                else:
                    geno.addTaxon(taxon, patient_id)

                # TODO: Add sex/gender  (as part of the karyotype?)
                # = row[col.index('')].strip()
                # #############    DEAL WITH THE DISEASES   #############
                omim_num = row[col.index('omim_num')].strip()

                # we associate the disease to the patient
                if affected == 'affected' and omim_num != '':
                    for disease in omim_num.split(';'):
                        if disease is not None and disease != '':
                            # if the omim number is in omim_map,
                            # then it is a gene not a pheno

                            # TEC - another place to use the mimTitle omim
                            # classifier omia & genereviews are using

                            if disease not in omim_map:
                                disease_id = 'OMIM:' + disease.strip()
                                # assume the label is taken care of in OMIM
                                model.addClassToGraph(disease_id, None)

                                # add the association:
                                #   the patient has the disease
                                assoc = G2PAssoc(
                                    graph, self.name,
                                    patient_id, disease_id)
                                assoc.add_association_to_graph()

                                # this line is a model of this disease
                                # TODO abstract out model into
                                # it's own association class?
                                graph.addTriple(
                                    cell_line_id,
                                    self.globaltt['is model of'],
                                    disease_id)
                            else:
                                LOG.info('drop gene %s from disease list', disease)

                # #############    ADD PUBLICATIONS   #############
                pubmed_ids = row[col.index('pubmed_ids')].strip()
                if pubmed_ids != '':
                    for pmid in pubmed_ids.split(';'):
                        pubmed_id = 'PMID:' + pmid.strip()
                        ref = Reference(graph, pubmed_id)
                        ref.setType(self.globaltt['journal article'])
                        ref.addRefToGraph()
                        graph.addTriple(
                            pubmed_id, self.globaltt['mentions'], cell_line_id)

                if not self.test_mode and (
                        limit is not None and line_counter > limit):
                    break
        return

    def _process_collection(self, collection_id, label, page):
        """
        This function will process the data supplied internally
        about the repository from Coriell.

        Triples:
            Repository a ERO:collection
            rdf:label Literal(label)
            foaf:page Literal(page)

        :param collection_id:
        :param label:
        :param page:
        :return:
        """
        # #############    BUILD THE CELL LINE REPOSITORY    #############
        for graph in [self.graph, self.testgraph]:
            # TODO: How to devise a label for each repository?
            model = Model(graph)
            reference = Reference(graph)
            repo_id = 'CoriellCollection:' + collection_id
            repo_label = label
            repo_page = page

            model.addIndividualToGraph(
                repo_id, repo_label, self.globaltt['collection'])
            reference.addPage(repo_id, repo_page)

        return

    @staticmethod
    def _get_affected_chromosomes_from_karyotype(karyotype):

        affected_chromosomes = set()
        # precompile these?
        chr_regex = r'(\d+|X|Y|M|\?);?'
        abberation_regex = r'(?:add|del|der|i|idic|inv|r|rec|t)\([\w;]+\)'
        sex_regex = r'(?:;)(X{2,}Y+|X?Y{2,}|X{3,}|X|Y)(?:;|$)'

        # first fetch the set of abberations
        abberations = re.findall(abberation_regex, karyotype)

        # iterate over them to get the chromosomes
        for abbv in abberations:
            chrs = re.findall(chr_regex, abbv)
            affected_chromosomes = affected_chromosomes.union(set(chrs))

        # remove the ? as a chromosome, since it isn't valid
        if '?' in affected_chromosomes:
            affected_chromosomes.remove('?')

        # check to see if there are any abnormal sex chromosomes
        mtch = re.search(sex_regex, karyotype)
        if mtch is not None:
            if re.search(r'X?Y{2,}', mtch.group(1)):
                # this is the only case where there is an extra Y chromosome
                affected_chromosomes.add('Y')
            else:
                affected_chromosomes.add('X')

        return affected_chromosomes

    @staticmethod
    def _is_normal_karyotype(karyotype):
        """
        This will default to true if no karyotype is provided.
        This is assuming human karyotypes.
        :param karyotype:
        :return:
        """

        is_normal = True
        if karyotype is not None:
            karyotype = karyotype.strip()
            if karyotype not in ['46;XX', '46;XY', '']:
                is_normal = False

        return is_normal

    def getTestSuite(self):
        import unittest
        from tests.test_coriell import CoriellTestCase
        # TODO add G2PAssoc, Genotype tests

        test_suite = unittest.TestLoader().loadTestsFromTestCase(CoriellTestCase)

        return test_suite
