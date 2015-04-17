import logging
import csv
import pysftp
import re
import unicodedata

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
    a karyotype (if specified) or a collection of variants.  Both the genotype and disease are linked
    to the patient, and the cell line is listed as derived from the patient.  The cell line is classified
    by it's cell type (which itself is linked to a tissue of origin).

    Unfortunately, the omim numbers listed in this file are both for genes and diseases; we have no way of knowing
    a priori if a designated omim number is a gene or disease; so we presently link the patient to any omim id
    via the has_phenotype relationship.

    Notice: The Coriell catalog is delivered to Monarch in a specific format, and requires ssh rsa fingerprint
    identification.  Other groups wishing to get this data in it's raw form will need to contact Coriell
    for credentials.

    """

    terms = {
        'cell_line_repository': 'CLO:0000008',
        'race': 'SIO:001015',
        'ethnic_group': 'EFO:0001799',
        'age' : 'EFO:0000246',
        'sampling_time': 'EFO:0000689',
    }



    files = {
        'ninds': {'file' : 'NINDS_latest.csv',
                  'id' : 'NINDS',
                  'label': 'NINDS Human Genetics DNA and Cell line Repository',
                  'page' : 'https://catalog.coriell.org/1/NINDS'},
        'nigms': {'file' : 'NIGMS_latest.csv',
                  'id' : 'NIGMS',
                  'label': 'NIGMS Human Genetic Cell Repository',
                  'page' : 'https://catalog.coriell.org/1/NIGMS'},
        'nia': {'file' : 'NIA_latest.csv',
                'id' : 'NIA',
                'label': 'NIA Aging Cell Repository',
                'page': 'https://catalog.coriell.org/1/NIA'},
        'nhgri': {'file' : 'NHGRI_latest.csv',
                'id' : 'NHGRI',
                'label': 'NHGRI Sample Repository for Human Genetic Research',
                'page': 'https://catalog.coriell.org/1/NIA'}
    }



    def __init__(self):
        Source.__init__(self, 'coriell')

        self.load_bindings()

        self.dataset = Dataset('coriell', 'Coriell', 'http://ccr.coriell.org/',None)

        #data-source specific warnings (will be removed when issues are cleared)
        #print()
        logger.warn('We assume that if a species is not provided, that it is a Human-derived cell line')
        logger.warn('We map all omim ids as a disease/phenotype entity, but should be fixed in the future')

        #check if config exists; if it doesn't, error out and let user know
        if (not (('keys' in config.get_config()) and ('coriell' in config.get_config()['keys']))):
            logger.error("ERROR: not configured with FTP user/password.")

#        self.properties = self.op.copy()


        return


    def fetch(self, is_dl_forced):
        """
        Be sure to have pg user/password connection details in your conf.json file, like:
        dbauth : {
        "coriell" : {"user" : "<username>", "password" : "<password>", "host" : <host>}
        }

        :param is_dl_forced:
        :return:
        """
        host1 = config.get_config()['keys']['coriell']['host']
        user1 = 'username=\''+ config.get_config()['keys']['coriell']['user']+'\''
        passwd1 = 'password=\''+ config.get_config()['keys']['coriell']['password']+'\''
        #print(host1,user1,passwd1)
        #ftp = FTP(config.get_config()['keys']['coriell']['host'],config.get_config()['keys']['coriell']['user'],config.get_config()['keys']['coriell']['password'],timeout=None)
        #ftp = FTP(host1,user1,passwd1,timeout=None)
        #ftp.login()
        #FIXME: Still resulting in login time out.
        #with pysftp.Connection(host1, user1, passwd1) as sftp:
            #with sftp.cd('public'):
                #print('success!')
                #sftp.get_r()
        #Will need to rename the files, or handle changing file names with the same beginning (NIGMS_..., etc.)

        return


    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")

        for f in ['ninds','nigms','nia','nhgri']:
            file = ('/').join((self.rawdir,self.files[f]['file']))
            self._process_repository(self.files[f]['id'],self.files[f]['label'],self.files[f]['page'])
            self._process_data(file, limit)

        logger.info("Finished parsing.")


        self.load_bindings()
        Assoc().loadAllProperties(self.graph)

        logger.info("Found %s nodes", len(self.graph))
        return

    def _process_data(self, raw, limit=None):
        """
        This function will process the data files from Coriell.



        Triples: (examples)

            :NIGMSrepository a CLO_0000008 #repository  ?
                label : NIGMS Human Genetic Cell Repository
                foaf:page https://catalog.coriell.org/0/sections/collections/NIGMS/?SsId=8

            line_id a CL_0000057,  #fibroblast line
                derives_from patient_id, uberon:Fibroblast
                part_of :NIGMSrepository
                #we also have the age_at_sampling type of property

            patient id a foaf:person, proband, OMIM:disease_id
                label: "fibroblast from patient 12345 with disease X"
                member_of family_id  #what is the right thing here?
                SIO:race EFO:caucasian  #subclass of EFO:0001799
                in_taxon NCBITaxon:9606
                dc:description Literal(remark)
                GENO:has_genotype genotype_id

            family_id a owl:NamedIndividual
                foaf:page "https://catalog.coriell.org/0/Sections/BrowseCatalog/FamilyTypeSubDetail.aspx?PgId=402&fam=2104&coll=GM"

            genotype_id a intrinsic_genotype
                GENO:has_variant_part allelic_variant_id
                #we don't necessarily know much about the genotype, other than the allelic variant.
                #also there's the sex here)

            pub_id mentions cell_line_id



        :param raw:
        :param limit:
        :return:
        """
        logger.info("Processing Data from %s",raw)
        gu = GraphUtils(curie_map.get())

        line_counter = 0
        geno = Genotype(self.graph)

        gu.loadProperties(self.graph,geno.object_properties,gu.OBJPROP)
        gu.loadAllProperties(self.graph)

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
#        with open(raw, 'r', encoding="utf-8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                if not row:
                    pass
                else:
                    line_counter += 1

                    (catalog_id,description,omim_number,sample_type,cell_line_available,
                    dna_in_stock,dna_ref,gender,age,race,ethnicity,affected,karyotype,
                    relprob,mutation,gene,family_id,collection,url,cat_remark,pubmed_ids,
                    family_member,variant_id,dbsnp_id,species) = row

                    #example:
                    #GM00003,HURLER SYNDROME,607014,Fibroblast,Yes,No,,Female,26 YR,Caucasian,,,,
                    # parent,,,39,NIGMS Human Genetic Cell Repository,
                    # http://ccr.coriell.org/Sections/Search/Sample_Detail.aspx?Ref=GM00003,
                    # 46;XX; clinically normal mother of a child with Hurler syndrome; proband not in Repository,,
                    # 2,,18343,Homo sapiens


                    ##############    BUILD REQUIRED VARIABLES    #############

                    #Make the cell line ID
                    cell_line_id = 'Coriell:'+catalog_id.strip()

                    #Map the cell/sample type
                    cell_type = self._map_cell_type(sample_type)

                    #Make a cell line label
                    line_label = collection.partition(' ')[0]+'-'+catalog_id.strip()

                    #Map the repository/collection
                    repository = self._map_collection(collection)

                    #patients are uniquely identified by one of:
                    #dbsnp id (which is == an individual haplotype)
                    #family id + family member (if present) OR   #probands are usually family member zero
                    #cell line id
                    #since some patients have >1 cell line derived from them, we must make sure that the genotype
                    #is attached to the patient, and can be inferred to the cell line
                    #examples of repeated patients are:  famid=1159, member=1; fam=152,member=1

                    #Make the patient ID

                    if family_id != '':
                        patient_id = 'MONARCH:Coriell'+family_id+'-'+family_member  #this won't resolve like this, but how to handle?
                    else:
                        patient_id = cell_line_id   #otherwise, just default to the cell line as the patient id


                    #properties of the individual patients:  sex, family id, member/relproband, description
                    #descriptions are really long and ugly SCREAMING text, so need to clean up
                    #the control cases are so odd with this labeling scheme; but we'll deal with it as-is for now.
                    short_desc = (description.split(';')[0]).capitalize()
                    if (affected=='Yes'):
                        affected = 'affected'
                    else:
                        affected = 'unaffected'
                    gender = gender.lower()
                    patient_label = (' ').join((affected,gender,relprob))
                    if (relprob=='proband'):
                        patient_label=(' ').join((patient_label,'with',short_desc))
                    else:
                        patient_label=(' ').join((patient_label,'of proband with',short_desc))

                    ##############    BUILD THE CELL LINE    #############

                    # Adding the cell line as a typed individual.
                    gu.addIndividualToGraph(self.graph,cell_line_id,line_label,cell_type)

                    # Cell line derives from patient
                    geno.addDerivesFrom(cell_line_id,patient_id)

                    # Cell line a member of repository
                    gu.addMember(self.graph,repository,cell_line_id)

                    if cat_remark !='':
                        gu.addDescription(self.graph,cell_line_id,cat_remark)


                    # Cell age_at_sampling
                    #TODO add in the age nodes once these are modeled properly in https://github.com/monarch-initiative/dipper/issues/78
                    #if (age != ''):
                        #this would give a BNode that is an instance of Age.  but i don't know how to connect
                        #the age node to the cell line?  we need to ask @mbrush.
                        #age_id = '_'+re.sub('\s+','_',age)
                        #gu.addIndividualToGraph(self.graph,age_id,age,self.terms['age'])
                        #gu.addTriple(self.graph,age_id,self.properties['has_measurement'],age,True)

                    #that the cell line is a model for a disease is dealt with below

                    ##############    BUILD THE PATIENT    #############

                    # Add the patient ID as an individual.
                    gu.addPerson(self.graph,patient_id,patient_label)
                    #TODO map relationship to proband as a class (what ontology?)

                    # Add race of patient
                    #FIXME: Adjust for subcategories based on ethnicity field
                    #EDIT: There are 743 different entries for ethnicity... Too many to map?
                    #Perhaps add ethnicity as a literal in addition to the mapped race?
                    #Need to adjust the ethnicity text (if using) to just initial capitalization as some entries:ALLCAPS

                    #TODO the race should go into the background of the individual.  this should be abstracted out to the Genotype class
                    #punting for now.
                    #if race != '':
                    #    mapped_race = self._map_race(race)
                    #    if mapped_race is not None:
                    #        gu.addTriple(self.graph,patient_id,self.terms['race'],mapped_race)
                    #        gu.addSubclass(self.graph,self.terms['ethnic_group'],mapped_race)


                    ##############    BUILD THE FAMILY    #############

                    # Add triples for family_id, if present.
                    # a family is just a small population
                    if family_id != '':
                        family_comp_id = 'CoriellFamily:'+family_id

                        # Add the family ID as a named individual
                        gu.addIndividualToGraph(self.graph,family_comp_id,None,geno.genoparts['population'])

                        #Add the patient as a member of the family
                        gu.addMemberOf(self.graph,patient_id,family_comp_id)


                    ##############    BUILD THE GENOTYPE   #############

                    #the important things to pay attention to here are:
                    #karyotype (used when there isn't a typical mutation), esp for chr rearrangements
                    #mutation = protein-level mutation as a label, often from omim
                    #gene = gene symbol - TODO get id
                    #variant_id = omim variant ids (; delimited)
                    #dbsnp_id = snp individual ids, is == patient id == haplotype?

                    #note GM00633 is a good example of chromosomal variation - do we have enough to capture this?

                    #make an assumption that if the taxon is blank, that it is human!
                    if (species is None or species == ''):
                        species = 'Homo sapiens'
                    taxon = self._map_species(species)

                    #if there's a dbSNP id, this is actually the individual's genotype
                    genotype_id = None
                    genotype_label = None
                    if (dbsnp_id != ''):
                        genotype_id = 'dbSNPIndividual:'+dbsnp_id.strip()

                    omim_map = {};
                    gvc_id = None;
                    if (karyotype.strip() != '' and karyotype != '46;XX' and karyotype != '46;XY'):
                        #if we have the karyotype that is abnormal, then let's make that the genotype
                        gvc_id = self.make_id(karyotype)

                        #some of the karyotypes are encoded with terrible hidden codes. remove them here
                        #i've seen a <98> character
                        gvc_label = self.remove_control_characters(karyotype)

                    elif (variant_id.strip() != ''):
                        #split the variants and add them as part of the genotype
                        #we don't know their zygosity, just that they are part of the genotype
                        #variant ids are from OMIM, so prefix as such
                        #we assume that the variants will be defined in OMIM rather than here.
                        gvc_id = self.make_id(variant_id)
                        gvc_label = mutation   #this is the full list of mutational consequences
                        #TODO sort the variant_id list, if the omim prefix is the same, then assume it's the locus
                        #make a hashmap of the omim id to variant id list; then build the genotype
                        #hashmap is also useful for removing the "genes" from the list of phenotypes
                        omim_map = {};
                        #print(variant_id)
                        locus_num = None
                        for v in variant_id.split(';'):
                            m = re.match('(\d+)\.+(.*)',v.strip())  #handle omim-style variant ids and odd ones like 610661.p.R401X

                            if (m is not None and len(m.groups()) == 2):
                                (locus_num,var_num) = m.groups()

                            if locus_num is not None and locus_num not in omim_map:
                                omim_map[locus_num] = [var_num]
                            else:
                                omim_map[locus_num] = omim_map[locus_num]+[var_num]

                        for o in omim_map:
                            #print(omim_map)
                            allele1_id = 'OMIM:'+o+'.'+omim_map.get(o)[0]
                            if (len(omim_map.get(o)) > 1):
                                #build a vslc
                                vslc_id = self.make_id(o+(' ').join(omim_map.get(o)))
                                allele2_id = 'OMIM:'+o+'.'+omim_map.get(o)[1]
                                geno.addPartsToVSLC(vslc_id,allele1_id,allele2_id)
                                geno.addVSLCtoParent(vslc_id,gvc_id)
                            else:
                                gu.addIndividualToGraph(self.graph,allele1_id,None,geno.genoparts['variant_locus'])
                                geno.addParts(allele1_id,gvc_id,geno.properties['has_alternate_part'])
                    #else:
                        #genotype_id = None
                        #if affected == 'affected':
                        #    logger.info('affected individual without recorded variants: %s, %s',cell_line_id,description)

                    if (affected == 'unaffected'):
                        #let's just say that this person is wildtype
                        gu.addType(self.graph,patient_id,geno.genoparts['wildtype'])
                    elif genotype_id is None:
                        #make an anonymous genotype id
                        genotype_id = '_geno'+cell_line_id


                    #add the gvc
                    if gvc_id is not None:
                        gu.addIndividualToGraph(self.graph,gvc_id,gvc_label,geno.genoparts['genomic_variation_complement'])

                        #add the gvc to the genotype
                        if genotype_id is not None:
                            if (affected == 'unaffected'):
                                rel = geno.object_properties['has_reference_part']
                            else:
                                rel = geno.object_properties['has_alternate_part']
                            geno.addParts(gvc_id,genotype_id,rel)
                        genotype_label = gvc_label


                    if genotype_id is not None:
                        geno.addGenotype(genotype_id,genotype_label,geno.genoparts['intrinsic_genotype'])
                        geno.addTaxon(taxon,genotype_id)
                        #add that the patient has the genotype
                        gu.addTriple(self.graph,patient_id,geno.properties['has_genotype'],genotype_id)
                    else:
                        geno.addTaxon(taxon,patient_id)

                    #TODO: Add sex/gender


                    ##############    DEAL WITH THE DISEASES   #############

                    #what has the disease, the genotype or patient?
                    if (affected == 'affected'):
                        if omim_number != '':
                            for d in omim_number.split(';'):
                                if d is not None and d != '':
                                    #if the omim number is in the omim_map, then it is a gene, and don't add it as a phenotype
                                    if d not in omim_map:
                                        disease_id = 'OMIM:'+d.strip()
                                        gu.addClassToGraph(self.graph,disease_id,None)   #assume the label is taken care of elsewhere

                                        #add the association: the patient has the disease
                                        assoc_id = self.make_id(patient_id+disease_id)
                                        assoc = G2PAssoc(assoc_id,patient_id,disease_id,None,None)
                                        assoc.addAssociationToGraph(self.graph)

                                        #also, this line is a model of this disease
                                        #TODO abstract out model into it's own association class?
                                        gu.addTriple(self.graph,cell_line_id,gu.properties['model_of'],disease_id)
                                    else:
                                        logger.info('removing %s from disease list since it is a gene',d)

                                        #BEWARE
                                        #Coriell currently seems to sometimes annotate variants to phenotype ids rather
                                        #than gene ids!  YIKES
                                        #we are checking with them now to find out why.
                                        #for example.  GM16467 is annotated to 278720.0006, but it should be 613208.0006
                                        #so this current solution is aggressive, and might mean missing data


                    ##############    ADD PUBLICATIONS   #############

                    if pubmed_ids != '':
                        for s in pubmed_ids.split(';'):
                            pubmed_id = 'PMID:'+s.strip()
                            gu.addTriple(self.graph,pubmed_id,gu.properties['mentions'],cell_line_id)


                    if (limit is not None and line_counter > limit):
                        break
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
        ##############    BUILD THE CELL LINE REPOSITORY    #############
        #FIXME: How to devise a label for each repository?
        gu = GraphUtils(curie_map.get())
        repo_id = 'CoriellCollection:'+id
        repo_label = label
        repo_page = page
        #gu.addClassToGraph(self.graph,repo_id,repo_label)
        gu.addIndividualToGraph(self.graph,repo_id,repo_label,self.terms['cell_line_repository'])
        gu.addPage(self.graph,repo_id,repo_page)

        return



    def _map_cell_type(self, sample_type):
        type = None
        type_map = {
            'Adipose stromal cell': 'CL:0002570', # FIXME: mesenchymal stem cell of adipose
            'Amniotic fluid-derived cell line': 'CL:0002323',  # FIXME: amniocyte?
            'B-Lymphocyte': 'CL:0000236',  # B cell
            'Chorionic villus-derived cell line': 'CL:0000000', # FIXME: No Match
            'Endothelial': 'CL:0000115',  # endothelial cell
            'Epithelial': 'CL:0000066',  # epithelial cell
            'Erythroleukemic cell line': 'CL:0000000',  # FIXME: No Match. "Abnormal precursor (virally transformed) of mouse erythrocytes that can be grown in culture and induced to differentiate by treatment with, for example, DMSO."
            'Fibroblast': 'CL:0000057',  # fibroblast
            'Keratinocyte': 'CL:0000312', # keratinocyte
            'Melanocyte': 'CL:0000148',  # melanocyte
            'Mesothelial': 'CL:0000077',
            'Microcell hybrid': 'CL:0000000',  # FIXME: No Match
            'Myoblast': 'CL:0000056',  # myoblast
            'Smooth muscle': 'CL:0000192',  # smooth muscle cell
            'Stem cell': 'CL:0000034',  # stem cell
            'T-Lymphocyte': 'CL:0000084',  # T cell
            'Tumor-derived cell line': 'CL:0002198'  # FIXME: No Match. "Cells isolated from a mass of neoplastic cells, i.e., a growth formed by abnormal cellular proliferation." Oncocyte? CL:0002198
        }
        if (sample_type.strip() in type_map):
            type = type_map.get(sample_type)
        else:
            logger.warn("ERROR: Cell type not mapped: %s", sample_type)

        return type


    def _map_race(self, race):
        type = None
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
        if (race.strip() in type_map):
            type = type_map.get(race)
        else:
            logger.warn("Race type not mapped: %s", race)

        return type


    def _map_species(self, species):
        type = None
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
            'Homo sapiens' : 'NCBITaxon:9606'
        }
        if (species.strip() in type_map):
            type = type_map.get(species)
        else:
            logger.warn("Species type not mapped: %s", species)

        return type


    def _map_collection(self, collection):
        type = None
        type_map = {
            'NINDS Repository': 'CoriellCollection:NINDS',
            'NIGMS Human Genetic Cell Repository': 'CoriellCollection:NIGMS',
            'NIA Aging Cell Culture Repository': 'CoriellCollection:NIA',
            'NHGRI Sample Repository for Human Genetic Research': 'CoriellCollection:NHGRI'
        }
        if (collection.strip() in type_map):
            type = type_map.get(collection)
        else:
            logger.warn("ERROR: Collection type not mapped: %s", collection)

        return type

    def remove_control_characters(self,s):
        return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")
