import logging
import csv

from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Dataset import Dataset
from dipper import config
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map



logger = logging.getLogger(__name__)

class Coriell(Source):
    """
    The Coriell Catalog provided to Monarch includes metadata and descriptions of NIGMS, NINDS, NHGRI, and
    NIA cell lines.  These lines are made available for research purposes.
    Here, we create annotations for the cell lines as models of the diseases from which
    they originate.

    Notice: The Coriell catalog is delivered to Monarch in a specific format, and requires ssh rsa fingerprint
    identification.  Other groups wishing to get this data in it's raw form will need to contact Coriell
    for credentials.

    """

    files = {
        'ninds': {'file' : 'NINDS_2014-02-03_13-32-24.csv',
                  'id' : 'NINDSrepository',
                  'label': 'NINDS Human Genetics DNA and Cell line Repository',
                  'page' : 'https://catalog.coriell.org/1/NINDS'},
        'nigms': {'file' : 'NIGMS_2014-02-03_13-31-42.csv',
                  'id' : 'NIGMSrepository',
                  'label': 'NIGMS Human Genetic Cell Repository',
                  'page' : 'https://catalog.coriell.org/1/NIGMS'},
        'nia': {'file' : 'NIA_2015-02-20_16-08-04.csv',
                'id' : 'NIArepository',
                'label': 'NIA Aging Cell Repository',
                'page': 'https://catalog.coriell.org/1/NIA'}
    }


    def __init__(self):
        Source.__init__(self, 'coriell')

        self.load_bindings()

        self.dataset = Dataset('coriell', 'Coriell', 'http://ccr.coriell.org/nigms/')

        #data-source specific warnings (will be removed when issues are cleared)
        #print()

        #check if config exists; if it doesn't, error out and let user know
        if (not (('keys' in config.get_config()) and ('coriell' in config.get_config()['keys']))):
            logger.error("ERROR: not configured with FTP user/password.")
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
        host1 = 'host=\''+ config.get_config()['keys']['coriell']['host']+'\''
        user1 = 'user=\''+ config.get_config()['keys']['coriell']['user']+'\''
        passwd1 = 'passwd=\''+ config.get_config()['keys']['coriell']['password']+'\''
        #print(host1,user1,passwd1)
        #ftp = FTP(config.get_config()['keys']['coriell']['host'],config.get_config()['keys']['coriell']['user'],config.get_config()['keys']['coriell']['password'],timeout=None)
        #ftp = FTP(host1,user1,passwd1,timeout=None)

        #ftp.login()

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:

        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        #pysed.replace("\r", '', ('/').join((self.rawdir,'dv.nlx_157874_1')))

        return

    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %s rows of each file", limit)

        logger.info("Parsing files...")

        for f in ['ninds','nigms','nia']:
            file = ('/').join((self.rawdir,self.files[f]['file']))
            #self._process_repository(self.files[f]['page'],self.files[f]['label'])
            self._process_data(file, limit)

        logger.info("Finished parsing.")


        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

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
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                if not row:
                    pass
                else:
                    line_counter += 1

                    (catalog_id,description,omim_number,sample_type,cell_line_available,
                    dna_in_stock,dna_ref,gender,age,race,ethnicity,affected,karyotype,
                    relprob,mutation,gene,fam,collection,url,cat_remark,pubmed_ids) = row

                    line_id = 'Coriell:'+catalog_id.strip()

                    cell_type = self._map_cell_type(sample_type)

                    #FIXME:Including a generic label for now.
                    line_label = collection.partition(' ')[0]+'-'+catalog_id.strip()

                    # Adding the cell line as a typed individual.
                    gu.addIndividualToGraph(self.graph,line_id,line_label,cell_type)



                    #if pubmed_ids != '':
                        #for s in pubmed_ids.split(';'):
                            #gu.addSynonym(self.graph,morphology_term_id,s.strip(), gu.relationships['hasRelatedSynonym'])



                    if (limit is not None and line_counter > limit):
                        break
        return

    def _process_repository(self, page, label):
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
        #FIXME: How to devise a label for each repository?
        gu = GraphUtils(curie_map.get())
        #repo_id = id
        repo_label = label
        repo_page = page
        #gu.addClassToGraph(self.graph,repo_id,repo_label)
        gu.addIndividualToGraph(self.graph,repo_page,repo_label,'CLO:0000008')
        #gu.addPage(self.graph,repo_id,repo_page)



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
            print("ERROR: Cell type (", sample_type, ") not mapped")

        return type