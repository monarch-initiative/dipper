import csv
import os
from datetime import datetime
from stat import *
import re
import logging
from Bio.Seq import Seq

from dipper.utils import pysed
from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Genotype import Genotype
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)

class ZFIN(Source):
    # TODO: Enter a description for the resource.
    """
    Notes/Description for ZFIN here.

    """

    files = {
        'geno': {'file': 'genotype_features.txt', 'url': 'http://zfin.org/downloads/genotype_features.txt'},
        'pheno': {'file': 'phenotype.txt', 'url': 'http://zfin.org/downloads/phenotype.txt'},
        'pubs': {'file': 'zfinpubs.txt', 'url': 'http://zfin.org/downloads/zfinpubs.txt'},
        'zpmap': {'file': 'zp-mapping.txt',
                  'url': 'https://phenotype-ontologies.googlecode.com/svn/trunk/src/ontology/zp/zp-mapping.txt'},
        'morph': {'file': 'Morpholinos.txt', 'url': 'http://zfin.org/downloads/Morpholinos.txt'},
        'enviro': {'file': 'pheno_environment.txt', 'url': 'http://zfin.org/Downloads/pheno_environment.txt'},
        'stage': {'file': 'stage_ontology.txt', 'url': 'http://zfin.org/Downloads/stage_ontology.txt'},
        'anatomy': {'file': 'anatomy_item.txt', 'url': 'http://zfin.org/Downloads/anatomy_item.txt'},
        'wild_expression': {'file' : 'wildtype-expression.txt',
                            'url': 'http://zfin.org/Downloads/wildtype-expression.txt'},
        'mappings': {'file': 'mappings.txt', 'url': 'http://zfin.org/downloads/mappings.txt'},
        'backgrounds': {'file': 'genotype_backgrounds.txt',
                        'url': 'http://zfin.org/downloads/genotype_backgrounds.txt'},
        'genbank': {'file': 'genbank.txt', 'url': 'http://zfin.org/downloads/genbank.txt'},
        'uniprot': {'file': 'uniprot.txt', 'url': 'http://zfin.org/downloads/uniprot.txt'},
        'gene': {'file': 'gene.txt', 'url': 'http://zfin.org/downloads/gene.txt'},
        'wild': {'file': 'wildtypes.txt', 'url': 'http://zfin.org/downloads/wildtypes.txt'},
        'human_orthos': {'file': 'human_orthos.txt', 'url': 'http://zfin.org/downloads/human_orthos.txt'},
        'features': {'file': 'features.txt', 'url': 'http://zfin.org/downloads/features.txt'},
        'feature_affected_gene': {'file': 'features-affected-genes.txt',
                                  'url' : 'http://zfin.org/downloads/features-affected-genes.txt'},
        'gene_marker_rel': {'file': 'gene_marker_relationship.txt',
                            'url': 'http://zfin.org/downloads/gene_marker_relationship.txt'},
        'crispr': {'file': 'CRISPR.txt', 'url': 'http://zfin.org/downloads/CRISPR.txt'},
        'talen': {'file': 'TALEN.txt', 'url': 'http://zfin.org/downloads/TALEN.txt'},
        'pub2pubmed': {'file': 'pub_to_pubmed_id_translation.txt',
                       'url': 'http://zfin.org/downloads/pub_to_pubmed_id_translation.txt'}
    }

    # I do not love putting these here; but I don't know where else to put them
    test_ids = {
        "genotype": ["ZDB-GENO-010426-2", "ZDB-GENO-010427-3", "ZDB-GENO-010427-4", "ZDB-GENO-050209-30",
                      "ZDB-GENO-051018-1", "ZDB-GENO-070209-80", "ZDB-GENO-070215-11", "ZDB-GENO-070215-12",
                      "ZDB-GENO-070228-3", "ZDB-GENO-070406-1", "ZDB-GENO-070712-5", "ZDB-GENO-070917-2",
                      "ZDB-GENO-080328-1", "ZDB-GENO-080418-2", "ZDB-GENO-080516-8", "ZDB-GENO-080606-609",
                      "ZDB-GENO-080701-2", "ZDB-GENO-080713-1", "ZDB-GENO-080729-2", "ZDB-GENO-080804-4",
                      "ZDB-GENO-080825-3", "ZDB-GENO-091027-1", "ZDB-GENO-091027-2", "ZDB-GENO-091109-1",
                      "ZDB-GENO-100325-3", "ZDB-GENO-100325-4", "ZDB-GENO-100325-5", "ZDB-GENO-100325-6",
                      "ZDB-GENO-100524-2", "ZDB-GENO-100601-2", "ZDB-GENO-100910-1", "ZDB-GENO-111025-3",
                      "ZDB-GENO-120522-18", "ZDB-GENO-121210-1", "ZDB-GENO-130402-5", "ZDB-GENO-980410-268",
                      "ZDB-GENO-080307-1", "ZDB-GENO-960809-7", "ZDB-GENO-990623-3"],
        "gene": ["ZDB-GENE-000616-6", "ZDB-GENE-000710-4", "ZDB-GENE-030131-2773", "ZDB-GENE-030131-8769",
                  "ZDB-GENE-030219-146", "ZDB-GENE-030404-2", "ZDB-GENE-030826-1", "ZDB-GENE-030826-2",
                  "ZDB-GENE-040123-1", "ZDB-GENE-040426-1309", "ZDB-GENE-050522-534", "ZDB-GENE-060503-719",
                  "ZDB-GENE-080405-1", "ZDB-GENE-081211-2", "ZDB-GENE-091118-129", "ZDB-GENE-980526-135",
                  "ZDB-GENE-980526-166", "ZDB-GENE-980526-196", "ZDB-GENE-980526-265", "ZDB-GENE-980526-299",
                  "ZDB-GENE-980526-41", "ZDB-GENE-980526-437", "ZDB-GENE-980526-44", "ZDB-GENE-980526-481",
                  "ZDB-GENE-980526-561", "ZDB-GENE-980526-89", "ZDB-GENE-990415-181", "ZDB-GENE-990415-72",
                  "ZDB-GENE-990415-75"],
        "allele": ["ZDB-ALT-010426-4", "ZDB-ALT-010427-8", "ZDB-ALT-011017-8", "ZDB-ALT-051005-2", "ZDB-ALT-051227-8",
                    "ZDB-ALT-060221-2", "ZDB-ALT-070314-1", "ZDB-ALT-070409-1", "ZDB-ALT-070420-6", "ZDB-ALT-080528-1",
                    "ZDB-ALT-080528-6", "ZDB-ALT-080827-15", "ZDB-ALT-080908-7", "ZDB-ALT-090316-1", "ZDB-ALT-100519-1",
                    "ZDB-ALT-111024-1", "ZDB-ALT-980203-1374", "ZDB-ALT-980203-412", "ZDB-ALT-980203-465",
                    "ZDB-ALT-980203-470", "ZDB-ALT-980203-605", "ZDB-ALT-980413-636"],
        "morpholino": ["ZDB-MRPHLNO-041129-1", "ZDB-MRPHLNO-041129-2", "ZDB-MRPHLNO-041129-3", "ZDB-MRPHLNO-050308-1",
                        "ZDB-MRPHLNO-050308-3", "ZDB-MRPHLNO-060508-2", "ZDB-MRPHLNO-070118-1", "ZDB-MRPHLNO-070522-3",
                        "ZDB-MRPHLNO-070706-1", "ZDB-MRPHLNO-070725-1", "ZDB-MRPHLNO-070725-2", "ZDB-MRPHLNO-071005-1",
                        "ZDB-MRPHLNO-071227-1", "ZDB-MRPHLNO-080307-1", "ZDB-MRPHLNO-080428-1", "ZDB-MRPHLNO-080430-1",
                        "ZDB-MRPHLNO-080919-4", "ZDB-MRPHLNO-081110-3", "ZDB-MRPHLNO-090106-5", "ZDB-MRPHLNO-090114-1",
                        "ZDB-MRPHLNO-090505-1", "ZDB-MRPHLNO-090630-11", "ZDB-MRPHLNO-090804-1", "ZDB-MRPHLNO-100728-1",
                        "ZDB-MRPHLNO-100823-6", "ZDB-MRPHLNO-101105-3", "ZDB-MRPHLNO-110323-3", "ZDB-MRPHLNO-111104-5",
                        "ZDB-MRPHLNO-130222-4"],
        "environment": ["ZDB-EXP-050202-1", "ZDB-EXP-071005-3", "ZDB-EXP-071227-14", "ZDB-EXP-080428-1",
                         "ZDB-EXP-080428-2", "ZDB-EXP-080501-1", "ZDB-EXP-080805-7", "ZDB-EXP-080806-5",
                         "ZDB-EXP-080806-8", "ZDB-EXP-080806-9", "ZDB-EXP-081110-3", "ZDB-EXP-090505-2",
                         "ZDB-EXP-100330-7", "ZDB-EXP-100402-1", "ZDB-EXP-100402-2", "ZDB-EXP-100422-3",
                         "ZDB-EXP-100511-5", "ZDB-EXP-101025-12", "ZDB-EXP-101025-13", "ZDB-EXP-110926-4",
                         "ZDB-EXP-110927-1", "ZDB-EXP-120809-5", "ZDB-EXP-120809-7", "ZDB-EXP-120809-9",
                         "ZDB-EXP-120913-5", "ZDB-EXP-130222-13", "ZDB-EXP-130222-7", "ZDB-EXP-130904-2"],
        "pub": ["PMID:11566854", "PMID:12588855", "PMID:12867027", "PMID:14667409", "PMID:15456722",
                 "PMID:16914492", "PMID:17374715", "PMID:17545503", "PMID:17618647", "PMID:17785424",
                 "PMID:18201692", "PMID:18358464", "PMID:18388326", "PMID:18638469", "PMID:18846223",
                 "PMID:19151781", "PMID:19759004", "PMID:19855021", "PMID:20040115", "PMID:20138861",
                 "PMID:20306498", "PMID:20442775", "PMID:20603019", "PMID:21147088", "PMID:21893049",
                 "PMID:21925157", "PMID:22718903", "PMID:22814753", "PMID:22960038", "PMID:22996643",
                 "PMID:23086717", "PMID:23203810", "PMID:23760954", "ZFIN:ZDB-PUB-140303-33",
                 "ZFIN:ZDB-PUB-140404-9"]

    }


    def __init__(self):
        Source.__init__(self, 'zfin')

        # update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('zfin', 'ZFIN', 'http://www.zfin.org',None,'http://zfin.org/warranty.html')

        # source-specific warnings.  will be cleared when resolved.
        logger.warn("We are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self, is_dl_forced):

        # fetch all the files
        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],
                                ('/').join((self.rawdir, file['file'])),
                                is_dl_forced)
            self.dataset.setFileAccessUrl(file['url'])
            # zfin versions are set by the date of download.
            st = os.stat(('/').join((self.rawdir, file['file'])))
        self.scrub()

        # this will set the version based on the last-ingested file.
        # TODO should be a date-stamp for each file?  how to track that prov?
        self.dataset.setVersion(datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d"))

        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:
        * remove oddities where there are "\" instead of empty strings
        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        pysed.replace("\\\\", '', ('/').join((self.rawdir,self.files['geno']['file'])))
        #FIXME: Trying to scrub "\" from the pheno_environment.txt file fails due to an oddity with the file type.
        #pysed.replace("\\\\", '', ('/').join((self.rawdir,self.files['enviro']['file'])))
        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            logger.info("Only parsing first %s rows of each file", limit)
        logger.info("Parsing files...")

        self._load_zp_mappings()

        if self.testOnly:
            self.testMode = True

        #self._process_mappings(limit)
        self._process_genotype_features(limit)
        self._process_g2p(limit)
        #self._process_genes(limit)
        #self._process_genbank_ids(limit)
        #self._process_uniprot_ids(limit)
        self._process_pubinfo(limit)
        #self._process_pub2pubmed(limit)
        #self._process_morpholinos(limit)
        #self._process_talens(limit)
        #self._process_crisprs(limit)
        #self._process_pheno_enviro(limit)

        logger.info("Finished parsing.")

        self.load_bindings()

        return

    def _process_genotype_features(self, limit=None):
        """
        We don't actually know the allele pairs, so we'll store some info in a hashmap for post-processing
        :param limit:
        :return:
        """
        raw = ('/').join((self.rawdir,self.files['geno']['file']))
        out = self.outfile

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno_hash = {}
        logger.info("Processing Genotypes")
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (genotype_id, genotype_name, genotype_unique_name, allele_id, allele_name, allele_ab,
                 allele_type, allele_disp_type, gene_symbol, gene_id, zygosity,
                 construct_name, construct_id, other) = row

                if self.testMode and genotype_id not in self.test_ids['genotype']:
                    continue

                genotype_id = 'ZFIN:' + genotype_id.strip()
                geno = Genotype(g)
                gt = geno.addGenotype(genotype_id, genotype_name)
                if genotype_id not in geno_hash:
                    geno_hash[genotype_id] = {};
                genoparts = geno_hash[genotype_id]

                # reassign the allele_type to a proper GENO or SO class
                allele_type = self._map_allele_type_to_geno(allele_type)

                allele_id = 'ZFIN:' + allele_id.strip()
                geno.addAllele(allele_id, allele_name, allele_type)

                if (gene_id is not None and gene_id.strip() != ''):
                    gene_id = 'ZFIN:' + gene_id.strip()
                    geno.addGene(gene_id, gene_symbol)

                    # if it's a transgenic construct, then we'll have to get the other bits
                    if (construct_id is not None and construct_id.strip() != ''):
                        construct_id = 'ZFIN:' + construct_id.strip()
                        geno.addDerivesFrom(allele_id, construct_id)


                    # allele to gene
                    geno.addAlleleOfGene(allele_id, gene_id)
                    if gene_id not in genoparts:
                        genoparts[gene_id] = [allele_id]
                    else:
                        genoparts[gene_id].append(allele_id)

                    if (zygosity == 'homozygous'):
                        genoparts[gene_id].append(allele_id)  #add the allele again
                    elif (zygosity == 'unknown'):
                        genoparts[gene_id].append('?')  #we'll just use this as a convention for unknown
                    #elif (zygosity == 'complex'):  #what are these?
                    #    genoparts[gene_id].append('complex')
                    geno_hash[genotype_id] = genoparts
                else:
                    # if the gene is not known, still need to add the allele to the genotype hash
                    # these will be added as sequence alterations.
                    genoparts[allele_id] = [allele_id]
                    if zygosity == 'homozygous':
                        genoparts[allele_id].append(allele_id)
                    elif zygosity == 'unknown':
                        genoparts[allele_id].append('?')
                    #elif zygosity == 'complex':  #not sure what to do with these?
                    #    genoparts[allele_id].append('complex')

                    geno_hash[allele_id] = genoparts

                if (limit is not None and line_counter > limit):
                    break

                # end loop through file

            # now loop through the geno_hash, and build the vslcs
            for gt in geno_hash:
                for gene_id in geno_hash.get(gt):
                    variant_locus_parts = geno_hash.get(gt).get(gene_id)
                    if gene_id in variant_locus_parts:
                        # reset the gene_id to none
                        gene_id = None

                    allele1_id = variant_locus_parts[0]
                    allele2_id = None
                    zygosity_id = None
                    # making the assumption that there are not more than 2 variant_locus_parts
                    if len(variant_locus_parts) > 1:
                        allele2_id = variant_locus_parts[1]
                    if allele2_id is not None:
                        if allele2_id == '?':
                            zygosity_id = geno.zygosity['indeterminate']
                            allele2_id = None
                        elif allele2_id == 'complex':
                            pass #not sure what to assign here
                        elif allele1_id != allele2_id:
                            zygosity_id = geno.zygosity['heterozygous']
                        elif allele1_id == allele2_id:
                            zygosity_id = geno.zygosity['homozygous']
                    else:
                        zygosity_id = geno.zygosity['indeterminate']

                    # create the vslc
                    if gene_id is None:
                        g = ''
                    else:
                        g = gene_id
                    if (allele2_id is None):
                        a2 = ''
                    else:
                        a2 = allele2_id

                    vslc_id = self.make_id(('-').join((g,allele1_id,a2)))
                    geno.addPartsToVSLC(vslc_id,allele1_id,allele2_id,zygosity_id)
                    geno.addVSLCtoParent(vslc_id,gt)


            logger.info("Done with genotypes")
        return

    def _map_allele_type_to_geno(self, allele_type):
        type = 'SO:0001059'  # default: sequence_alteration
        type_map = {
            'complex_substitution': 'SO:1000005',  # complex substitution
            'deficiency': 'SO:1000029',  # incomplete chromosome
            'deletion': 'SO:0000159',  # deletion
            'indel': 'SO:1000032',  # indel
            'insertion': 'SO:0000667',  # insertion
            'point_mutation': 'SO:1000008',  # point_mutation
            'sequence_variant': 'SO:0001060',  # sequence variant
            'transgenic_insertion': 'SO:0001218',  # transgenic insertion
            'transgenic_unspecified': 'SO:0000781',  # transgenic unspecified
            'transloc': 'SO:0000199',  # translocation
            'unspecified' : 'SO:0001059'  # sequence alteration
        }
        if (allele_type.strip() in type_map):
            type = type_map.get(allele_type)
        else:
            # TODO add logging
            logger.error("Allele Type (%s) not mapped", allele_type)

        return type

    def _process_g2p(self, limit=None):
        """
        This module currently filters for only wild-type environments, which clearly excludes application
        of morpholinos.  Very stringent filter.  To be updated at a later time.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        """
        logger.info("Processing G2P")
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        missing_zpids = list()
        mapped_zpids = list()

        # hardcode
        eco_id = "ECO:0000059"  #experimental_phenotypic_evidence
        raw = ('/').join((self.rawdir,self.files['pheno']['file']))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (genotype_id, genotype_name,
                 start_stage_id, start_stage_name,
                 end_stage_id, end_stage_name,
                 subterm1_id, subterm1_name,
                 postcomp1_rel_id, postcomp1_rel_name,
                 superterm1_id, superterm1_name,
                 quality_id, quality_name, modifier,
                 subterm2_id, subterm2_name,
                 postcomp2_rel_id, postcomp2_rel_name,
                 superterm2_id, superterm2_name,
                 pub_id, env_id, empty) = row

                if self.testMode and genotype_id not in self.test_ids['genotype']\
                        and env_id not in self.test_ids['environment']:
                    continue


                # deal with environments
                #FIXME i am only dealing with 'wild-type' environments for now
                if (not re.match('ZDB-EXP-041102-1', env_id)):
                    logger.info("Skipping non-wildtype environment %s for %s", env_id, genotype_id)
                    continue

                genotype_id = 'ZFIN:' + genotype_id.strip()

                geno = Genotype(g)
                geno.addGenotype(genotype_id,genotype_name)
                # because we are using only w.t. environments, the genotype is just intrinsic.

                phenotype_id = self._map_sextuple_to_phenotype(superterm1_id, subterm1_id, quality_id,
                                                               superterm2_id, subterm2_id, modifier)

                if (phenotype_id is None):
                    #add this phenotype to a set to report at the end
                    missing_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier])
                    continue
                else:
                    mapped_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier])

                # add abnormal phenotypes
                if (not re.match('^normal', modifier)):
                    assoc_id = self.make_id((genotype_id+env_id+phenotype_id+pub_id))
                    pub_id = 'ZFIN:' + pub_id.strip()
                    assoc = G2PAssoc(assoc_id, genotype_id, phenotype_id, pub_id, eco_id)
                    assoc.addAssociationNodeToGraph(g)
                else:
                    # add normal phenotypes
                    logger.warn("Found normal phenotype; skipping for now")

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        myset = set([','.join(x) for x in mapped_zpids])
        logger.info("%d phenotype-sextuples were mapped",len(myset))
        myset = set([','.join(x) for x in missing_zpids])
        logger.warn("The following %d phenotype-sextuples were not mapped:\n,%s",len(myset),str(myset))

        Assoc().loadAllProperties(g)

        return

    def _process_genes(self, limit=None):
        """

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
        raw = ('/').join((self.rawdir,self.files['gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, ncbi_gene_id, empty) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue


                #FIXME: Is my approach with geno here correct?
                geno = Genotype(g)

                gene_id = 'ZFIN:'+gene_id.strip()
                ncbi_gene_id = 'NCBIGene:'+ncbi_gene_id.strip()

                geno.addGene(gene_id,gene_symbol)
                gu.addEquivalentClass(g, gene_id, ncbi_gene_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with genes")
        return


    def _process_genbank_ids(self, limit=None):
        """
        This file contains BACs, cDNAs, engineered foreign genes, ESTs, engineered plasmids, Fosmids, pseudogenes,
        engineered plasmids, P1 artificial chromosomes, SSLPs, and STS's in addition to genes, maps all to GenBank.
        :param limit:
        :return:
        """
        #TODO: Test the output, make sure the GenBank URI resolves for all construct types.
        # (It does, although ESTs redirect to http://www.ncbi.nlm.nih.gov/nucest/)

        # FIXME: Is this method unnecessary once the ZFIN gene ID has been mapped to the NCBIGene ID in process_genes?
        logger.info("Processing GenBank IDs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['genbank']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (zfin_id, so_id, symbol, genbank_id, empty) = row

                if self.testMode and zfin_id not in self.test_ids['gene']:
                    continue


                # FIXME: Is my approach with geno here correct?
                geno = Genotype(g)

                zfin_id = 'ZFIN:'+zfin_id.strip()
                genbank_id = 'GenBank:'+genbank_id.strip()
                if re.match('ZFIN:ZDB-GENE.*', zfin_id):
                    geno.addGene(zfin_id, symbol)
                    gu.addClassToGraph(g, genbank_id, symbol, so_id)
                    gu.addEquivalentClass(g,zfin_id, genbank_id)
                else:
                    geno.addConstruct(zfin_id, symbol, so_id)
                    gu.addIndividualToGraph(g, genbank_id, symbol, so_id)
                    gu.addSameIndividual(g, zfin_id, genbank_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with GenBank IDs")
        return



    def _process_uniprot_ids(self, limit=None):
        """

        :param limit:
        :return:
        """
        #FIXME: Is this method unnecessary once the ZFIN gene ID has been mapped to the NCBIGene ID in process_genes?
        logger.info("Processing UniProt IDs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir,self.files['uniprot']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id,gene_so_id,gene_symbol,uniprot_id,empty) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue


                #FIXME: Is my approach with geno here correct?
                geno = Genotype(g)

                gene_id = 'ZFIN:'+gene_id.strip()
                uniprot_id = 'UniProtKB:'+uniprot_id.strip()


                # FIXME: Need to lookup with a hash whether or not the gene already exists in the graph?
                # Or just create the gene as a class, although it would be redundant?
                geno.addGene(gene_id,gene_symbol)
                # Need to add some type of 'has_gene_product' relationship here.
                # TODO: Abstract to one of the model utilities
                gu.addIndividualToGraph(g,uniprot_id,None,'SO:0000104')
                gu.addTriple(g,gene_id,gu.properties['has_gene_product'],uniprot_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with UniProt IDs")
        return



    def _process_pubinfo(self, limit=None):
        '''
        This will pull the zfin internal publication information, and map them to their equivalent
        pmid, and make labels.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        '''
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir,self.files['pubs']['file']))
        with open(raw, 'r', encoding="latin-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pub_id, pubmed_id, authors, title, journal, year, vol, pages, empty) = row

                if self.testMode and (pub_id.replace('ZFIN:','') not in self.test_ids['pub']
                    or pubmed_id.replace('PMID:','') not in self.test_ids['pub']):
                    continue


                pub_id = 'ZFIN:'+pub_id.strip()
                pub_label = ('; ').join((authors, title, journal, year, vol, pages))
                gu.addIndividualToGraph(g,pub_id,pub_label)


                if (pubmed_id != '' and pubmed_id is not None):
                    pubmed_id = 'PMID:'+pubmed_id.strip()
                    gu.addIndividualToGraph(g,pubmed_id,None)
                    gu.addSameIndividual(g,pub_id,pubmed_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_pub2pubmed(self, limit=None):
        '''
        This will pull the zfin internal publication to pubmed mappings. Somewhat redundant with the
        process_pubinfo method, but this mapping includes additional internal pub to pubmed mappings.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        '''
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir,self.files['pub2pubmed']['file']))
        with open(raw, 'r', encoding="latin-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pub_id, pubmed_id,empty) = row

                if self.testMode and (pub_id.replace('ZFIN:', '') not in self.test_ids['pub']
                    or pubmed_id.replace('PMID:', '') not in self.test_ids['pub']):
                    continue

                pub_id = 'ZFIN:'+pub_id.strip()
                gu.addIndividualToGraph(g, pub_id, None)

                if (pubmed_id != '' and pubmed_id is not None):
                    pubmed_id = 'PMID:'+pubmed_id.strip()
                    gu.addIndividualToGraph(g, pubmed_id, None)
                    gu.addSameIndividual(g, pub_id, pubmed_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return


        # TODO: The G2P function is only dealing with wild-type environments, meaning just intrinsic genotypes
        # If mapping in these extrinsic modifiers, will need to adjust the G2P function as used above.

        #TODO: We have the sequence information for each of the targeting reagents. How to model?
    def _process_morpholinos(self, limit=None):
        """
        Morpholinos are knockdown reagents.
        Only the morpholino sequence is provided, so the target sequence is calculated using biopython.
        :param limit:
        :return:
        """

        logger.info("Processing Morpholinos")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['morph']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, morpholino_id, morpholino_so_id,
                morpholino_symbol, morpholino_sequence, publication, note) = row

                if self.testMode and morpholino_id not in self.test_ids['morpholino']:
                    continue

                #FIXME: Is my approach with geno here correct?
                geno = Genotype(g)
                morpholino_id = 'ZFIN:'+morpholino_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()

                # TODO: map target sequence to morpholino.
                # FIXME: Is this correct?
                # Is the reverse complement of the morpholino sequence the target sequence or, like miRNAs, is there
                # a seed sequence that is the target sequence and it is not the full reverse complement of the sequence?
                # Also, does the morpholino require the exact sequence match or can there be mismatches?
                # Take the morpholino sequence and get the reverse complement as the target sequence
                seq = Seq(morpholino_sequence)
                target_sequence = seq.reverse_complement()
                #print(seq)
                #print(target_sequence)
                #print(morpholino_id)

                geno.addGeneTargetingReagent(morpholino_id,morpholino_symbol,morpholino_so_id)
                geno.addReagentTargetedGene(morpholino_id,gene_id, gene_id)

                # Add publication
                if(publication != ''):
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g,pub_id,None)
                    gu.addTriple(g,pub_id,gu.properties['mentions'],morpholino_id)

                # Add comment?
                if(note != ''):
                    gu.addComment(g,morpholino_id,note)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        logger.info("Done with Morpholinos")
        return


    def _process_talens(self, limit=None):
        """
        TALENs are knockdown reagents
        :param limit:
        :return:
        """

        logger.info("Processing TALENs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['talen']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, talen_id, talen_so_id,
                talen_symbol, talen_target_sequence_1, talen_target_sequence_2, publication,note) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue


                #FIXME: Is my approach with geno here correct?
                geno = Genotype(g)
                talen_id = 'ZFIN:'+talen_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()

                geno.addGeneTargetingReagent(talen_id, talen_symbol, talen_so_id)
                geno.addReagentTargetedGene(talen_id, gene_id, gene_id)

                # Add publication
                if(publication != ''):
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g,pub_id,None)
                    gu.addTriple(g, pub_id, gu.properties['mentions'], talen_id)

                # Add comment?
                if(note != ''):
                    gu.addComment(g, talen_id, note)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        logger.info("Done with TALENS")
        return


    def _process_crisprs(self, limit=None):
        """
        CRISPRs are knockdown reagents.
        :param limit:
        :return:
        """
        logger.info("Processing CRISPRs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['crispr']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, crispr_id, crispr_so_id,
                crispr_symbol, crispr_target_sequence, publication, note) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue

                # FIXME: Is my approach with geno here correct?
                geno = Genotype(g)
                crispr_id = 'ZFIN:'+crispr_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()


                geno.addGeneTargetingReagent(crispr_id, crispr_symbol, crispr_so_id)
                geno.addReagentTargetedGene(crispr_id, gene_id, gene_id)

                # Add publication
                if(publication != ''):
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g, pub_id, None)
                    gu.addTriple(g, pub_id, gu.properties['mentions'], crispr_id)


                # Add comment?
                if(note != ''):
                    gu.addComment(g, crispr_id, note)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        logger.info("Done with CRISPRs")
        return


    # FIXME: It seems that the pheno_environment.txt file has an additional column that is empty.
    # Note that ZFIN indicates this file as only having 3 columns, but 6 columns exist plus the empty column.
    # Consider scrubbing the extra column?
    def _process_pheno_enviro(self, limit=None):
        """
        The pheno_environment.txt file ties experimental conditions to an environment ID.
        An environment ID may have one or more associated conditions.
        Condition groups present: chemical, CRISPR, morpholino, pH, physical, physiological, salinity, TALEN,
        temperature, and Generic-control.
        The condition column may contain knockdown reagent IDs or mixed text.
        :param limit:
        :return:
        """

        logger.info("Processing phenotype environments")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['enviro']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (environment_id,condition_group,condition,values,units,comment,empty) = row

                if self.testMode and environment_id not in self.test_ids['environment']:
                    continue

                # FIXME: For now just using the general "Environment" geno ID (GENO:0000099)
                # There are a few specific environments available in GENO, including some standard
                # zfin environments (standard salinity and temperature, heat shock (37C), etc), which
                # includes the zfin ID instead of a GENO ID for those environments.

                environment_id = 'ZFIN:'+environment_id.strip()
                if re.match('ZDB.*',condition):
                    condition = 'ZFIN:'+condition.strip()
                    gu.addIndividualToGraph(g,environment_id,condition,gu.datatype_properties['environment'],condition_group)
                else:
                    gu.addIndividualToGraph(g,environment_id,None,gu.datatype_properties['environment'],condition_group)

                # TODO: Need to wrangle a better description, alternative parsing of variables
                # (condition_group, condition, values, units, comment). Leaving as condition group for now.
                # Data is problematic with differing values (numeric values, N/A's, blanks).
                #if(comment !=''):
                    #enviro_description = condition_group+': '+condition+' at '+values+' '+units+comment
                #else:
                    #enviro_description = condition_group+': '+condition+' at '+values+' '+units
                #print(enviro_description)
                #gu.addIndividualToGraph(g,environment_id,None,gu.datatype_properties['environment'],condition_group)


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with phenotype environments")

        return


    def _process_mappings(self, limit=None):
        """
        This function imports linkage mappings of various entities to genetic locations in cM or cR.
        Entities include sequence variants, BAC ends, cDNA, ESTs, genes, PAC ends, RAPDs, SNPs, SSLPs, and  STSs.
        :param limit:
        :return:
        """

        logger.info("Processing chromosome mappings")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['mappings']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (zfin_id, symbol, so_id, panel_symbol, chromosome, location, metric, empty) = row

                if self.testMode and zfin_id not in self.test_ids['gene']+self.test_ids['allele']:
                    continue

                #FIXME: Is my approach with geno here correct?
                geno = Genotype(g)

                zfin_id = 'ZFIN:'+zfin_id.strip()
                if re.match('ZFIN:ZDB-GENE.*', zfin_id):
                    geno.addGene(zfin_id, symbol)
                    gu.addClassToGraph(g,zfin_id, symbol, so_id)
                elif re.match('ZFIN:ZDB-ALT.*', zfin_id):
                    # FIXME: Is this correct?
                    # The mappings.txt file has these typed as SO:0001060=sequence_variant, while Genotype.py so far
                    # only has addSequenceAlteration. Sequence variant is SO:0001059. Sequence alteration is below
                    # sequence feature, while sequence variant is at the top of a different tree in the hierarchy.
                    # If I'm correct to not be using the addSequenceAlteration method, this will need to be abstracted
                    # to an addSequenceVariant method in Genotype.py.
                    gu.addIndividualToGraph(g, zfin_id, symbol, so_id)
                else:
                    geno.addConstruct(zfin_id, symbol, so_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with chromosome mappings")
        return


    def _map_sextuple_to_phenotype(self, superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier):
        """
        This will take the 6-part EQ-style annotation used by ZFIN and return the ZP id.
        Currently relies on an external mapping file, but the method may be swapped out in the future
        :param superterm1_id:
        :param subterm1_id:
        :param quality_id:
        :param superterm2_id:
        :param subterm2_id:
        :param modifier:
        :return: ZP id
        """
        zp_id = None
        #FIXME hardcode
        mod_id=modifier
        # zfin uses free-text modifiers, but we need to convert them to proper PATO classes for the mapping
        modifiers = {
            'abnormal' : 'PATO:0000460',
            'normal' : 'PATO:0000461'
        }
        if (modifier in modifiers.keys()):
            mod_id = modifiers.get(modifier)

        key = self._make_zpkey(superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, mod_id)
        mapping = self.zp_map.get(key)

        if (mapping is None):
            pass
            # logger.warn("Couldn't map ZP id to %s",("_").join(
            #    (superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, mod_id)))
        else:
            zp_id = mapping['zp_id']

        return zp_id


    def _load_zp_mappings(self):
        """
        Given a file that defines the mapping between ZFIN-specific EQ definitions and the automatically
        derived ZP ids, create a mapping here.
        This may be deprecated in the future
        :return:
        """
        self.zp_map = {}
        logger.info("Loading ZP-to-EQ mappings")
        line_counter = 0
        file=('/').join((self.rawdir, self.files['zpmap']['file']))
        with open(file, 'r', encoding="utf-8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (zp_id, zp_label, superterm1_id, subterm1_id,
                 quality_id, modifier, superterm2_id, subterm2_id) = row
                key = self._make_zpkey(superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier)
                self.zp_map[key] = {
                    'zp_id': zp_id,
                    'label': zp_label,
                    'superterm1_id': superterm1_id,
                    'subterm1_id': subterm1_id,
                    'quality_id': quality_id,
                    'modifier': modifier,
                    'superterm2_id': superterm2_id,
                    'subterm2_id': subterm2_id,
                }
        logger.info("Loaded %s zp terms", self.zp_map.__len__())

        return

    def _make_zpkey(self, superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier):
        key = self.make_id(('_').join((superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier)))
        return key

    def getTestSuite(self):
        import unittest
        from tests.test_zfin import ZFINTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(ZFINTestCase)

        return test_suite
