import csv
import re
import logging
# from Bio.Seq import Seq

from dipper.utils import pysed
from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Genotype import Genotype
from dipper.models.OrthologyAssoc import OrthologyAssoc
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


class ZFIN(Source):
    """
    This is the parser for the [Zebrafish Model Organism Database (ZFIN)](http://www.zfin.org),
    from which we process genotype and phenotype data for laboratory zebrafish.
    Genotypes leverage the GENO genotype model and includes both intrinsic and extrinsic genotypes.

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
        'wild_expression': {'file': 'wildtype-expression.txt',
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
                                  'url': 'http://zfin.org/downloads/features-affected-genes.txt'},
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
                     "ZDB-GENO-080307-1", "ZDB-GENO-960809-7", "ZDB-GENO-990623-3", "ZDB-GENO-130603-1",
                     "ZDB-GENO-001127-3", "ZDB-GENO-001129-1", "ZDB-GENO-090203-8"],
        "gene": ["ZDB-GENE-000616-6", "ZDB-GENE-000710-4", "ZDB-GENE-030131-2773", "ZDB-GENE-030131-8769",
                 "ZDB-GENE-030219-146", "ZDB-GENE-030404-2", "ZDB-GENE-030826-1", "ZDB-GENE-030826-2",
                 "ZDB-GENE-040123-1", "ZDB-GENE-040426-1309", "ZDB-GENE-050522-534", "ZDB-GENE-060503-719",
                 "ZDB-GENE-080405-1", "ZDB-GENE-081211-2", "ZDB-GENE-091118-129", "ZDB-GENE-980526-135",
                 "ZDB-GENE-980526-166", "ZDB-GENE-980526-196", "ZDB-GENE-980526-265", "ZDB-GENE-980526-299",
                 "ZDB-GENE-980526-41", "ZDB-GENE-980526-437", "ZDB-GENE-980526-44", "ZDB-GENE-980526-481",
                 "ZDB-GENE-980526-561", "ZDB-GENE-980526-89", "ZDB-GENE-990415-181", "ZDB-GENE-990415-72",
                 "ZDB-GENE-990415-75", "ZDB-GENE-980526-44", "ZDB-GENE-030421-3", "ZDB-GENE-980526-196",
                 "ZDB-GENE-050320-62", "ZDB-GENE-061013-403", "ZDB-GENE-041114-104", "ZDB-GENE-030131-9700",
                 "ZDB-GENE-031114-1"],
        "allele": ["ZDB-ALT-010426-4", "ZDB-ALT-010427-8", "ZDB-ALT-011017-8", "ZDB-ALT-051005-2", "ZDB-ALT-051227-8",
                   "ZDB-ALT-060221-2", "ZDB-ALT-070314-1", "ZDB-ALT-070409-1", "ZDB-ALT-070420-6", "ZDB-ALT-080528-1",
                   "ZDB-ALT-080528-6", "ZDB-ALT-080827-15", "ZDB-ALT-080908-7", "ZDB-ALT-090316-1", "ZDB-ALT-100519-1",
                   "ZDB-ALT-111024-1", "ZDB-ALT-980203-1374", "ZDB-ALT-980203-412", "ZDB-ALT-980203-465",
                   "ZDB-ALT-980203-470", "ZDB-ALT-980203-605", "ZDB-ALT-980413-636", "ZDB-ALT-021021-2",
                   "ZDB-ALT-080728-1", "ZDB-ALT-100729-1", "ZDB-ALT-980203-1560", "ZDB-ALT-001127-6",
                   "ZDB-ALT-001129-2"],
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
                        "ZDB-EXP-120913-5", "ZDB-EXP-130222-13", "ZDB-EXP-130222-7", "ZDB-EXP-130904-2",
                        "ZDB-EXP-041102-1"],
        "pub": ["PMID:11566854", "PMID:12588855", "PMID:12867027", "PMID:14667409", "PMID:15456722",
                "PMID:16914492", "PMID:17374715", "PMID:17545503", "PMID:17618647", "PMID:17785424",
                "PMID:18201692", "PMID:18358464", "PMID:18388326", "PMID:18638469", "PMID:18846223",
                "PMID:19151781", "PMID:19759004", "PMID:19855021", "PMID:20040115", "PMID:20138861",
                "PMID:20306498", "PMID:20442775", "PMID:20603019", "PMID:21147088", "PMID:21893049",
                "PMID:21925157", "PMID:22718903", "PMID:22814753", "PMID:22960038", "PMID:22996643",
                "PMID:23086717", "PMID:23203810", "PMID:23760954", "ZFIN:ZDB-PUB-140303-33",
                "ZFIN:ZDB-PUB-140404-9", "ZFIN:ZDB-PUB-080902-16"]
    }

    def __init__(self):
        Source.__init__(self, 'zfin')

        # update the dataset object with details about this resource
        self.dataset = Dataset('zfin', 'ZFIN', 'http://www.zfin.org', None,
                               'http://zfin.org/warranty.html')

        return

    def fetch(self, is_dl_forced=False):

        # fetch all the files
        # zfin versions are set by the date of download.
        self.get_files(is_dl_forced)
        self.scrub()

        return

    def scrub(self):
        """
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:
          * remove oddities where there are "\" instead of empty strings
        :return: None
        """

        # scrub file of the oddities where there are "\" instead of empty strings
        pysed.replace("\\\\", '', '/'.join((self.rawdir, self.files['geno']['file'])))

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows of each file", limit)
        logger.info("Parsing files...")

        self._load_zp_mappings()
        
        if self.testOnly:
            self.testMode = True

        self.wildtype_hash = {}  # to hold the w.t. background id:symbol mappings

        self.kd_reagent_hash = {} # to hold the reagent symbol and list of target ids

        self.label_hash = {'gene_label': {}, 'allele_label': {}, 'construct_label': {}, 'genotype_label': {},
                           'background_label': {}, 'morpholino_label': {}}
        self.id_label_map = {}
        self.genotype_backgrounds = {}  # to hold the mappings between genotype and background
        self.extrinsic_id_to_enviro_id_hash = {}
        self.variant_loci_genes = {}   #to hold the genes variant due to a seq alt

        # basic information on classes and instances
        self._process_genes(limit)   #REVIEWED-COMPLETE
        self._process_stages(limit)  #REVIEWED-COMPLETE
        self._process_pubinfo(limit)  #REVIEWED-COMPLETE
        self._process_pub2pubmed(limit)  #REVIEWED-COMPLETE

        # The knockdown reagents
        for t in ['morph','crispr','talen']:
            self._process_targeting_reagents(t, limit)  #REVIEWED-COMPLETE

        self._process_gene_marker_relationships(limit)  #REVIEW-IN PROCESS
        self._process_features(limit)  #REVIEWED-COMPLETE
        self._process_feature_affected_genes(limit)  #REVIEWED-COMPLETE

        # These must be processed in a specific order
        # Must be processed before wildtype_expression
        self._process_wildtypes(limit)    #REVIEWED-COMPLETE
        self._process_genotype_backgrounds(limit)  #REVIEWED-COMPLETE
        self._process_genotype_features(limit)

        # materialize the genotype parts (vslcs, gvcs, etc.)
        #self._build_genotype_parts(limit)  #TODO

        # NOTE: TALENs/CRISPRs are not currently added to the extrinsic genotype, but are added to the env label.
        #self._process_pheno_enviro(limit)  # Must be processed after morpholinos/talens/crisprs

        # once the genotypes and environments are processed, we can associate these with the phenotypes
        #self._process_g2p(limit)

        #self._process_wildtype_expression(limit)
        #self._process_genbank_ids(limit)  #HOLDING OFF FOR NOW
        #self._process_uniprot_ids(limit)  #HOLDING OFF FOR NOW
        #self._process_human_orthos(limit)
        #self._process_mappings(limit)  #NEEDS WORK-DO NOT INCLUDE YET

        logger.info("Finished parsing.")

        self.load_bindings()

        return

    def _process_genotype_features(self, limit=None):
        """
        Here we process the genotype_features file, which lists genotypes together with
        any intrinsic sequence alterations, their zygosity, and affected gene.
        We don't get allele pair (VSLC) ids in a single row, so we iterate through the
        file and build up a hash that contains all of a genotype's partonomy.
        We can then assemble a genotype based on that partonomy.
        This table does not list the background genotype/strain: that is listed elsewhere.

        ZFIN "ALT" objects are mapped to sequence alterations in our system.

        :param limit:
        :return:
        """
        raw = '/'.join((self.rawdir, self.files['geno']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        gu = GraphUtils(curie_map.get())

        geno_hash = {}  # This is used to store the genotype partonomy
        gvc_hash = {}

        logger.info("Processing Genotypes")
        line_counter = 0
        geno = Genotype(g)
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (genotype_num, genotype_name, genotype_unique_name, allele_num, allele_name, allele_ab,
                 allele_type, allele_disp_type, gene_symbol, gene_num, zygosity,
                 construct_name, construct_num, other) = row

                if self.testMode and genotype_num not in self.test_ids['genotype']:
                    continue

                # add the genotype to the graph
                # not adding the genotype label here, since it doesn't include the background
                # that will be done in another method
                genotype_id = 'ZFIN:' + genotype_num.strip()
                geno.addGenotype(genotype_id, None)

                # add the given name and uniquename as synonyms
                gu.addSynonym(g, genotype_id, genotype_name)
                gu.addSynonym(g, genotype_id, genotype_unique_name)

                if genotype_id not in geno_hash:
                    geno_hash[genotype_id] = {}

                genoparts = geno_hash[genotype_id]

                # reassign the allele_type to a proper GENO or SO class
                allele_type = self._map_allele_type_to_geno(allele_type)

                allele_id = 'ZFIN:' + allele_num.strip()

                if allele_num != '':
                    self.id_label_map[allele_id] = allele_name

                # alleles in zfin are really sequence alterations in our system
                geno.addSequenceAlteration(allele_id, allele_name, allele_type)
                gu.addSynonym(g, allele_id, allele_ab)

                # here, we assemble the items into a genotype hash
                # we need to do this because each row only holds a single allele of a gene
                # but a genotype may have many alleles and therefore many rows
                # so we loop through the file once to build a hash of genotype components
                if gene_num is not None and gene_num.strip() != '':
                    # add the gene to the graph, along with it's symbol as the primary label
                    gene_id = 'ZFIN:' + gene_num.strip()
                    geno.addGene(gene_id, gene_symbol)
                    self.id_label_map[gene_id] = gene_symbol

                    # if it's a transgenic construct, then we'll have to get the other bits
                    if construct_num is not None and construct_num.strip() != '':
                        construct_id = 'ZFIN:' + construct_num.strip()
                        geno.addDerivesFrom(allele_id, construct_id)
                        self.id_label_map[construct_id] = construct_name

                    # allele to gene
                    if allele_id not in self.variant_loci_genes:
                        self.variant_loci_genes[allele_id] = [gene_id]
                    else:
                        if gene_id not in self.variant_loci_genes[allele_id]:
                            self.variant_loci_genes[allele_id] += [gene_id]

                    if gene_id not in genoparts:
                        genoparts[gene_id] = [allele_id]
                    else:
                        genoparts[gene_id].append(allele_id)

                    other_allele = self._get_other_allele_by_zygosity(allele_id, zygosity)
                    genoparts[gene_id].append(other_allele)

                    # geno_hash[genotype_id] = genoparts
                else:
                    # if the gene is not known, still need to add the allele to the genotype hash
                    # these will be added as sequence alterations.
                    genoparts[allele_id] = [allele_id]
                    other_allele = self._get_other_allele_by_zygosity(allele_id, zygosity)
                    genoparts[allele_id].append(other_allele)

                    # geno_hash[allele_id] = genoparts    #CHECK should this be genotype_id instead of allele_id?

                geno_hash[genotype_id] = genoparts

                # fetch the other affected genes, and make sure they are in the geno hash
                # we have to do this because some of the affected genes are not listed in this file
                genes_from_hash = None
                if allele_id in self.variant_loci_genes:
                    genes_from_hash = self.variant_loci_genes[allele_id]
                else:
                    logger.info('no gene found for %s', allele_id)

                if genes_from_hash is not None and genes_from_hash != [gene_id] and gene_id not in genes_from_hash:
                    logger.info("***Found other genes not in genotype_features for %s: %s", allele_id, genes_from_hash)
                    for gh in genes_from_hash:
                        if gh not in genoparts:
                            genoparts[gh] = [allele_id]
                        else:
                            genoparts[gh] += [allele_id]

                        other_allele = self._get_other_allele_by_zygosity(allele_id, zygosity)
                        genoparts[gene_id].append(other_allele)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

                # end loop through file
        csvfile.close()


        # ############## BUILD THE INTRINSIC GENOTYPES ###############
        # using the geno_hash, build the genotype parts, and add them to the graph
        # the hash is organized like:
        # genotype_id : { gene_id : [list, of, alleles], # for located things
        #                 allele_id : [list, of, alleles] # for unlocated things
        #               }
        # now loop through the geno_hash, and build the vslcs
        for gt in geno_hash:

            if self.testMode and re.sub('ZFIN:','',gt) not in self.test_ids['genotype']:
                continue

            if gt not in gvc_hash:
                gvc_hash[gt] = {}
            gvcparts = gvc_hash[gt]

            for locus_id in geno_hash.get(gt):
                locus_label = self.id_label_map[locus_id]
                variant_locus_parts = geno_hash.get(gt).get(locus_id)

                # if the locus == part, then it isn't a gene location
                if locus_id in variant_locus_parts:
                    # set the gene_id to none
                    gene_id = None
                else:
                    gene_id = locus_id

                allele1_id = variant_locus_parts[0]
                if allele1_id not in self.id_label_map:
                    allele1_label = allele1_id
                    logger.error('allele1 %s not in hash', allele1_id)
                else:
                    allele1_label = self.id_label_map[allele1_id]
                allele2_id = None
                allele2_label = None
                zygosity_id = None

                if len(variant_locus_parts) > 2:
                    logger.error("There may be a problem. >2 parts for this locus: %s", variant_locus_parts)
                elif len(variant_locus_parts) > 1:
                    allele2_id = variant_locus_parts[1]
                    if allele2_id not in ['0','?']:
                        allele2_label = self.id_label_map[allele2_id]
                    else:
                        allele2_label = allele2_id
                if allele2_id is not None:
                    if allele2_id == '?':
                        zygosity_id = geno.zygosity['indeterminate']
                        allele2_id = None
                    elif allele2_id == '0':
                        zygosity_id = geno.zygosity['hemizygous']
                        allele2_id = None
                    elif allele1_id != allele2_id:
                        zygosity_id = geno.zygosity['complex_heterozygous']
                        # print('heterozygous pair='+allele1_id+'_'+allele2_id)
                    elif allele1_id == allele2_id:
                        zygosity_id = geno.zygosity['homozygous']
                else:
                    zygosity_id = geno.zygosity['simple_heterozygous']

                # make variant_loci
                vloci2 = vloci2_label = None
                if gene_id is not None:
                    vloci1 = self.make_id('-'.join((gene_id, allele1_id)))
                    vloci1_label = geno.make_variant_locus_label(locus_label, allele1_label)
                    geno.addSequenceAlterationToVariantLocus(allele1_id, vloci1)
                    geno.addAlleleOfGene(allele1_id, gene_id, geno.object_properties['has_alternate_part'])
                    gu.addIndividualToGraph(g, vloci1, vloci1_label, geno.genoparts['variant_locus'])
                    if allele2_id is not None:
                        vloci2 = self.make_id('-'.join((gene_id, allele2_id)))
                        vloci2_label = geno.make_variant_locus_label(locus_label, allele2_label)
                        geno.addSequenceAlterationToVariantLocus(allele2_id, vloci2)
                        gu.addIndividualToGraph(g, vloci2, vloci2_label, geno.genoparts['variant_locus'])
                        geno.addAlleleOfGene(allele2_id, gene_id, geno.object_properties['has_alternate_part'])
                else:
                    vloci1 = allele1_id
                    vloci1_label = allele1_label
                    vloci2 = allele2_id
                    vloci2_label = allele2_label

                # create the vslc
                gene_label = ''
                if gene_id is None:
                    gn = ''
                else:
                    gn = gene_id
                    gene_label = self.id_label_map[gene_id]
                if allele2_id is None:
                    a2 = ''
                else:
                    a2 = allele2_id

                # TODO MAKE THIS ANONYMOUS (add a _ prefix)
                # TODO also consider adding this to Genotype.py
                vslc_id = self.make_id('-'.join((gn, allele1_id, a2)))
                vslc_label = geno.make_vslc_label(gene_label, allele1_label, allele2_label)

                #add to global hash
                self.id_label_map[vslc_id] = vslc_label

                gu.addIndividualToGraph(g, vslc_id, vslc_label, geno.genoparts['variant_single_locus_complement'])
                geno.addPartsToVSLC(vslc_id, vloci1, vloci2, zygosity_id,
                                    geno.object_properties['has_alternate_part'],
                                    geno.object_properties['has_alternate_part'])

                if vslc_id not in gvcparts:
                    gvcparts = [vslc_id]
                else:
                    gvcparts += [vslc_id]

                gvc_hash[gt] = gvcparts
                # end loop through geno_hash

        # now loop through the gvc_hash, and build the gvc
        for gt in gvc_hash:
            if self.testMode and re.sub('ZFIN:','',gt) not in self.test_ids['genotype']:
                continue

            gvc_parts = gvc_hash.get(gt)

            # only need to make a gvc specifically if there's >1 vslc
            if len(gvc_parts) > 1:
                gvc_labels = []
                gvc_parts.sort()  # put these in order so they will always make the same id
                gvc_id = self.make_id('-'.join(gvc_parts))
                for vslc_id in gvc_parts:
                    # add the vslc to the gvc
                    geno.addVSLCtoParent(vslc_id, gvc_id)

                    # build the gvc label
                    vslc_label = self.id_label_map[vslc_id]
                    if vslc_label is not None:
                        gvc_labels += [vslc_label]
                    else:
                        gvc_labels += [vslc_id]

                #todo make this anonymous

                #todo sort
                gvc_label = '; '.join(gvc_labels)

                # add the gvc, and add to the genotype
                gu.addIndividualToGraph(g, gvc_id, gvc_label, geno.genoparts['genomic_variation_complement'])
            elif len(gvc_parts) == 1:
                # assign the vslc to be also a gvc
                vslc_id = gvc_parts[0]
                gvc_id = vslc_id
                gvc_label = self.id_label_map[vslc_id]
                gu.addType(g, vslc_id, geno.genoparts['genomic_variation_complement'])
            else:
                logger.error("No GVC parts for %s", gt)

            background_id = self.genotype_backgrounds[genotype_id]
            if background_id is not None:
                background_label = self.id_label_map[background_id]
                if background_label is None:
                    background_label = background_id
                    logger.error("We don't have the label for %s stored", background_id)
            else:
                background_label = '[n.s.]'

            genotype_name = gvc_label + ' [' + background_label + ']'

            geno.addGenotype(gt, genotype_name)

            self.id_label_map[gt] = genotype_name

            # Add the GVC to the genotype
            geno.addParts(gvc_id, gt, geno.object_properties['has_alternate_part'])

            # end of gvc loop

        # end of genotype loop

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
            'unspecified': 'SO:0001059'  # sequence alteration
        }
        if allele_type.strip() in type_map:
            type = type_map.get(allele_type)
        else:
            logger.error("Allele Type (%s) not mapped", allele_type)

        return type

    def _process_genotype_backgrounds(self, limit=None):
        """
        This table provides a mapping of genotypes to their background genotypes.
        Note that the background_id is also a genotype_id.

        Makes these triples:
        <ZFIN:genotype_id> GENO:has_reference_part <ZFIN:background_id>
        <ZFIN:background_id> a GENO:genomic_background
        <ZFIN:background_id> in_taxon <taxon_id>
        <taxon_id> a class
        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        logger.info("Processing genotype backgrounds")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['backgrounds']['file']))
        geno = Genotype(g)

        # Add the taxon as a class
        taxon_id = 'NCBITaxon:7955'  # Danio rerio
        gu.addClassToGraph(g, taxon_id, None)

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (genotype_id, genotype_name, background_id, empty) = row

                if self.testMode and genotype_id not in self.test_ids['genotype']:
                    continue

                genotype_id = 'ZFIN:' + genotype_id.strip()
                background_id = 'ZFIN:' + background_id.strip()

                # store this in the hash for later lookup when building fish genotypes
                self.genotype_backgrounds[genotype_id] = background_id

                # add the background into the graph, in case we haven't seen it before
                geno.addGenomicBackground(background_id, None)

                # hang the taxon from the background
                geno.addTaxon(taxon_id, background_id)

                # add the intrinsic genotype to the graph
                # we DO NOT ADD THE LABEL here as it doesn't include the background
                geno.addGenotype(genotype_id, None, geno.genoparts['intrinsic_genotype'])

                # Add background to the intrinsic genotype
                geno.addGenomicBackgroundToGenotype(background_id, genotype_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with genotype backgrounds")
        return


    def _process_wildtypes(self, limit=None):
        """
        This table provides the genotype IDs, name, and abbreviation of the wildtype genotypes.
        These are the typical genomic backgrounds...there's about 20 of them.

        Triples created:
        <genotype id> a GENO:wildtype
        <genotype id> rdfs:label genotype_abbreviation
        <genotype id> dc:description genotype_name

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        logger.info("Processing wildtype genotypes")
        line_counter = 0
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['wild']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (genotype_id, genotype_name, genotype_abbreviation, empty) = row

                genotype_id = 'ZFIN:' + genotype_id.strip()

                # Add genotype to graph with label and description, as a genomic_background genotype
                geno.addGenotype(genotype_id, genotype_abbreviation,
                                 geno.genoparts['genomic_background'], genotype_name)

                # Build the hash for the wild type genotypes.
                self.wildtype_hash[genotype_id] = genotype_abbreviation
                self.id_label_map[genotype_id] = genotype_abbreviation

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with wildtype genotypes")
        return

    def _process_wildtype_expression(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing wildtype expression")
        line_counter = 0

        raw = '/'.join((self.rawdir, self.files['wild_expression']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (gene_id, gene_symbol, genotype_name, super_structure_id, super_structure_name, sub_structure_id,
                 sub_structure_name, start_stage, end_stage, assay, publication_id, probe_id, antibody_id, empty) = row

                #genotype_id = self.wildtype_hash['id'][genotype_name]
                #genotype_id = 'ZFIN:' + genotype_id.strip()

                # TODO: Consider how to model wildtype genotypes with genes and associated expression.
                gene_id = 'ZFIN:' + gene_id.strip()
                geno.addGene(gene_id, gene_symbol)

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with wildtype expression")
        return

    def _process_stages(self, limit=None):
        """
        This table provides mappings between ZFIN stage IDs and ZFS terms,
        and includes the starting and ending hours for the developmental stage.
        Currently only processing the mapping from the ZFIN stage ID to the ZFS ID.

        Triples created:
        <begin_hour_id> an individual
        <begin_hour_id> rdf:type uo:hours
        <begin_hour_id> rdfs:label values+units

        <end_hour_id> an individual
        <end_hour_id> rdf:type uo:hours
        <end_hour_id> rdfs:label values+units

        <stage_id> an individual
        <stage_id> rdf:type zfs:stage_obo_id (ZFS:1234567)
        <stage_id> rdfs:label values+units

        <stage_id> uberon:existence_starts_at begin_hour_id
        <stage_id> uberon:existence_ends_at end_hour_id
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        logger.info("Processing stages")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['stage']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (stage_id, stage_obo_id, stage_name, begin_hours, end_hours, empty) = row

                # Add the stage as a class, and it's obo equivalent
                stage_id = 'ZFIN:' + stage_id.strip()
                gu.addClassToGraph(g, stage_id, stage_name)
                gu.addEquivalentClass(g, stage_id, stage_obo_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with stages")
        return

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
        gu = GraphUtils(curie_map.get())
        # hardcode
        eco_id = "ECO:0000059"  # experimental_phenotypic_evidence
        raw = '/'.join((self.rawdir, self.files['pheno']['file']))
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
                # FIXME i am only dealing with 'wild-type' environments for now
                #if (not re.match('ZDB-EXP-041102-1', env_id)):
                    #logger.info("Skipping non-wildtype environment %s for %s", env_id, genotype_id)
                    #continue

                genotype_id = 'ZFIN:' + genotype_id.strip()
                env_id = 'ZFIN:' + env_id.strip()
                extrinsic_geno_id = self.make_id(env_id)
                start_stage_id = 'ZFIN:' + start_stage_id.strip()
                end_stage_id = 'ZFIN:' + end_stage_id.strip()

                geno = Genotype(g)
                geno.addGenotype(genotype_id, genotype_name)
                # because we are using only w.t. environments, the genotype is just intrinsic.

                # FIXME: Switch to make_id after QA testing.
                # make an ID for the effective genotype
                # effective_genotype_id = self.make_id(genotype_id+env_id)
                effective_genotype_id = self.make_id(genotype_id+'_'+env_id)

                # FIXME: Need to pass in labels for the intrinsic/extrinsic genotypes to make the effective labels.
                # Only adding the effective genotype if there is an extrinsic genotype present.
                if self.label_hash['genotype_label'].get(extrinsic_geno_id) is not None and self.label_hash['genotype_label'].get(extrinsic_geno_id) != '':
                    intrinsic_genotype_label = self.label_hash['genotype_label'].get(genotype_id)
                    extrinsic_genotype_label = self.label_hash['genotype_label'].get(extrinsic_geno_id)
                    if intrinsic_genotype_label is not None and extrinsic_genotype_label is not None:
                        effective_genotype_label = intrinsic_genotype_label+'; '+extrinsic_genotype_label
                    elif intrinsic_genotype_label is None and extrinsic_genotype_label is not None:
                        effective_genotype_label = extrinsic_genotype_label
                    elif intrinsic_genotype_label is not None and extrinsic_genotype_label is None:
                        effective_genotype_label = intrinsic_genotype_label
                    else:
                        #print(intrinsic_genotype_label)
                        #print(genotype_id)
                        logger.error('No effective genotype label created.')
                        effective_genotype_label = ''
                        #effective_genotype_label = '<empty>'
                    #if intrinsic_genotype_label is not None:
                        #print(intrinsic_genotype_label)
                    #print(effective_genotype_label)
                    geno.addGenotype(effective_genotype_id, effective_genotype_label,
                                     geno.genoparts['effective_genotype'])
                    if extrinsic_genotype_label is not None and extrinsic_genotype_label != '':
                        geno.addParts(extrinsic_geno_id, effective_genotype_id)
                    if intrinsic_genotype_label is not None and intrinsic_genotype_label != '':
                        geno.addParts(genotype_id, effective_genotype_id)

                phenotype_id = self._map_sextuple_to_phenotype(superterm1_id, subterm1_id, quality_id,
                                                               superterm2_id, subterm2_id, modifier)

                if phenotype_id is None:
                    # add this phenotype to a set to report at the end
                    missing_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier])
                    continue
                else:
                    mapped_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier])

                # add abnormal phenotypes
                if not re.match('^normal', modifier):
                    assoc_id = self.make_id((genotype_id+env_id+phenotype_id+pub_id))
                    pub_id = 'ZFIN:' + pub_id.strip()
                    # FIXME: Assuming we change from the intrinsic genotype_id to the effective genotype ID.
                    assoc = G2PAssoc(assoc_id, effective_genotype_id, phenotype_id, pub_id, eco_id)
                    assoc.set_environment(env_id)

                    # TODO look up zfin stage id to ZFS id.
                    if start_stage_id != '':
                        start_stage_id = 'ZFIN:'+start_stage_id
                    if end_stage_id != '':
                        end_stage_id = 'ZFIN:'+end_stage_id
                    assoc.set_stage(start_stage_id, end_stage_id)
                    assoc.addAssociationNodeToGraph(g)

                else:
                    # TODO add normal phenotypes
                    logger.warn("Found normal phenotype; skipping for now")

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        myset = set([','.join(x) for x in mapped_zpids])
        logger.info("%d phenotype-sextuples were mapped", len(myset))
        myset = set([','.join(x) for x in missing_zpids])
        logger.warn("The following %d phenotype-sextuples were not mapped:\n,%s", len(myset), str(myset))

        Assoc().loadAllProperties(g)

        return

    def _process_genes(self, limit=None):
        """
        This table provides the ZFIN gene id, the SO type of the gene, the gene symbol, and the NCBI Gene ID.

        Triples created:
        <gene id> a class
        <gene id> rdfs:label gene_symbol
        <gene id> equivalent class <ncbi_gene_id>
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
        raw = '/'.join((self.rawdir, self.files['gene']['file']))
        geno = Genotype(g)
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, ncbi_gene_id, empty) = row
                
                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue

                gene_id = 'ZFIN:'+gene_id.strip()
                ncbi_gene_id = 'NCBIGene:'+ncbi_gene_id.strip()

                geno.addGene(gene_id, gene_symbol)
                gu.addEquivalentClass(g, gene_id, ncbi_gene_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with genes")
        return

    def _process_features(self, limit=None):
        """
        This module provides information for the intrinsic and extrinsic genotype features of zebrafish.
        All items here are 'alterations', and are therefore instances.

        sequence alteration ID, SO type, abbreviation, and relationship to
        the affected gene, with the gene's ID, symbol, and SO type (gene/pseudogene).

        Triples created:
        <gene id> a class:
        :param limit:
        :return:
        """
        if self.testMode:
            g = self.graph
        else:
            g = self.testgraph

        logger.info("Processing features")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['features']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (genomic_feature_id, feature_so_id, genomic_feature_abbreviation, genomic_feature_name,
                 genomic_feature_type, mutagen, mutagee, construct_id, construct_name, construct_so_id, empty) = row

                genomic_feature_id = 'ZFIN:' + genomic_feature_id.strip()

                gu.addIndividualToGraph(g, genomic_feature_id, genomic_feature_name, feature_so_id)
                gu.addSynonym(g, genomic_feature_id,genomic_feature_abbreviation)
                if construct_id is not None and construct_id != '':
                    construct_id = 'ZFIN:' + construct_id.strip()
                    geno.addConstruct(construct_id, construct_name, construct_so_id)
                    geno.addDerivesFrom(genomic_feature_id, construct_id)

                # Note, we don't really care about how the variant was derived.  so we skip that.

                # add to the id-label map
                self.id_label_map[genomic_feature_id] = genomic_feature_abbreviation
                self.id_label_map[construct_id] = construct_name

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with features")
        return

    def _process_feature_affected_genes(self, limit=None):
        """
        This table lists (intrinsic) genomic sequence alterations and their affected gene(s).
        It provides the sequence alteration ID, SO type, abbreviation, and relationship to
        the affected gene, with the gene's ID, symbol, and SO type (gene/pseudogene).

        Triples created:
        <gene id> a class:
        <gene id> rdfs:label gene_symbol
        <gene id> subclass of gene/pseudogene

        <variant locus id> is a GENO:allele
        <variant locus id> rdfs:label <variant_locus_label>
        <variant locus id> is an allele of <gene id>
        <variant locus id> has alternate part <sequence alteration id>

        <sequence alteration id> is an allele of <gene id>
        <sequence alteration id> rdf:type <sequence alteration type>

        :param limit:
        :return:
        """

        # FIXME THIS NEEDS MAJOR WORK [nlw]
        # can use this to process and build the variant locus.  but will need to process through some kind of
        # helper hash, just like we do with the genotype file.  that's because each gene is on a separate line
        # for example, ZDB-ALT-021021-2 is a deficiency that affects 4 genes
        # that case is when the relationship is != 'is allele of'

        logger.info("Processing feature affected genes")
        line_counter = 0
        raw = '/'.join((self.rawdir, self.files['feature_affected_gene']['file']))
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (genomic_feature_id, feature_so_id, genomic_feature_abbreviation, gene_symbol, gene_id, gene_so_id,
                 genomic_feature_marker_relationship, empty) = row

                # Sequence alteration types present in file: SO:0000159 - deletion, SO:0000199 - translocation,
                # SO:0000667 - insertion, SO:0001059 - sequence_alteration,
                # SO:0001060 - sequence_variant, SO:0001218 - transgenic insertion, SO:1000005 - complex_substitution,
                # SO:1000008 - point_mutation, SO:1000029 - chromosomal_deletion, SO:1000032 - indel

                genomic_feature_id = 'ZFIN:' + genomic_feature_id.strip()
                gene_id = 'ZFIN:' + gene_id.strip()

                self.id_label_map[genomic_feature_id] = genomic_feature_abbreviation
                self.id_label_map[gene_id] = gene_symbol

                geno.addGene(gene_id, gene_symbol, gene_so_id)

                # add the gene to the list of things altered by this thing
                if genomic_feature_id not in self.variant_loci_genes:
                    self.variant_loci_genes[genomic_feature_id] = [gene_id]
                else:
                    if gene_id not in self.variant_loci_genes[genomic_feature_id]:
                        self.variant_loci_genes[genomic_feature_id] += [gene_id]

                sequence_alteration_type = feature_so_id
                # Add the sequence alteration id, label, and type
                geno.addSequenceAlteration(genomic_feature_id, genomic_feature_abbreviation, sequence_alteration_type)

                if sequence_alteration_type == 'is allele of':
                    vl_label = geno.make_variant_locus_label(gene_symbol, genomic_feature_abbreviation)
                    vl_id = self.make_id('-'.join((gene_id, genomic_feature_id)))

                    self.id_label_map[vl_id] = vl_label

                    # create the variant locus, add it's parts and relationship to the gene
                    geno.addSequenceAlterationToVariantLocus(genomic_feature_id, vl_id)
                    geno.addAlleleOfGene(vl_id, gene_id, geno.object_properties['has_alternate_part'])
                    gu.addIndividualToGraph(g, vl_id, vl_label, geno.genoparts['variant_locus'])

                    # note that deficiencies or translocations that affect only one gene are
                    # considered alleles here by zfin, which is appropriate.
                    # i don't yet see duplications
                else:
                    # don't make the variant loci for the other things
                    # which include deficiencies, translocations, transgenes
                    geno.addAlleleOfGene(genomic_feature_id, gene_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with feature affected genes")
        return

    def _process_gene_marker_relationships(self, limit=None):
        """
        Gene-marker relationships include: clone contains gene, coding sequence of, contains polymorphism,
        gene contains small segment, gene encodes small segment, gene has artifact,
        gene hybridized by small segment, gene produces transcript, gene product recognized by antibody,
        knockdown reagent targets gene, promoter of, transcript targets gene

        We only take a fraction of these here... we are interested in the knockdown reagents, promoters, and
        the transgenic constructs with coding bits.

        :param limit:
        :return:
        """
        # TODO need to review relationships with @mbrush
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        logger.info("Processing gene marker relationships")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['gene_marker_rel']['file']))
        geno = Genotype(g)

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol,
                 marker_id, marker_so_id, marker_symbol, relationship, empty) = row

                # Gene-marker relationships: clone contains gene, coding sequence of, contains polymorphism,
                # gene contains small segment, gene encodes small segment, gene has artifact,
                # gene hybridized by small segment, gene produces transcript, gene product recognized by antibody,
                # knockdown reagent targets gene, promoter of, transcript targets gene

                # only process a fraction of these
                if relationship in ['knockdown reagent targets gene', 'coding sequence of',
                                    'gene product recognized by antibody', 'promoter of',
                                    'transcript targets gene']:
                    gene_id = 'ZFIN:'+gene_id.strip()
                    geno.addGene(gene_id, gene_symbol, gene_so_id)

                    marker_id = 'ZFIN:'+marker_id.strip()
                    if relationship == 'knockdown reagent targets gene':
                        geno.addGeneTargetingReagent(marker_id, marker_symbol, marker_so_id)
                        geno.addReagentTargetedGene(marker_id,gene_id)
                        if marker_id not in self.kd_reagent_hash:
                            self.kd_reagent_hash[marker_id] = {'label': marker_symbol,
                                                               'targets': [gene_id]}
                        else:
                            if gene_id not in self.kd_reagent_hash[marker_id]['targets']:
                                self.kd_reagent_hash[marker_id]['targets'] += [gene_id]
                    elif relationship == 'coding sequence of':
                        # transgenic constructs with coding regions
                        # but we don't know if they are wild-type or mutant, so just has_part for now
                        geno.addConstruct(marker_id, marker_symbol, marker_so_id)
                        geno.addParts(gene_id, marker_id, geno.object_properties['has_part'])
                    elif relationship == 'gene product recognized by antibody':
                        # TODO for ticket #32
                        pass
                    elif relationship == 'promoter of':
                        # transgenic constructs with promoters regions
                        # we are making the assumption that they are wild-type promoters
                        geno.addConstruct(marker_id, marker_symbol, marker_so_id)
                        geno.addParts(gene_id, marker_id, geno.object_properties['has_reference_part'])
                    elif relationship == 'transcript targets gene':  #miRNAs
                        # TODO should this be an interaction instead of this special relationship?
                        gu.addIndividualToGraph(g, marker_id, marker_symbol, marker_so_id)
                        gu.addTriple(g, marker_id, geno.object_properties['targets_instance_of'], gene_id)
                    else:
                        pass  #skip

                self.id_label_map[marker_id] = marker_symbol
                self.id_label_map[gene_id] = gene_symbol  # just in case we haven't seen it before

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with gene marker relationships")
        return

    def _process_pubinfo(self, limit=None):
        """
        This will pull the zfin internal publication information, and map them to their equivalent
        pmid, and make labels.

        Triples created:
        <pub_id> is an individual
        <pub_id> rdfs:label <pub_label>
        <pubmed_id> is an individual
        <pubmed_id> rdfs:label <pub_label>

        <pub_id> sameIndividual <pubmed_id>
        :param limit:
        :return:
        """
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['pubs']['file']))
        with open(raw, 'r', encoding="latin-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pub_id, pubmed_id, authors, title, journal, year, vol, pages, empty) = row

                if self.testMode and (pub_id.replace('ZFIN:', '') not in self.test_ids['pub']
                                      or pubmed_id.replace('PMID:', '') not in self.test_ids['pub']):
                    continue

                pub_id = 'ZFIN:'+pub_id.strip()
                pub_label = '; '.join((authors, title, journal, year, vol, pages))
                gu.addIndividualToGraph(g, pub_id, pub_label)

                if pubmed_id != '' and pubmed_id is not None:
                    pubmed_id = 'PMID:'+pubmed_id.strip()
                    gu.addIndividualToGraph(g, pubmed_id, None)
                    gu.addSameIndividual(g, pub_id, pubmed_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_pub2pubmed(self, limit=None):
        """
        This will pull the zfin internal publication to pubmed mappings. Somewhat redundant with the
        process_pubinfo method, but this mapping includes additional mappings.

        <pub_id> is an individual
        <pub_id> rdfs:label <pub_label>
        <pubmed_id> is an individual
        <pubmed_id> rdfs:label <pub_label>

        <pub_id> sameIndividual <pubmed_id>
        :param limit:
        :return:
        """
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['pub2pubmed']['file']))
        with open(raw, 'r', encoding="latin-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pub_id, pubmed_id, empty) = row

                if self.testMode and (pub_id.replace('ZFIN:', '') not in self.test_ids['pub']
                                      or pubmed_id.replace('PMID:', '') not in self.test_ids['pub']):
                    continue

                pub_id = 'ZFIN:'+pub_id.strip()
                gu.addIndividualToGraph(g, pub_id, None)

                if pubmed_id != '' and pubmed_id is not None:
                    pubmed_id = 'PMID:'+pubmed_id.strip()
                    gu.addIndividualToGraph(g, pubmed_id, None)
                    gu.addSameIndividual(g, pub_id, pubmed_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_targeting_reagents(self, reagent_type, limit=None):
        """
        This method processes the gene targeting knockdown reagents, such as morpholinos,
        talens, and crisprs.  We create triples for the reagents and pass the
        data into a hash map for use in the pheno_enviro method.

        Morpholinos work similar to RNAi.
        TALENs are artificial restriction enzymes that can be used for genome editing in situ.
        CRISPRs are knockdown reagents, working similar to RNAi but at the transcriptional level instead of mRNA level.

        You can read more about TALEN and CRISPR techniques in review [Gaj et al](http://www.cell.com/trends/biotechnology/abstract/S0167-7799%2813%2900087-5)

        TODO add sequences

        Triples created:
        <reagent_id> is a gene_targeting_reagent
        <reagent_id> rdfs:label <reagent_symbol>
        <reagent_id> has type <reagent_so_id>
        <reagent_id> has comment <note>

        <publication_id> is an individual
        <publication_id> mentions <morpholino_id>
        :param reagent_type: should be one of: morph, talen, crispr
        :param limit:
        :return:
        """

        logger.info("Processing Gene Targeting Ragents")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        geno = Genotype(g)

        if reagent_type not in ['morph','talen','crispr']:
            logger.error("You didn't specify the right kind of file type.")
            return

        raw = '/'.join((self.rawdir, self.files[reagent_type]['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                if reagent_type in ['morph','crispr']:
                    (gene_id, gene_so_id, gene_symbol, reagent_id, reagent_so_id,
                     reagent_symbol, reagent_sequence, publication, note) = row
                elif reagent_type == 'talen':
                    (gene_id, gene_so_id, gene_symbol, reagent_id, reagent_so_id,
                     reagent_symbol, reagent_sequence, reagent_sequence2,
                     publication, note) = row
                else:
                    # should not get here
                    return

                if self.testMode and reagent_id not in self.test_ids['morpholino']:
                    continue

                reagent_id = 'ZFIN:'+reagent_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()

                # FIXME: This is incorrect, as it requires the concentration if available, and is related to the extrinsic genotype.
                # Should add the morpholino as a typed individual instead. Same for TALENs/CRISPRs.
                geno.addGeneTargetingReagent(reagent_id, reagent_symbol, reagent_so_id)
                # Now adding the reagent targeted gene in the pheno_environment processing function.
                #geno.addReagentTargetedGene(morpholino_id,gene_id, gene_id)

                # Add publication
                if publication != '':
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g, pub_id, None)
                    gu.addTriple(g, pub_id, gu.properties['mentions'], reagent_id)

                # Add comment?
                if note != '':
                    gu.addComment(g, reagent_id, note)

                # Build the hash for the reagents and the gene targets
                if reagent_id not in self.kd_reagent_hash:
                    self.kd_reagent_hash[reagent_id] = {'label': reagent_symbol,
                                                        'targets': [gene_id]}
                else:
                    self.kd_reagent_hash[reagent_id]['targets']+=[gene_id]

                self.id_label_map[reagent_id] = reagent_symbol

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with Reagent type %s", reagent_type)
        return


    def _process_pheno_enviro(self, limit=None):
        """
        The pheno_environment.txt file ties experimental conditions to an environment ID.
        An environment ID may have one or more associated conditions.
        Condition groups present: chemical, CRISPR, morpholino, pH, physical, physiological, salinity, TALEN,
        temperature, and Generic-control.
        The condition column may contain knockdown reagent IDs or mixed text.

        This method also constructs the extrinsic genotype and sub-parts, currently just using the morpholinos.
        TALENs and CRISPRs will be added after modeling is completed.

        Triples created:
        <environment_id> is an Individual
        <environment_id> rdfs:label <environment_label>
        <environment_id> has type <environment>


        <extrinsic_id> is an Individual
        <extrinsic_id> rdfs:label <extrinsic_genotype>
        <extrinsic_id> has type <extrinsic_genotype>

        <targeted_gene_subregion_id> is an Individual
        <targeted_gene_subregion_id> rdfs:label <targeted_gene_subregion_label>
        <targeted_gene_subregion_id> has type <targeted_gene_subregion>
        <targeted_gene_subregion_id> has part <morpholino_id>

        <targeted_gene_variant_id> is an Individual
        <targeted_gene_variant_id> rdfs:label <targeted_gene_variant_label>
        <targeted_gene_variant_id> has type <reagent_targeted_gene>
        <targeted_gene_variant_id> has part <targeted_gene_subregion_id>

        <targeted_gene_complement_id> is an Individual
        <targeted_gene_complement_id> rdfs:label <targeted_gene_complement_label>
        <targeted_gene_complement_id> has type <targeted_gene_complement>
        <targeted_gene_complement_id> has part <targeted_gene_variant_id>

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
        kd_reagent_conc_hash = {}
        kd_reagent_conc_label_hash= {}
        extrinsic_part_hash = {}
        enviro_label_hash = {}
        geno = Genotype(self.graph)
        # condition_
        raw = '/'.join((self.rawdir, self.files['enviro']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (environment_id, condition_group, condition, values, units, comment, empty) = row

                if re.match("\\\\", values):
                    values = ''
                environment_id = 'ZFIN:'+environment_id.strip()
                if self.testMode and environment_id not in self.test_ids['environment']:
                    continue

                # FIXME: Can the environment serve as the extrinsic genotype ID?  ***NO***
                # For now, making an ID from environment ID only, but may want to revisit.
                extrinsic_geno_id = self.make_id(environment_id)
                self.extrinsic_id_to_enviro_id_hash[extrinsic_geno_id] = environment_id

                # We can start to build the extrinsic genotype using this file.
                # Requires creating a hash similar to what is used for genotypes to get the VSLCs and GVCs.
                # TODO: For now adding Morpholinos and not working with temp/chemical/physical/etc, aside from labels
                # FIXME: For now using the general "Environment" geno ID (GENO:0000099)
                # There are a few specific environments available in GENO, including some standard
                # zfin environments (standard salinity and temperature, heat shock (37C), etc), which
                # includes the zfin ID instead of a GENO ID for those environments.

                # Clean up the units
                if units == 'N/A' or units == '':
                    units = None

                # Clean up the values
                if values == '' or values == 'N/A':
                    values = None

                if comment == 'NULL' or comment == '':
                    comment = None

                # Use this regex match if using all knockdown reagents.
                #if re.match('ZDB.*',condition):
                # Use this regex match if using only morpholino knockdown reagents.
                if re.match('ZDB-MRPHLNO.*', condition):

                    condition = 'ZFIN:'+condition.strip()
                    # TODO: In the future may want to subdivide from the general "environment" type
                    # and map to specific environment types (altered chemical environment,
                    # elevated temperature environment, altered light environment, etc.)
                    gu.addIndividualToGraph(self.graph, environment_id, None, gu.datatype_properties['environment'])
                    geno.addGenotype(extrinsic_geno_id, None, geno.genoparts['extrinsic_genotype'])

                    # Clean up the units
                    if units is not None and re.match('.*\/.*', units):
                        units = re.sub(r"/", '_', units)

                    # Clean up the values
                    if values == '' or values == 'N/A':
                        values = None
                    if values is not None:
                        values = values.replace(' ', '_')
                        # FIXME: Better way to indicate > and < ?
                        values = values.replace('<', 'less_than_')
                        values = values.replace('>', 'greater_than_')

                    #if units is not None and values is not None:
                        #print(values+units)

                    # Create the targeted sequence id
                    if units is not None and values is not None:
                        targeted_sequence_id = condition+'_'+values+units
                        conc_label = '('+values+' '+units+')'
                    else:
                        # FIXME: Better way to indicate that the concentration is not provided?
                        targeted_sequence_id = condition+'_ns'
                        conc_label = '(n.s.)'
                    #print(targeted_sequence_id)

                    if extrinsic_geno_id not in extrinsic_part_hash:
                        extrinsic_part_hash[extrinsic_geno_id] = [condition]
                        #extrinsic_parts = extrinsic_geno_hash[extrinsic_geno_id]
                        #print(extrinsic_parts)

                    if condition not in extrinsic_part_hash[extrinsic_geno_id]:
                        extrinsic_part_hash[extrinsic_geno_id].append(condition)

                    if extrinsic_geno_id not in kd_reagent_conc_hash:
                        kd_reagent_conc_hash[extrinsic_geno_id] = {}
                    # TODO:Change to a make_id after testing.
                    targeted_sequence_key = extrinsic_geno_id+condition

                    if condition not in kd_reagent_conc_hash[extrinsic_geno_id]:
                        kd_reagent_conc_hash[extrinsic_geno_id][condition] = targeted_sequence_id

                    if condition not in self.kd_reagent_hash['kd_reagent_label']:
                        kd_reagent_label = ''
                    else:
                        # targeted gene subregion label will come from hash
                        kd_reagent_label = self.kd_reagent_hash['kd_reagent_label'][condition]
                    #print(kd_reagent_label)
                    targeted_gene_subregion_label = '<'+kd_reagent_label+' '+conc_label+'>'
                    #print(targeted_gene_subregion_label)

                    if comment is None or comment == '':
                        environment_label = condition_group+'['+condition+': '+kd_reagent_label+' '+conc_label+']'
                    else:
                        environment_label = condition_group+'['+condition+': '+kd_reagent_label+' '+\
                                            conc_label+' ('+comment+')]'

                    if environment_id not in enviro_label_hash:
                        enviro_label_hash[environment_id] = [environment_label]
                    else:
                        enviro_label_hash[environment_id].append(environment_label)

                    if extrinsic_geno_id not in kd_reagent_conc_label_hash:
                        kd_reagent_conc_label_hash[extrinsic_geno_id] = {}

                    if condition not in kd_reagent_conc_label_hash[extrinsic_geno_id]:
                        kd_reagent_conc_label_hash[extrinsic_geno_id][condition] = targeted_gene_subregion_label

                    #print(kd_reagent_conc_hash[extrinsic_geno_id][condition])
                    #print(kd_reagent_conc_hash[extrinsic_geno_id])

                    #if condition not in extrinsic_part_hash[extrinsic_geno_id]:

                    #gvc_hash[vslc_id] = {};
                    #if vslc_counter == 0:
                        #gvcparts[gt_vslc] = [vslc_id]
                    #elif vslc_id not in gvcparts:
                        #gvcparts[gt_vslc].append(vslc_id)
                    #vslc_counter += 1

                    #if extrinsic_geno_id[morpholino] not in extrinsic_parts:
                            #extrinsic_geno_hash[extrinsic_geno_id][morpholino] = {};
                            #extrinsic_parts = extrinsic_geno_hash[extrinsic_geno_id][morpholino]
                            #extrinsic_parts[enviro_con].append(condition)
                    #except KeyError:
                        #extrinsic_parts[enviro_con] = [condition]


                # Add the TALENs/CRISPRs to the environment but not the extrinsic genotype
                elif re.match('ZDB.*', condition):
                    condition = 'ZFIN:'+condition.strip()
                    gu.addIndividualToGraph(self.graph, environment_id, None, gu.datatype_properties['environment'])

                    # Clean up the units
                    if units is not None and re.match('.*\/.*', units):
                        units = re.sub(r"/", '_', units)

                    # Clean up the values
                    if values == '' or values == 'N/A':
                        values = None
                    if values is not None:
                        values = values.replace(' ', '_')
                        # FIXME: Better way to indicate > and < ?
                        values = values.replace('<', 'less_than_')
                        values = values.replace('>', 'greater_than_')

                    # if units is not None and values is not None:
                        #print(values+units)

                    # Create the targeted sequence id
                    if units is not None and values is not None:
                        targeted_sequence_id = condition+'_'+values+units
                        conc_label = '('+values+' '+units+')'
                    else:
                        # FIXME: Better way to indicate that the concentration is not provided?
                        targeted_sequence_id = condition+'_ns'
                        conc_label = '(n.s.)'
                    #print(targeted_sequence_id

                    # targeted gene subregion label will come from hash
                    kd_reagent_label = self.kd_reagent_hash['kd_reagent_label'][condition]

                    if comment is None or comment == '':
                        environment_label = condition_group+'['+condition+': '+kd_reagent_label+' '+conc_label+']'
                    else:
                        environment_label = condition_group+'['+condition+': '+kd_reagent_label+' '+conc_label+' ('+comment+')]'

                    if environment_id not in enviro_label_hash:
                        enviro_label_hash[environment_id] = [environment_label]
                    else:
                        enviro_label_hash[environment_id].append(environment_label)

                # FIXME: Can remove this if we don't want to deal with any other abnormal environments.
                elif not re.match('ZDB.*', condition):
                    # FIXME:Need to adjust label for non-knockdown reagent environments

                    if values is not None and units is not None and comment is not None:
                        enviro_label = condition_group+'['+condition+': '+values+units+' ('+comment+')]'
                    elif values is not None and units is not None and comment is None:
                        enviro_label = condition_group+'['+condition+': '+values+units+']'
                    elif values is not None and units is None and comment is not None:
                        enviro_label = condition_group+'['+condition+': '+values+' ('+comment+')]'
                    elif values is not None and units is None and comment is None:
                        enviro_label = condition_group+'['+condition+': '+values+']'
                    elif values is None and units is None and comment is None:
                        enviro_label = condition_group+'['+condition+']'
                    elif values is None and units is None and comment is not None:
                        enviro_label = condition_group+'['+condition+' ('+comment+')]'
                    elif values is None and units is not None and comment is None:
                        enviro_label = condition_group+'['+condition+': '+units+']'
                    elif values is None and units is not None and comment is not None:
                        enviro_label = condition_group+'['+condition+': '+units+' ('+comment+')]'
                    else:
                        logger.warn('No environment label created for environment %s.', environment_id)
                        #enviro_label = '<empty>'
                        enviro_label = ''
                    #print(enviro_label)

                    if environment_id not in enviro_label_hash:
                        enviro_label_hash[environment_id] = [enviro_label]
                    else:
                        enviro_label_hash[environment_id].append(enviro_label)

                    # Adding this results in additional environmental variables being added
                    # to the morpholino environment.
                    #gu.addIndividualToGraph(self.graph,environment_id,None,gu.datatype_properties['environment'],condition_group)

                # TODO: Need to wrangle a better description, alternative parsing of variables
                # (condition_group, condition, values, units, comment). Leaving as condition group for now.
                # Data is problematic with differing values (numeric values, N/A's, blanks).
                #if(comment !=''):
                    #enviro_description = condition_group+': '+condition+' at '+values+' '+units+comment
                #else:
                    #enviro_description = condition_group+': '+condition+' at '+values+' '+units
                #print(enviro_description)
                #gu.addIndividualToGraph(g,environment_id,None,gu.datatype_properties['environment'],condition_group)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

                # End of loop

            # Now process through the enviro_label_hash to add labels to the environments.

            for environment_id in enviro_label_hash:
                environment_label = '; '.join(enviro_label_hash[environment_id])
                #print(environment_id+': '+environment_label)
                gu.addIndividualToGraph(self.graph, environment_id, environment_label,
                                        gu.datatype_properties['environment'])

            # Now process through the extrinsic_part_hash to produce targeted_gene_subregion and targeted_gene_variant
            #print(extrinsic_part_hash)
            tgc_hash = {}
            for extrinsic_geno_id in extrinsic_part_hash:
                #print(extrinsic_part_hash[extrinsic_geno_id])

                geno = Genotype(self.graph)
                ex_geno = geno.addGenotype(extrinsic_geno_id, None, geno.genoparts['extrinsic_genotype'])
                for condition in extrinsic_part_hash[extrinsic_geno_id]:
                    kd_reagent_conc_id = kd_reagent_conc_hash[extrinsic_geno_id][condition]
                    if condition not in self.kd_reagent_hash['kd_reagent_id']:
                        logger.error("%s not in kd_reagent_id")
                        kd_reagent_gene_ids = []
                    else:
                        kd_reagent_gene_ids = self.kd_reagent_hash['kd_reagent_id'][condition]
                        #print(kd_reagent_gene_ids)

                    # Make the tgs id and label, add tgs to graph
                    targeted_gene_subregion_label = kd_reagent_conc_label_hash[extrinsic_geno_id][condition]
                    # TODO: Change to makeID after testing.
                    #targeted_gene_subregion_id = kd_reagent_conc_id+ ('_').join(kd_reagent_gene_ids)
                    targeted_gene_subregion_id = self.make_id(kd_reagent_conc_id+ '_'.join(kd_reagent_gene_ids))
                    #print(targeted_gene_subregion_label)
                    geno.addTargetedGeneSubregion(targeted_gene_subregion_id, targeted_gene_subregion_label)
                    geno.addParts(condition, targeted_gene_subregion_id)

                    for i in kd_reagent_gene_ids:
                        # TODO: Change to a makeID after testing.
                        #targeted_gene_variant_id = i+'_'+kd_reagent_conc_id
                        targeted_gene_variant_id = self.make_id(i+'_'+kd_reagent_conc_id)
                        # FIXME: What about for reagents that target more than one gene? Concatenated or separate?
                        #print(targeted_gene_variant_id)
                        kd_reagent_gene_label = self.kd_reagent_hash['gene_label'][i]
                        kd_reagent_conc_label = self.kd_reagent_hash['kd_reagent_id'][condition]
                        targeted_gene_variant_label = kd_reagent_gene_label+targeted_gene_subregion_label
                        #print('tgv_id='+targeted_gene_variant_id)
                        #print('tgv_label='+targeted_gene_variant_label)
                        geno.addReagentTargetedGene(condition, i, targeted_gene_variant_id, targeted_gene_variant_label)
                        geno.addParts(targeted_gene_subregion_id, targeted_gene_variant_id)

                        if extrinsic_geno_id not in tgc_hash:
                            tgc_hash[extrinsic_geno_id] = {}

                        if targeted_gene_variant_id not in tgc_hash[extrinsic_geno_id]:
                            tgc_hash[extrinsic_geno_id][targeted_gene_variant_id] = targeted_gene_variant_label

                # End of loop
            # Now process through the tgc_hash to produce the targeted_gene_variant_complement
            for extrinsic_geno_id in tgc_hash:
                tgc_ids = []
                tgc_labels = []
                geno = Genotype(self.graph)

                for targeted_gene_variant_id in tgc_hash[extrinsic_geno_id]:
                    if targeted_gene_variant_id not in tgc_ids:
                        tgc_ids.append(targeted_gene_variant_id)
                    if tgc_hash[extrinsic_geno_id][targeted_gene_variant_id] not in tgc_labels:
                        tgc_labels.append(tgc_hash[extrinsic_geno_id][targeted_gene_variant_id])
                # FIXME:Change to MakeID after QA testing.
                #targeted_gene_complement_id = ('_').join(tgc_ids)
                targeted_gene_complement_id = self.make_id('_'.join(tgc_ids))
                targeted_gene_complement_label = ('; ').join(tgc_labels)
                # FIXME: For now just using the TGC label as the extrinsic genotype label
                ex_geno = geno.addGenotype(extrinsic_geno_id, targeted_gene_complement_label,
                                           geno.genoparts['extrinsic_genotype'])
                geno.addTargetedGeneComplement(targeted_gene_complement_id, targeted_gene_complement_label)
                if self.label_hash['genotype_label'].get(extrinsic_geno_id) is None:
                    self.label_hash['genotype_label'][extrinsic_geno_id] = targeted_gene_complement_label
                # TODO: Abstract adding TGC to Genotype.
                # Add the TGC to the genotype.
                geno.addParts(targeted_gene_complement_id, extrinsic_geno_id)
                # TODO: Abstract adding TGVs to TGCs.
                for targeted_gene_variant_id in tgc_hash[extrinsic_geno_id]:
                    geno.addParts(targeted_gene_variant_id, targeted_gene_complement_id)

        #print(extrinsic_part_hash)
        logger.info("Done with phenotype environments")

        return


    def _process_mappings(self, limit=None):
        """
        This function imports linkage mappings of various entities to genetic locations in cM or cR.
        Entities include sequence variants, BAC ends, cDNA, ESTs, genes, PAC ends, RAPDs, SNPs, SSLPs, and  STSs.
        Status: NEEDS REVIEW
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
        raw = '/'.join((self.rawdir, self.files['mappings']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (zfin_id, symbol, so_id, panel_symbol, chromosome, location, metric, empty) = row

                if self.testMode and zfin_id not in self.test_ids['gene']+self.test_ids['allele']:
                    continue

                geno = Genotype(g)

                zfin_id = 'ZFIN:'+zfin_id.strip()
                if re.match('ZFIN:ZDB-GENE.*', zfin_id):
                    geno.addGene(zfin_id, symbol)
                    gu.addClassToGraph(g, zfin_id, symbol, so_id)
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
                taxon_id = 'NCBITaxon:7955'
                taxon_num = '7955'
                taxon_label = 'Danio rerio'
                geno.addChromosome(str(chromosome), taxon_id, taxon_label)

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with chromosome mappings")
        return

    def _process_genbank_ids(self, limit=None):
        """
        This file contains BACs, cDNAs, engineered foreign genes, ESTs, engineered plasmids, Fosmids, pseudogenes,
        engineered plasmids, P1 artificial chromosomes, SSLPs, and STS's in addition to genes, maps all to GenBank IDs.

        Triples created:

        :param limit:
        :return:
        """
        # TODO: Test the output, make sure the GenBank URI resolves for all construct types.
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

                geno = Genotype(g)

                zfin_id = 'ZFIN:'+zfin_id.strip()
                genbank_id = 'GenBank:'+genbank_id.strip()
                if re.match('ZFIN:ZDB-GENE.*', zfin_id):
                    geno.addGene(zfin_id, symbol)
                    gu.addClassToGraph(g, genbank_id, symbol, so_id)
                    gu.addEquivalentClass(g, zfin_id, genbank_id)
                else:
                    geno.addConstruct(zfin_id, symbol, so_id)
                    gu.addIndividualToGraph(g, genbank_id, symbol, so_id)
                    gu.addSameIndividual(g, zfin_id, genbank_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with GenBank IDs")
        return
    
    def _process_uniprot_ids(self, limit=None):
        """
        This method processes the the mappings from ZFIN gene IDs to UniProtKB IDs.

        Triples created:
        <zfin_gene_id> a class
        <zfin_gene_id> rdfs:label gene_symbol

        <uniprot_id> is an Individual
        <uniprot_id> has type <polypeptide>

        <zfin_gene_id> has_gene_product <uniprot_id>
        :param limit:
        :return:
        """

        logger.info("Processing UniProt IDs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['uniprot']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, uniprot_id, empty) = row
                
                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue

                geno = Genotype(g)

                gene_id = 'ZFIN:'+gene_id.strip()
                uniprot_id = 'UniProtKB:'+uniprot_id.strip()

                geno.addGene(gene_id, gene_symbol)
                # Need to add some type of 'has_gene_product' relationship here.
                # TODO: Abstract to one of the model utilities
                gu.addIndividualToGraph(g, uniprot_id, None, geno.genoparts['polypeptide'])
                gu.addTriple(g, gene_id, gu.properties['has_gene_product'], uniprot_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break
                    
        logger.info("Done with UniProt IDs")
        return

    def _process_human_orthos(self, limit=None):
        """
        This table provides ortholog mappings between zebrafish and humans, including OMIM IDs

        Triples created:
        <zfin gene id> a class
        <zfin gene id> rdfs:label gene_symbol
        <zfin gene id> dc:description gene_name

        <human gene id> a class
        <human gene id> rdfs:label gene_symbol
        <human gene id> dc:description gene_name
        <human gene id> equivalent class <omim id>

        <zfin gene id> orthology association <human gene id>
        :param limit:
        :return:
        """

        # Is this file necessary if we can get human orthologs through PANTHER?
        # Are the ZFIN genes mapped to an NCBI Gene ID in other files?

        logger.info("Processing human orthos")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['human_orthos']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (zfin_id, zfin_symbol, zfin_name, human_symbol, human_name, omim_id, gene_id, empty) = row

                #genotype_id = 'ZFIN:' + genotype_id.strip()

                # Add the zebrafish gene.
                zfin_id = 'ZFIN:' + zfin_id.strip()
                geno.addGene(zfin_id, zfin_symbol, None, zfin_name)

                # Add the human gene.
                gene_id = 'NCBIGene:' + gene_id.strip()
                geno.addGene(gene_id, human_symbol, None, human_name)

                # TODO: Need to add the ortholog relationship between the zebrafish gene and the human gene
                # Is this the correct handling of the relationship?
                assoc_id = self.make_id(''.join((zfin_id, gene_id)))
                assoc = OrthologyAssoc(assoc_id, zfin_id, gene_id, None, None)
                assoc.setRelationship('RO:HOM0000017')
                #assoc.loadAllProperties(self.graph)    # FIXME inefficient
                assoc.addAssociationToGraph(self.graph)

                # Add the OMIM gene ID as an equivalent class for the human gene.
                omim_id = 'OMIM:' + omim_id.strip()
                gu.addEquivalentClass(self.graph, gene_id, omim_id)

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with human orthos")
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
        # FIXME hardcode
        mod_id = modifier
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
        file = '/'.join((self.rawdir, self.files['zpmap']['file']))
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
        key = self.make_id('_'.join((superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier)))
        return key

    def _get_other_allele_by_zygosity(self, allele_id, zygosity):
        """
        A helper function to switch on the zygosity, and return the appropriate allele id, or symbol.
        :param allele_id:
        :param zygosity:
        :return:
        """
        other_allele = None
        if zygosity == 'homozygous':
            other_allele = allele_id
        elif zygosity == 'hemizygous':
            other_allele = '0'
        elif zygosity == 'unknown':  # we'll use this as a convention
            other_allele = '?'
        elif zygosity == 'complex':  # transgenics
            other_allele = '0'
        elif zygosity == 'heterozygous':
            # passing on hets until a different fxn
            pass
        else:
            logger.warn("Unconfigured zygosity: %s", zygosity)

        return other_allele



    def getTestSuite(self):
        import unittest
        from tests.test_zfin import ZFINTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(ZFINTestCase)

        return test_suite

