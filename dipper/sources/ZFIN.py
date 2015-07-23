import csv
import re
import logging

from dipper.utils import pysed
from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Genotype import Genotype
from dipper.models.assoc.OrthologyAssoc import OrthologyAssoc
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Environment import Environment
from dipper.models.GenomicFeature import makeChromID
from dipper.models.GenomicFeature import Feature
from dipper.models.Reference import Reference
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map


logger = logging.getLogger(__name__)


class ZFIN(Source):
    """
    This is the parser for the [Zebrafish Model Organism Database (ZFIN)](http://www.zfin.org),
    from which we process genotype and phenotype data for laboratory zebrafish.

    We generate the zfin graph to include the following information:
    * genes
    * sequence alterations (includes SNPs/del/ins/indel and large chromosomal rearrangements)
    * transgenic constructs
    * morpholinos, talens, crisprs as expression-affecting reagents
    * genotypes, and their components
    * publications (and their mapping to PMIDs, if available)
    * genotype-to-phenotype associations (including environments and stages at which they are assayed)
    * environmental components
    * orthology to human genes
    * genetic positional information for genes and sequence alterations

    Genotypes leverage the GENO genotype model and includes both intrinsic and extrinsic genotypes.  Where necessary,
    we create anonymous nodes of the genotype partonomy (such as for variant single locus complements,
    genomic variation complements, variant loci, extrinsic genotypes, and extrinsic genotype parts).
    """

    files = {
        'geno': {'file': 'genotype_features.txt', 'url': 'http://zfin.org/downloads/genotype_features.txt'},
        'pheno': {'file': 'phenotype.txt', 'url': 'http://zfin.org/downloads/phenotype.txt'},
        'pubs': {'file': 'zfinpubs.txt', 'url': 'http://zfin.org/downloads/zfinpubs.txt'},
        'zpmap': {'file': 'zp-mapping.txt',
                  'url': 'http://compbio.charite.de/hudson/job/zp-owl/lastSuccessfulBuild/artifact/zp.annot_sourceinfo'},
        # old     'url': 'https://phenotype-ontologies.googlecode.com/svn/trunk/src/ontology/zp/zp-mapping.txt',},
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
                     "ZDB-GENO-001127-3", "ZDB-GENO-001129-1", "ZDB-GENO-090203-8", "ZDB-GENO-070209-1",
                     "ZDB-GENO-070118-1", "ZDB-GENO-140529-1", "ZDB-GENO-070820-1", "ZDB-GENO-071127-3",
                     "ZDB-GENO-000209-20"],
        "gene": ["ZDB-GENE-000616-6", "ZDB-GENE-000710-4", "ZDB-GENE-030131-2773", "ZDB-GENE-030131-8769",
                 "ZDB-GENE-030219-146", "ZDB-GENE-030404-2", "ZDB-GENE-030826-1", "ZDB-GENE-030826-2",
                 "ZDB-GENE-040123-1", "ZDB-GENE-040426-1309", "ZDB-GENE-050522-534", "ZDB-GENE-060503-719",
                 "ZDB-GENE-080405-1", "ZDB-GENE-081211-2", "ZDB-GENE-091118-129", "ZDB-GENE-980526-135",
                 "ZDB-GENE-980526-166", "ZDB-GENE-980526-196", "ZDB-GENE-980526-265", "ZDB-GENE-980526-299",
                 "ZDB-GENE-980526-41", "ZDB-GENE-980526-437", "ZDB-GENE-980526-44", "ZDB-GENE-980526-481",
                 "ZDB-GENE-980526-561", "ZDB-GENE-980526-89", "ZDB-GENE-990415-181", "ZDB-GENE-990415-72",
                 "ZDB-GENE-990415-75", "ZDB-GENE-980526-44", "ZDB-GENE-030421-3", "ZDB-GENE-980526-196",
                 "ZDB-GENE-050320-62", "ZDB-GENE-061013-403", "ZDB-GENE-041114-104", "ZDB-GENE-030131-9700",
                 "ZDB-GENE-031114-1", "ZDB-GENE-990415-72", "ZDB-GENE-030131-2211", "ZDB-GENE-030131-3063",
                 "ZDB-GENE-030131-9460", "ZDB-GENE-980526-26", "ZDB-GENE-980526-27", "ZDB-GENE-980526-29",
                 "ZDB-GENE-071218-6", "ZDB-GENE-070912-423"],
        "allele": ["ZDB-ALT-010426-4", "ZDB-ALT-010427-8", "ZDB-ALT-011017-8", "ZDB-ALT-051005-2", "ZDB-ALT-051227-8",
                   "ZDB-ALT-060221-2", "ZDB-ALT-070314-1", "ZDB-ALT-070409-1", "ZDB-ALT-070420-6", "ZDB-ALT-080528-1",
                   "ZDB-ALT-080528-6", "ZDB-ALT-080827-15", "ZDB-ALT-080908-7", "ZDB-ALT-090316-1", "ZDB-ALT-100519-1",
                   "ZDB-ALT-111024-1", "ZDB-ALT-980203-1374", "ZDB-ALT-980203-412", "ZDB-ALT-980203-465",
                   "ZDB-ALT-980203-470", "ZDB-ALT-980203-605", "ZDB-ALT-980413-636", "ZDB-ALT-021021-2",
                   "ZDB-ALT-080728-1", "ZDB-ALT-100729-1", "ZDB-ALT-980203-1560", "ZDB-ALT-001127-6",
                   "ZDB-ALT-001129-2", "ZDB-ALT-980203-1091", "ZDB-ALT-070118-2", "ZDB-ALT-991005-33",
                   "ZDB-ALT-020918-2", "ZDB-ALT-040913-6"],
        "morpholino": ["ZDB-MRPHLNO-041129-1", "ZDB-MRPHLNO-041129-2", "ZDB-MRPHLNO-041129-3", "ZDB-MRPHLNO-050308-1",
                       "ZDB-MRPHLNO-050308-3", "ZDB-MRPHLNO-060508-2", "ZDB-MRPHLNO-070118-1", "ZDB-MRPHLNO-070522-3",
                       "ZDB-MRPHLNO-070706-1", "ZDB-MRPHLNO-070725-1", "ZDB-MRPHLNO-070725-2", "ZDB-MRPHLNO-071005-1",
                       "ZDB-MRPHLNO-071227-1", "ZDB-MRPHLNO-080307-1", "ZDB-MRPHLNO-080428-1", "ZDB-MRPHLNO-080430-1",
                       "ZDB-MRPHLNO-080919-4", "ZDB-MRPHLNO-081110-3", "ZDB-MRPHLNO-090106-5", "ZDB-MRPHLNO-090114-1",
                       "ZDB-MRPHLNO-090505-1", "ZDB-MRPHLNO-090630-11", "ZDB-MRPHLNO-090804-1", "ZDB-MRPHLNO-100728-1",
                       "ZDB-MRPHLNO-100823-6", "ZDB-MRPHLNO-101105-3", "ZDB-MRPHLNO-110323-3", "ZDB-MRPHLNO-111104-5",
                       "ZDB-MRPHLNO-130222-4", "ZDB-MRPHLNO-080430", "ZDB-MRPHLNO-100823-6", "ZDB-MRPHLNO-140822-1",
                       "ZDB-MRPHLNO-100520-4", "ZDB-MRPHLNO-100520-5", "ZDB-MRPHLNO-100920-3", "ZDB-MRPHLNO-050604-1"],
        "environment": ["ZDB-EXP-050202-1", "ZDB-EXP-071005-3", "ZDB-EXP-071227-14", "ZDB-EXP-080428-1",
                        "ZDB-EXP-080428-2", "ZDB-EXP-080501-1", "ZDB-EXP-080805-7", "ZDB-EXP-080806-5",
                        "ZDB-EXP-080806-8", "ZDB-EXP-080806-9", "ZDB-EXP-081110-3", "ZDB-EXP-090505-2",
                        "ZDB-EXP-100330-7", "ZDB-EXP-100402-1", "ZDB-EXP-100402-2", "ZDB-EXP-100422-3",
                        "ZDB-EXP-100511-5", "ZDB-EXP-101025-12", "ZDB-EXP-101025-13", "ZDB-EXP-110926-4",
                        "ZDB-EXP-110927-1", "ZDB-EXP-120809-5", "ZDB-EXP-120809-7", "ZDB-EXP-120809-9",
                        "ZDB-EXP-120913-5", "ZDB-EXP-130222-13", "ZDB-EXP-130222-7", "ZDB-EXP-130904-2",
                        "ZDB-EXP-041102-1", "ZDB-EXP-140822-13", "ZDB-EXP-041102-1", "ZDB-EXP-070129-3",
                        "ZDB-EXP-110929-7", "ZDB-EXP-100520-2", "ZDB-EXP-100920-3", "ZDB-EXP-100920-5",
                        "ZDB-EXP-090601-2"
                        ],
        "pub": ["PMID:11566854", "PMID:12588855", "PMID:12867027", "PMID:14667409", "PMID:15456722",
                "PMID:16914492", "PMID:17374715", "PMID:17545503", "PMID:17618647", "PMID:17785424",
                "PMID:18201692", "PMID:18358464", "PMID:18388326", "PMID:18638469", "PMID:18846223",
                "PMID:19151781", "PMID:19759004", "PMID:19855021", "PMID:20040115", "PMID:20138861",
                "PMID:20306498", "PMID:20442775", "PMID:20603019", "PMID:21147088", "PMID:21893049",
                "PMID:21925157", "PMID:22718903", "PMID:22814753", "PMID:22960038", "PMID:22996643",
                "PMID:23086717", "PMID:23203810", "PMID:23760954", "ZFIN:ZDB-PUB-140303-33",
                "ZFIN:ZDB-PUB-140404-9", "ZFIN:ZDB-PUB-080902-16", "ZFIN:ZDB-PUB-101222-7", "ZFIN:ZDB-PUB-140614-2",
                "ZFIN:ZDB-PUB-120927-26", "ZFIN:ZDB-PUB-100504-5"]
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

        # pubs has control characters!
        self._remove_backslash_r('/'.join((self.rawdir, self.files['pubs']['file'])), 'latin-1')

        return

    def parse(self, limit=None):
        if limit is not None:
            logger.info("Only parsing first %s rows of each file", limit)
        logger.info("Parsing files...")

        self._load_zp_mappings()

        if self.testOnly:
            self.testMode = True

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        self.id_label_map = {}  # to hold any label for a given id
        self.genotype_backgrounds = {}  # to hold the mappings between genotype and background
        self.extrinsic_id_to_enviro_id_hash = {}
        self.variant_loci_genes = {}  # to hold the genes variant due to a seq alt
        self.environment_hash = {}  # to hold the parts of an environment
        self.wildtype_genotypes = []

        # basic information on classes and instances
        self._process_genes(limit)
        self._process_stages(limit)
        self._process_pubinfo(limit)
        self._process_pub2pubmed(limit)

        # The knockdown reagents
        for t in ['morph', 'crispr', 'talen']:
            self._process_targeting_reagents(t, limit)

        self._process_gene_marker_relationships(limit)
        self._process_features(limit)
        self._process_feature_affected_genes(limit)
        self._process_mappings(limit)  # only adds features on chromosomes, not positions

        # These must be processed before G2P and expression
        self._process_wildtypes(limit)
        self._process_genotype_backgrounds(limit)
        self._process_genotype_features(limit)  # REVIEWED - NEED TO REVIEW LABELS ON Deficiencies

        self._process_pheno_enviro(limit)  # Must be processed after morpholinos/talens/crisprs id/label

        # once the genotypes and environments are processed, we can associate these with the phenotypes
        self._process_g2p(limit)

        # zfin-curated orthology calls to human genes
        self._process_human_orthos(limit)

        # FOR THE FUTURE - needs verification
        # self._process_wildtype_expression(limit)
        # self._process_uniprot_ids(limit)

        logger.info("Finished parsing.")

        self.load_bindings()
        gu = GraphUtils(curie_map.get())
        gu.loadAllProperties(g)
        gu.loadObjectProperties(g, Genotype.object_properties)

        logger.info("Found %d nodes in graph", len(self.graph))
        logger.info("Found %d nodes in testgraph", len(self.testgraph))

        return

    def _process_genotype_features(self, limit=None):
        """
        Here we process the genotype_features file, which lists genotypes together with
        any intrinsic sequence alterations, their zygosity, and affected gene.
        Because we don't necessarily get allele pair (VSLC) ids in a single row,
        we therefore iterate through the file and build up a hash that contains all
        of a genotype's partonomy.
        We then assemble a genotype based on that partonomy.
        This table does not list the background genotype/strain: that is listed elsewhere.

        ZFIN "ALT" objects are mapped to sequence alterations in our system.

        By the end of this method, we have built up the intrinsic genotype, with Monarch-style labels.
        All ZFIN labels are added as synonyms (including the "sup" html tags).

        We make assumptions here that any variants that affect the same locus are in trans.  All of the
        genotype parts are created as BNodes at the moment, to avoid minting new Monarch ids, which means
        for anyone consuming this data they are inherently unstable.  This may change in the future.


        :param limit:
        :return:
        """
        raw = '/'.join((self.rawdir, self.files['geno']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        gu = GraphUtils(curie_map.get())
        taxon_id = 'NCBITaxon:7955'  # hardcode to zebrafish

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
                        genoparts[gene_id] += [allele_id]

                    other_allele = self._get_other_allele_by_zygosity(allele_id, zygosity)
                    if other_allele is not None:
                        genoparts[gene_id] += [other_allele]

                else:
                    # if the gene is not known, still need to add the allele to the genotype hash
                    # these will be added as sequence alterations.
                    genoparts[allele_id] = [allele_id]
                    other_allele = self._get_other_allele_by_zygosity(allele_id, zygosity)
                    if other_allele is not None:
                        genoparts[allele_id] += [other_allele]

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
                        if other_allele is not None:
                            genoparts[gh].append(other_allele)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

                    # end loop through file
        csvfile.close()
        logger.info("Finished parsing file")

        # ############## BUILD THE INTRINSIC GENOTYPES ###############
        # using the geno_hash, build the genotype parts, and add them to the graph
        # the hash is organized like:
        # genotype_id : { gene_id : [list, of, alleles], # for located things
        # allele_id : [list, of, alleles] # for unlocated things
        #               }
        # now loop through the geno_hash, and build the vslcs

        logger.info("Building intrinsic genotypes from partonomy")
        for gt in geno_hash:
            if self.testMode and re.sub('ZFIN:', '', gt) not in self.test_ids['genotype']:
                print('skipping ', gt)
                continue

            if gt not in gvc_hash:
                gvc_hash[gt] = []
            gvcparts = gvc_hash[gt]

            for locus_id in geno_hash[gt]:
                # logger.info("locus id %s",locus_id)
                locus_label = self.id_label_map[locus_id]
                variant_locus_parts = geno_hash.get(gt).get(locus_id)
                # logger.info('vl parts: %s',pprint.pformat(variant_locus_parts))
                # if the locus == part, then it isn't a gene, rather a variant not in a specific gene
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
                    logger.error("There may be a problem. >2 parts for this locus (%s): %s", locus_id,
                                 variant_locus_parts)
                elif len(variant_locus_parts) > 1:
                    allele2_id = variant_locus_parts[1]
                    if allele2_id not in ['0', '?']:
                        allele2_label = self.id_label_map[allele2_id]
                    else:
                        allele2_label = allele2_id
                if allele2_id is not None:
                    if allele2_id == '?':
                        zygosity_id = geno.zygosity['indeterminate']
                        allele2_id = 'UN'
                    elif allele2_id == '0':
                        zygosity_id = geno.zygosity['hemizygous']
                    elif allele1_id != allele2_id:
                        zygosity_id = geno.zygosity['complex_heterozygous']
                    elif allele1_id == allele2_id:
                        zygosity_id = geno.zygosity['homozygous']
                else:
                    zygosity_id = geno.zygosity['simple_heterozygous']
                    allele2_label = '+'
                    allele2_id = 'WT'

                # make variant_loci
                vloci2 = vloci2_label = None
                if gene_id is not None:
                    vloci1 = self._make_variant_locus_id(gene_id, allele1_id)
                    vloci1_label = geno.make_variant_locus_label(locus_label, allele1_label)
                    geno.addSequenceAlterationToVariantLocus(allele1_id, vloci1)
                    geno.addAlleleOfGene(vloci1, gene_id)
                    gu.addIndividualToGraph(g, vloci1, vloci1_label, geno.genoparts['variant_locus'])
                    if allele2_id is not None and allele2_id not in ['WT', '0', 'UN']:
                        vloci2 = self._make_variant_locus_id(gene_id, allele2_id)
                        vloci2_label = geno.make_variant_locus_label(locus_label, allele2_label)
                        geno.addSequenceAlterationToVariantLocus(allele2_id, vloci2)
                        gu.addIndividualToGraph(g, vloci2, vloci2_label, geno.genoparts['variant_locus'])
                        geno.addAlleleOfGene(vloci2, gene_id)
                else:
                    vloci1 = allele1_id
                    vloci1_label = allele1_label
                    vloci2 = None
                    if allele2_id not in ['WT', '0', 'UN']:
                        vloci2 = allele2_id
                    vloci2_label = allele2_label

                # create the vslc
                gene_label = ''
                if gene_id is None:
                    gn = 'UN'
                else:
                    gn = gene_id
                    gene_label = self.id_label_map[gene_id]

                # TODO also consider adding this to Genotype.py
                vslc_id = '-'.join((gn, allele1_id, allele2_id))
                vslc_id = '_' + re.sub('(ZFIN)?:', '', vslc_id)

                if self.nobnodes is True:
                    vslc_id = ':' + vslc_id
                vslc_label = geno.make_vslc_label(gene_label, allele1_label, allele2_label)

                # add to global hash
                self.id_label_map[vslc_id] = vslc_label

                gu.addIndividualToGraph(g, vslc_id, vslc_label, geno.genoparts['variant_single_locus_complement'])
                geno.addPartsToVSLC(vslc_id, vloci1, vloci2, zygosity_id,
                                    geno.object_properties['has_alternate_part'],
                                    geno.object_properties['has_alternate_part'])

                gvcparts += [vslc_id]

            gvc_hash[gt] = gvcparts

        # end loop through geno_hash

        logger.info('Finished finding all the intrinsic genotype parts')
        logger.info('Build pretty genotype labels')
        # now loop through the gvc_hash, and build the gvc
        for gt in gvc_hash:
            if self.testMode and re.sub('ZFIN:', '', gt) not in self.test_ids['genotype']:
                continue

            gvc_parts = gvc_hash[gt]

            # only need to make a gvc specifically if there's >1 vslc
            if len(gvc_parts) > 1:
                gvc_labels = []
                gvc_parts.sort()  # put these in order so they will always make the same id

                gvc_id = '-'.join(gvc_parts)
                gvc_id = re.sub('(ZFIN)?:', '', gvc_id)
                gvc_id = '_' + gvc_id
                if self.nobnodes:
                    gvc_id = ':' + gvc_id

                for vslc_id in gvc_parts:
                    # add the vslc to the gvc
                    geno.addVSLCtoParent(vslc_id, gvc_id)

                    # build the gvc label
                    vslc_label = self.id_label_map[vslc_id]
                    if vslc_label is not None:
                        gvc_labels += [vslc_label]
                    else:
                        gvc_labels += [vslc_id]

                gvc_labels.sort()
                gvc_label = '; '.join(gvc_labels)

                # add the gvc to the id-label hash
                self.id_label_map[gvc_id] = gvc_label

                # add the gvc
                gu.addIndividualToGraph(g, gvc_id, gvc_label, geno.genoparts['genomic_variation_complement'])
            elif len(gvc_parts) == 1:
                # assign the vslc to be also a gvc
                vslc_id = gvc_parts[0]
                gvc_id = vslc_id
                gvc_label = self.id_label_map[vslc_id]
                gu.addType(g, vslc_id, geno.genoparts['genomic_variation_complement'])
            else:
                gvc_id = None
                gvc_label = ''
                logger.error("No GVC parts for %s", gt)

            if gt in self.genotype_backgrounds:
                background_id = self.genotype_backgrounds[gt]
                if background_id in self.id_label_map:
                    background_label = self.id_label_map[background_id]
                else:
                    background_label = background_id
                    logger.error("We don't have the label for %s stored", background_id)

                geno.addTaxon(taxon_id, background_id)  # this should already have been done
            else:
                background_id = None
                background_label = 'n.s.'
                # there is no background for this genotype; need to add the taxon to this one!
                geno.addTaxon(taxon_id, gt)

            genotype_name = gvc_label + ' [' + background_label + ']'

            geno.addGenotype(gt, genotype_name)

            self.id_label_map[gt] = genotype_name

            # Add the GVC to the genotype
            geno.addParts(gvc_id, gt, geno.object_properties['has_alternate_part'])

            # end of gvc loop

        # end of genotype loop

        # TODO this is almost complete; deficiencies with >1 locus deleted are still not right
        logger.info("Finished building genotype labels")

        logger.info("Done with genotypes")

        return

    def _map_allele_type_to_geno(self, allele_type):
        t = 'SO:0001059'  # default: sequence_alteration
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
            t = type_map.get(allele_type)
        else:
            logger.error("Allele Type (%s) not mapped", allele_type)

        return t

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
                self.id_label_map[genotype_id] = genotype_abbreviation

                # store these in a special hash to look up later
                self.wildtype_genotypes += [genotype_id]

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with wildtype genotypes")
        return

    # def _process_wildtype_expression(self, limit=None):
    #     """
    #
    #     :param limit:
    #     :return:
    #     """
    #
    #     logger.info("Processing wildtype expression")
    #     line_counter = 0
    #
    #     raw = '/'.join((self.rawdir, self.files['wild_expression']['file']))
    #     with open(raw, 'r', encoding="iso-8859-1") as csvfile:
    #         filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
    #         for row in filereader:
    #             line_counter += 1
    #             geno = Genotype(self.graph)
    #             (gene_id, gene_symbol, genotype_name, super_structure_id, super_structure_name, sub_structure_id,
    #              sub_structure_name, start_stage, end_stage, assay,
    #              publication_id, probe_id, antibody_id, empty) = row
    #
    #             # TODO: Consider how to model wildtype genotypes with genes and associated expression.
    #             gene_id = 'ZFIN:' + gene_id.strip()
    #             geno.addGene(gene_id, gene_symbol)
    #
    #             if limit is not None and line_counter > limit:
    #                 break
    #
    #     logger.info("Done with wildtype expression")
    #     return

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

        # FIXME are these equivalences already in the ZFS?  if so remove this
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
        Here, we process the fish-to-phenotype associations, which also include environmental perturbations.
        The phenotypes may also be recorded as observed at specific stages.
        We create association objects with as much of the information as possible.

        ZFIN combines morpholinos into an "environment" under a single identifier, but this doesn't exactly fit
        into the Monarch genotype model.  We retain the morpholinos within the environmental cocktail,
        but we also create extrinsic genotypes which contain the genome/transcriptome-affecting bits.

        For example, if an "environment" maps to a morpholino plus extra salt, we will both attach
        the environment to the association, AND we will create an extrinsic genotype that contains the
        "gene targeting reagent" at the specified concentration, linked to the gene that it is affecting.

        We now output all associations, regardless of the environmental contribution.  (For a time, we used
        to only output "wildtype" environments.)

        When a fish is created and the phenotype observed, it might have the following as the "subject" of the
        association:

        1.  intrinsic genotype (this could be just a standard wild-type strain like AB, or it might
               have genomic variation)
        2.  intrinsic genotype + extrinsic genotype ([1] plus it includes expression-affecting reagents,
               like morpholinos)
        3.  intrinsic genotype + environment ([1] where the environment is currently a bag of perturbations, including
                temp, salinity, chemicals)
        4.  intrinsic + extrinsic + environment (as in [2-3], including morpholinos/talens/crisprs)

        Any time there is a combination of intrinsic + extrinsic genotypes, this nets an "effective" genotype,
        which contains the intrinsic and extrinsic genotypes as parts.  The assembly of the intrinsic
        and extrinsic genotype parts is done in other methods (therefore this method must be called after
        most others).  For now, we mint new identifiers for the effective genotypes, though eventually these
        will be overwritten by identifiers generated from ZFIN.

        A full association object may include:

        assoc hasSubject effective genotype
        assoc hasObject phenotype
        assoc hasSource publication (PMID or ZFIN pub id)
        assoc hasEnvironment environment
        assoc hasStartStage start_stage_id
        assoc hasEndStage end_stage_id

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

        geno = Genotype(g)
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

                if self.testMode and (genotype_id not in self.test_ids['genotype']
                                      or env_id not in self.test_ids['environment']):
                    continue

                # deal with environments

                genotype_id = 'ZFIN:' + genotype_id.strip()
                env_id = 'ZFIN:' + env_id.strip()

                # the intrinsic genotype
                intrinsic_genotype_label = self.id_label_map[genotype_id]
                if genotype_id not in self.wildtype_genotypes:
                    geno.addGenotype(genotype_id, intrinsic_genotype_label)
                else:
                    pass
                    # wildtype (genomic backgrounds) have already been added

                # get the extrinsic genotype
                # find the extrinsic parts in the hash, and add them here
                extrinsic_genotype_id = None
                extrinsic_genotype_label = ''
                if env_id in self.extrinsic_id_to_enviro_id_hash:
                    extrinsic_genotype_id = self.extrinsic_id_to_enviro_id_hash[env_id]
                    if extrinsic_genotype_id in self.id_label_map:
                        extrinsic_genotype_label = self.id_label_map[extrinsic_genotype_id]

                # build an effective genotype, if relevant
                effective_genotype_id = genotype_id
                effective_genotype_label = intrinsic_genotype_label
                if extrinsic_genotype_id is not None:
                    effective_genotype_id = self._make_effective_genotype_id(genotype_id, extrinsic_genotype_id)

                    effective_genotype_label = '; '.join((intrinsic_genotype_label, extrinsic_genotype_label))
                    self.id_label_map[effective_genotype_id] = effective_genotype_label

                    if genotype_id in self.wildtype_genotypes:
                        intrinsic_part_rel = geno.object_properties['has_reference_part']
                    else:
                        intrinsic_part_rel = geno.object_properties['has_alternate_part']
                    geno.addParts(genotype_id, effective_genotype_id, intrinsic_part_rel)
                    geno.addParts(extrinsic_genotype_id, effective_genotype_id,
                                  geno.object_properties['has_alternate_part'])

                    geno.addGenotype(effective_genotype_id, effective_genotype_label,
                                     geno.genoparts['effective_genotype'])

                # logger.debug("added: %s", effective_genotype_label)

                # ########### PHENOTYPES ##########
                phenotype_id = self._map_sextuple_to_phenotype(superterm1_id, subterm1_id, quality_id,
                                                               superterm2_id, subterm2_id, modifier)

                if phenotype_id is None:
                    # check to see if it is a "normal" phenotype;
                    # if so, then check to see if the "abnormal" version is found
                    # if the abnormal version is not found, then report it
                    if modifier == 'normal':
                        p2 = self._map_sextuple_to_phenotype(superterm1_id, subterm1_id, quality_id,
                                                             superterm2_id, subterm2_id, 'abnormal')
                        if p2 is None:
                            missing_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id,
                                                  subterm2_id, modifier])
                        else:
                            pass
                            # logger.info("Normal phenotype found, and abnormal version exists")
                    else:
                        missing_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id,
                                              subterm2_id, modifier])

                else:
                    mapped_zpids.append([superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, modifier])

                # add abnormal phenotypes
                # logger.info("modifier for %s = %s",effective_genotype_label,modifier)
                if not re.match('^normal', modifier):
                    if phenotype_id is None:
                        continue
                    if start_stage_id != '':
                        start_stage_id = 'ZFIN:' + start_stage_id.strip()
                    if end_stage_id != '':
                        end_stage_id = 'ZFIN:' + end_stage_id.strip()
                    if pub_id != '':
                        pub_id = 'ZFIN:' + pub_id.strip()

                    assoc = G2PAssoc(self.name, effective_genotype_id, phenotype_id)
                    assoc.set_environment(env_id)
                    assoc.set_stage(start_stage_id, end_stage_id)
                    assoc.add_evidence(eco_id)
                    assoc.add_source(pub_id)
                    assoc.add_association_to_graph(g)
                    assoc.load_all_properties(g)

                else:
                    # TODO add normal phenotypes as associations #134
                    # when https://github.com/sba1/bio-ontology-zp/issues/9 is finished, we can use these
                    # add normal phenotypes as a comment on the genotype for now
                    clist = []
                    for x in [superterm1_name, subterm1_name, quality_name, superterm2_name, subterm2_name, modifier]:
                        if x != '':
                            clist += [x]

                    c = '+'.join(clist)
                    c = ' '.join(("Normal phenotype observed:", c, "(" + pub_id + ")"))
                    gu.addComment(g, effective_genotype_id, c)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        myset = set([','.join(x) for x in mapped_zpids])
        myset2 = set([','.join(x) for x in missing_zpids])
        logger.info("Phenotype-sextuples: %d mapped : %d unmapped", len(myset), len(myset2))
        self._write_missing_zp_report(missing_zpids)

        return

    def _write_missing_zp_report(self, missing_zpids, include_normal=True):
        """
        This will write the sextuples of anatomy+quality to a file if they do not map to any current
        ZP definition.
        Set include_normal to False if you do not want to log the unmatched "normal" phenotypes.
        :param missing_zpids:
        :param include_normal:
        :return:
        """

        f = '/'.join((self.outdir, 'missing_zps.txt'))

        myset = set([','.join(x) for x in missing_zpids])
        # missing_zpids = set(missing_zpids)  # make it a unique set
        with open(f, 'w', newline='\n') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            # write header
            h = ['superterm1_id', 'subterm1_id', 'quality_id', 'superterm2_id', 'subterm2_id', 'modifier']
            writer.writerow(h)
            for x in myset:
                writer.writerow(x.split(','))
        csvfile.close()

        logger.info("Wrote %d missing zp defs to %s", len(myset), f)

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

                gene_id = 'ZFIN:' + gene_id.strip()
                ncbi_gene_id = 'NCBIGene:' + ncbi_gene_id.strip()

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
            g = self.testgraph
        else:
            g = self.graph

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

                if self.testMode and (genomic_feature_id not in self.test_ids['allele']):
                    continue

                genomic_feature_id = 'ZFIN:' + genomic_feature_id.strip()
                gu.addIndividualToGraph(g, genomic_feature_id, genomic_feature_name, feature_so_id)

                gu.addSynonym(g, genomic_feature_id, genomic_feature_abbreviation)
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

                if self.testMode and (re.sub('ZFIN:', '', gene_id) not in self.test_ids['gene'] and
                                      re.sub('ZFIN:', '', genomic_feature_id) not in self.test_ids['allele']):
                    continue

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

                    vl_id = self._make_variant_locus_id(gene_id, genomic_feature_id)

                    self.id_label_map[vl_id] = vl_label

                    # create the variant locus, add it's parts and relationship to the gene
                    geno.addSequenceAlterationToVariantLocus(genomic_feature_id, vl_id)
                    gu.addIndividualToGraph(g, vl_id, vl_label, geno.genoparts['variant_locus'])
                    geno.addAlleleOfGene(vl_id, gene_id)

                    # note that deficiencies or translocations that affect only one gene are
                    # considered alleles here by zfin, which is appropriate.
                    # i don't yet see duplications
                else:
                    # don't make the variant loci for the other things
                    # which include deficiencies, translocations, transgenes
                    # TODO review this
                    pass

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

        Here, we only process the following:
        knockdown reagent targets gene, coding sequence of, promoter of, transcript targets gene

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

                if self.testMode and ('ZFIN:' + gene_id not in self.test_ids['gene']
                                      and 'ZFIN:' + marker_id not in self.test_ids['allele']
                                      and 'ZFIN:' + marker_id not in self.test_ids['morpholino']):
                    continue

                # there are many relationships, but we only take a few for now
                if relationship in ['knockdown reagent targets gene', 'coding sequence of',
                                    'gene product recognized by antibody', 'promoter of',
                                    'transcript targets gene']:
                    gene_id = 'ZFIN:' + gene_id.strip()
                    geno.addGene(gene_id, gene_symbol, gene_so_id)

                    marker_id = 'ZFIN:' + marker_id.strip()
                    if relationship == 'knockdown reagent targets gene':
                        geno.addGeneTargetingReagent(marker_id, marker_symbol, marker_so_id, gene_id)
                        # waiting to add the reagent_targeted_gene until processing environments

                    elif relationship == 'coding sequence of':
                        # TODO we know that these need to be modeled as expression variants!
                        # see https://github.com/monarch-initiative/GENO-ontology/issues/20
                        # transgenic constructs with coding regions
                        # but we don't know if they are wild-type or mutant, so just has_part for now
                        geno.addConstruct(marker_id, marker_symbol, marker_so_id)
                        # add a trangene part - TODO move to Genotype.py
                        transgene_part_id = '_' + '-'.join((marker_id, gene_id, re.sub('\W+', '-', relationship)))
                        gu.addIndividualToGraph(g, transgene_part_id, 'Tg(' + relationship + ' ' + gene_symbol + ')',
                                                geno.genoparts['coding_transgene_feature'])
                        geno.addParts(transgene_part_id, marker_id, geno.object_properties['has_alternate_part'])
                        gu.addTriple(g, transgene_part_id, geno.object_properties['derives_sequence_from_gene'],
                                     gene_id)
                    elif relationship == 'gene product recognized by antibody':
                        # TODO for ticket #32
                        pass
                    elif relationship == 'promoter of':
                        # transgenic constructs with promoters regions
                        # we are making the assumption that they are wild-type promoters
                        geno.addConstruct(marker_id, marker_symbol, marker_so_id)

                        # add a trangene part
                        transgene_part_id = '_' + '-'.join((marker_id, gene_id, re.sub('\W+', '-', relationship)))
                        gu.addIndividualToGraph(g, transgene_part_id, 'Tg(' + relationship + ' ' + gene_symbol + ')',
                                                geno.genoparts['regulatory_transgene_feature'])
                        geno.addParts(transgene_part_id, marker_id)
                        gu.addTriple(g, transgene_part_id, geno.object_properties['derives_sequence_from_gene'],
                                     gene_id)
                    elif relationship == 'transcript targets gene':  # miRNAs
                        # TODO should this be an interaction instead of this special relationship?
                        gu.addIndividualToGraph(g, marker_id, marker_symbol, marker_so_id)
                        gu.addTriple(g, marker_id, geno.object_properties['targets_instance_of'], gene_id)
                    else:
                        pass

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

                if self.testMode and ('ZFIN:' + pub_id not in self.test_ids['pub']
                                      and 'PMID:' + pubmed_id not in self.test_ids['pub']):
                    continue

                pub_id = 'ZFIN:' + pub_id.strip()
                # trim the author list for ease of reading
                alist = re.split(',', authors)
                if len(alist) > 1:
                    astring = ' '.join((alist[0].strip(), 'et al'))
                else:
                    astring = authors

                pub_label = '; '.join((astring, title, journal, year, vol, pages))
                r = Reference(pub_id)
                r.setShortCitation(pub_label)
                r.setYear(year)
                r.setTitle(title)

                if pubmed_id is not None and pubmed_id != '':
                    # let's make an assumption that if there's a pubmed id, that it is a journal article
                    r.setType(Reference.ref_types['journal_article'])

                    pubmed_id = 'PMID:' + pubmed_id.strip()
                    rpm = Reference(pubmed_id, Reference.ref_types['journal_article'])
                    rpm.addRefToGraph(g)

                    gu.addSameIndividual(g, pub_id, pubmed_id)

                r.addRefToGraph(g)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_pub2pubmed(self, limit=None):
        """
        This will pull the zfin internal publication to pubmed mappings. Somewhat redundant with the
        process_pubinfo method, but this includes additional mappings.

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

                if self.testMode and ('ZFIN:' + pub_id not in self.test_ids['pub']
                                      and 'PMID:' + pubmed_id not in self.test_ids['pub']):
                    continue

                pub_id = 'ZFIN:' + pub_id.strip()
                gu.addIndividualToGraph(g, pub_id, None)

                if pubmed_id != '' and pubmed_id is not None:
                    pubmed_id = 'PMID:' + pubmed_id.strip()
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

        logger.info("Processing Gene Targeting Reagents")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        geno = Genotype(g)

        if reagent_type not in ['morph', 'talen', 'crispr']:
            logger.error("You didn't specify the right kind of file type.")
            return

        raw = '/'.join((self.rawdir, self.files[reagent_type]['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                if reagent_type in ['morph', 'crispr']:
                    (gene_num, gene_so_id, gene_symbol, reagent_num, reagent_so_id,
                     reagent_symbol, reagent_sequence, publication, note) = row
                elif reagent_type == 'talen':
                    (gene_num, gene_so_id, gene_symbol, reagent_num, reagent_so_id,
                     reagent_symbol, reagent_sequence, reagent_sequence2,
                     publication, note) = row
                else:
                    # should not get here
                    return

                reagent_id = 'ZFIN:' + reagent_num.strip()
                gene_id = 'ZFIN:' + gene_num.strip()

                self.id_label_map[reagent_id] = reagent_symbol

                if self.testMode and (reagent_num not in self.test_ids['morpholino']
                                      and gene_num not in self.test_ids['gene']):
                    continue

                geno.addGeneTargetingReagent(reagent_id, reagent_symbol, reagent_so_id, gene_id)
                # The reagent targeted gene is added in the pheno_environment processing function.

                # Add publication
                # note that the publications can be comma-delimited, like: ZDB-PUB-100719-4,ZDB-PUB-130703-22
                if publication != '':
                    pubs = re.split(',', publication.strip())
                    for p in pubs:
                        pub_id = 'ZFIN:' + p.strip()
                        gu.addIndividualToGraph(g, pub_id, None)
                        gu.addTriple(g, pub_id, gu.properties['mentions'], reagent_id)

                # Add comment?
                if note != '':
                    gu.addComment(g, reagent_id, note)

                # use the variant hash for reagents to list the affected genes
                if reagent_id not in self.variant_loci_genes:
                    self.variant_loci_genes[reagent_id] = [gene_id]
                else:
                    if gene_id not in self.variant_loci_genes[reagent_id]:
                        self.variant_loci_genes[reagent_id] += [gene_id]

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

        First, we build a nice human-readable label of all of the components of the environment.  This is
        added to our global id-label hash.  Eventually, we will add each component of the environment to
        an environmental object, but needs to be modeled first

        This method also constructs the extrinsic genotype and sub-parts.
        Since the only way to get to the regents that comprise the extrinsic genotype, we build up an
        independent hash that contains all the extrinsic genotype components, and also add them to the graph.
        These are then accessible outside of this function from the global hash, and can be referenced
        by the genotype-to-phenotype method.

        Triples created:
        <environment_id> is an Individual
        <environment_id> rdfs:label <environment_label>
        <environment_id> has type GENO:environment


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

        logger.info("Processing environments")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        extrgeno_hash_by_env_id = {}
        env_hash = {}
        enviro_label_hash = {}
        geno = Genotype(g)
        envo = Environment(g)
        # pp = pprint.PrettyPrinter(indent=4)

        raw = '/'.join((self.rawdir, self.files['enviro']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (environment_num, condition_group, condition, values, units, comment, empty) = row

                if re.match("\\\\", values):
                    values = ''
                environment_id = 'ZFIN:' + environment_num.strip()
                if self.testMode and environment_num not in self.test_ids['environment']:
                    continue

                # We can start to build the extrinsic genotype using this file.
                # A single environment can have >1 row in the file, so we build a hash to store all the components
                # of the environment first.  Then we can build a label containing all the parts.
                # Using a strategy similar to what is used for genotypes to get the VSLCs and GVCs.

                # Clean up the units
                units = units.strip()
                if re.match('n/a', units.strip(), re.I) or units.strip() == '':
                    units = None
                if units is not None:
                    units = re.sub('[\(\)]','', units)

                # Clean up the values
                values = values.strip()
                if values == '' or re.match('n/a', values.strip(), re.I) or values == '?':
                    values = None
                else:
                    values = values.replace('<', 'less than')
                    values = values.replace('>', 'greater than')

                if values is None and units is not None:
                    # sometimes they just annotate to "ug" to indicate that scale, but no value
                    logger.warn("A unit without a value is indicated for env %s (condition %s)",
                                environment_id, condition)
                    values = ''

                if values is not None and units is None:
                    logger.warn("A value without a unit is indicated for env %s (condition %s)",
                                environment_id, condition)
                    units = ''

                # IDify the values and units (get rid of spaces and awkward chars)
                value_id = unit_id = None
                if values is not None:
                    value_id = values.replace(' ', '-')
                if units is not None:
                    if re.match('.*\/.*', units):
                        unit_id = re.sub(r"/", '-per-', units)
                    elif re.search('%', units):
                        unit_id = re.sub(r"%", 'percent', units)
                    else:
                        unit_id = units.replace(' ', '-')

                if comment == 'NULL' or comment == '':
                    comment = None

                if environment_id not in env_hash:
                    env_hash[environment_id] = []
                if environment_id not in extrgeno_hash_by_env_id:
                    extrgeno_hash_by_env_id[environment_id] = []

                if re.match('ZDB-(MRPHLNO|CRISPR|TALEN)', condition):
                    # this is the actual use of the knockdown reagent in the environment, including
                    # concentration.
                    # we build a node that contains this "application" of the reagent
                    morph_id = 'ZFIN:' + condition.strip()
                    if morph_id in self.id_label_map:
                        morph_label = self.id_label_map[morph_id]
                    else:
                        morph_label = morph_id
                        logger.warn('morph label not found %s', morph_id)

                    # Create the targeted sequence id, which includes units
                    if units is not None and values is not None:
                        applied_morph_id = '-'.join((condition.strip(), value_id, unit_id))
                        if values == '':
                            conc_label = '(' + units + ')'
                        else:
                            conc_label = '(' + values + ' ' + units + ')'
                    else:
                        applied_morph_id = condition + '_ns'
                        conc_label = '(n.s.)'

                    applied_morph_id = '_' + applied_morph_id
                    if self.nobnodes:
                        applied_morph_id = ':' + applied_morph_id
                    applied_morph_label = ' '.join((morph_label, conc_label))

                    # link the morpholino to the genes that it affects
                    ag = self.variant_loci_genes.get(morph_id)
                    # logger.info("%s affected genes %s", morph_id, pp.pformat(ag))
                    list_of_targeted_genes = []
                    if ag is None:
                        logger.warn("No affected genes for %s", morph_id)
                    else:
                        # create variant gene(s) that have been targeted by the reagent
                        for gid in ag:
                            if gid not in self.id_label_map:
                                # should not happen, except maybe in testing
                                logger.error("%s not in id-label-hash", gid)
                                glabel = gid
                            else:
                                glabel = self.id_label_map[gid]

                            targeted_gene_id = '-'.join((gid, applied_morph_id))
                            # these are not zfin resolvable, so make BNodes
                            targeted_gene_id = re.sub('(ZFIN)?:', '', targeted_gene_id)
                            targeted_gene_id = '_' + targeted_gene_id
                            if self.nobnodes:
                                targeted_gene_id = ':' + targeted_gene_id
                            targeted_gene_label = glabel + '<' + applied_morph_label + '>'

                            geno.addReagentTargetedGene(morph_id, gid, targeted_gene_id, targeted_gene_label)
                            self.id_label_map[targeted_gene_id] = targeted_gene_label
                            list_of_targeted_genes += [targeted_gene_id]

                    env_hash[environment_id] += [applied_morph_id]

                    self.id_label_map[applied_morph_id] = applied_morph_label

                    # add to the environment as a part
                    if re.match('ZDB-MRPHLNO', condition):
                        reagent_prefix = 'morpholino'
                    elif re.match('ZDB-CRISPR', condition):
                        reagent_prefix = 'crispr'
                    elif re.match('ZDB-TALEN', condition):
                        reagent_prefix = 'talen'
                    else:
                        reagent_prefix = 'reagent'

                    applied_morph_label = reagent_prefix+'['+applied_morph_label+']'
                    envo.addEnvironmentalCondition(applied_morph_id, applied_morph_label)
                    envo.addComponentAttributes(applied_morph_id, morph_id)  # TODO add value/unit
                    if comment is not None:
                        gu.addComment(g, morph_id, comment)
                    self.id_label_map[applied_morph_id] = applied_morph_label

                    # also add these reagent-targeted-genes to the extrinsic geno hash to be accessed by env id
                    extrgeno_hash_by_env_id[environment_id] += list_of_targeted_genes

                elif not re.match('ZDB.*', condition):
                    # create environmental components, and add to the hash
                    # cleanup the "condition" to remove non-id-friendly chars
                    cond_id = condition.strip()
                    cond_id = re.sub('\W+', '-', cond_id)

                    if units is not None and values is not None:
                        env_component_id = '-'.join((condition_group.strip(), cond_id.strip(), value_id, unit_id))
                        if values == '':
                            conc_label = '(' + units + ')'
                        else:
                            conc_label = '(' + values + ' ' + units + ')'
                    else:
                        env_component_id = '-'.join((condition_group.strip(), cond_id.strip()))
                        conc_label = ''

                    # make them blank nodes
                    env_component_id = '_' + env_component_id
                    if self.nobnodes:
                        env_component_id = ':' + env_component_id

                    env_condition = condition.strip()
                    if conc_label != '':
                        env_condition += ' ' + conc_label

                    env_component_label = condition_group + '[' + env_condition + ']'

                    self.id_label_map[env_component_id] = env_component_label
                    env_hash[environment_id] += [env_component_id]

                    if environment_id not in enviro_label_hash:
                        enviro_label_hash[environment_id] = [env_component_id]
                    else:
                        enviro_label_hash[environment_id].append(env_component_id)

                    # add each component to the environment as a part
                    envo.addEnvironmentalCondition(env_component_id, env_component_label)
                    if comment is not None:
                        gu.addComment(g, env_component_id, comment)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

                    # End of loop

        csvfile.close()

        logger.info("Building complex environments from components")

        self.environment_hash = env_hash

        # iterate through the env hash to build the full environment label
        for env_id in env_hash:
            environment_labels = []
            env_hash[env_id].sort()
            env_component_list = env_hash[env_id]
            for env_comp_id in env_component_list:
                env_comp_label = self.id_label_map[env_comp_id]
                environment_labels += [env_comp_label]
                envo.addComponentToEnvironment(env_id, env_comp_id)
            environment_labels.sort()
            env_label = 'Environment that includes: '+'; '.join(environment_labels)
            envo.addEnvironment(env_id, env_label)

        # iterate over the env hash to build the extrinsic genotypes that are lumped
        # together by the environment id.  these can then be accessed outside of this
        # function via the EXP id
        # for clarity i do this in the separate loop below, rather than combining it with the above

        logger.info("Building extrinsic genotype using environmental id groupings")

        for env_id in extrgeno_hash_by_env_id:
            extrgeno_labels = []
            extrgeno_hash_by_env_id[env_id].sort()

            # build an identifier for the extrinsic genotype, based on it's components
            extr_id = '-'.join(extrgeno_hash_by_env_id[env_id])

            # replace any colons (for example, when we materialize the b nodes for the components)
            extr_id = re.sub(':', '', extr_id)
            if self.nobnodes:
                extr_id = ':' + extr_id

            env_component_list = extrgeno_hash_by_env_id[env_id]
            for env_comp_id in env_component_list:
                env_comp_label = self.id_label_map[env_comp_id]
                extrgeno_labels += [env_comp_label]
                if len(env_component_list) > 1:
                    geno.addGeneTargetingReagentToGenotype(env_comp_id, extr_id)

            extrgeno_labels.sort()
            if len(extrgeno_labels) > 0:
                extr_label = '; '.join(extrgeno_labels)
                gu.addIndividualToGraph(g, extr_id, extr_label, geno.genoparts['extrinsic_genotype'])
                self.extrinsic_id_to_enviro_id_hash[env_id] = extr_id
                self.id_label_map[extr_id] = extr_label

                # logger.info('extr %s | %s', extr_id, extr_label)

        logger.info("Done with environments")

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
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['mappings']['file']))

        taxon_num = '7955'
        taxon_id = 'NCBITaxon:' + taxon_num
        taxon_label = 'Danio rerio'
        # genome_id = geno.makeGenomeID(taxon_id)
        geno.addGenome(taxon_id, taxon_label)

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (zfin_num, symbol, so_id, panel_symbol, chromosome, location, metric, empty) = row

                if self.testMode and zfin_num not in self.test_ids['gene'] + self.test_ids['allele']:
                    continue

                zfin_id = 'ZFIN:' + zfin_num.strip()
                if re.match('ZDB-GENE.*', zfin_num):
                    gu.addClassToGraph(g, zfin_id, None)  # assume type and label get added elsewhere
                    geno.addTaxon(taxon_id, zfin_id)
                elif re.match('ZDB-ALT.*', zfin_num):
                    gu.addIndividualToGraph(g, zfin_id, None)  # assume type and label get added elsewhere
                    geno.addTaxon(taxon_id, zfin_id)
                else:
                    continue
                    # skip any of the others

                # make the chromosome class
                chr_id = makeChromID(chromosome, taxon_id, 'CHR')
                # chr_label = makeChromLabel(chromosome, taxon_label)
                geno.addChromosomeClass(chromosome, taxon_id, taxon_label)

                pinfo = self._get_mapping_panel_info(panel_symbol)
                panel_label = ' '.join((panel_symbol,pinfo['type'],'map'))
                if pinfo is not None:
                    # add the panel as a genome build
                    panel_id = 'ZFIN:' + pinfo['id']
                    geno.addReferenceGenome(panel_id, panel_label, taxon_id)
                    gu.addSynonym(g, panel_id, panel_symbol)
                    gu.addDescription(g, panel_id, pinfo['name'])

                    # add the mapping-panel chromosome
                    chr_inst_id = makeChromID(chromosome, panel_id, 'MONARCH')
                    geno.addChromosomeInstance(chromosome, panel_id, panel_label, chr_id)
                    # add the feature to the mapping-panel chromosome
                    f = Feature(zfin_id, None, None)
                    f.addSubsequenceOfFeature(g, chr_inst_id)
                    f.loadAllProperties(g)
                    # TODO add the coordinates see https://github.com/JervenBolleman/FALDO/issues/24
                else:
                    logger.error("There's a panel (%s) we don't have info for", panel_symbol)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with chromosome mappings")
        return

    def _process_uniprot_ids(self, limit=None):
        """
        This method processes the mappings from ZFIN gene IDs to UniProtKB IDs.

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
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, self.files['uniprot']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, uniprot_id, empty) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue

                gene_id = 'ZFIN:' + gene_id.strip()
                uniprot_id = 'UniProtKB:' + uniprot_id.strip()

                geno.addGene(gene_id, gene_symbol)
                # TODO: Abstract to one of the model utilities
                gu.addIndividualToGraph(g, uniprot_id, None, geno.genoparts['polypeptide'])
                gu.addTriple(g, gene_id, gu.properties['has_gene_product'], uniprot_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with UniProt IDs")
        return

    def _process_human_orthos(self, limit=None):
        """
        This table provides ortholog mappings between zebrafish and humans.
        ZFIN has their own process of creating orthology mappings, that we take in addition to other
        orthology-calling sources (like PANTHER).
        We ignore the omim ids, and only use the gene_id.

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

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        logger.info("Processing human orthos")
        line_counter = 0
        geno = Genotype(g)

        raw = '/'.join((self.rawdir, self.files['human_orthos']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (zfin_id, zfin_symbol, zfin_name, human_symbol, human_name, omim_id, gene_id, empty) = row

                if self.testMode and zfin_id not in self.test_ids['gene']:
                    continue

                # Add the zebrafish gene.
                zfin_id = 'ZFIN:' + zfin_id.strip()
                geno.addGene(zfin_id, zfin_symbol, None, zfin_name)

                # Add the human gene.
                gene_id = 'NCBIGene:' + gene_id.strip()
                geno.addGene(gene_id, human_symbol, None, human_name)

                # make the association
                assoc = OrthologyAssoc(self.name, zfin_id, gene_id)
                # we don't know anything about the orthology type, so we just use the default
                assoc.load_all_properties(g)
                assoc.add_association_to_graph(g)

                # FIXME we have requested that ZFIN add evidence codes and papers to the orthology calls
                # TODO we can get this from zfin mine
                # (This data is in their web front-end, but not in the downloads)

                if not self.testMode and limit is not None and line_counter > limit:
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

        # zfin uses free-text modifiers, but we need to convert them to proper PATO classes for the mapping
        mod_id = modifier

        modifiers = {
            'abnormal': 'PATO:0000460',
            'normal': 'PATO:0000461'
        }
        if modifier in modifiers.keys():
            mod_id = modifiers.get(modifier)

        key = self._make_zpkey(superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, mod_id)
        mapping = self.zp_map.get(key)

        if mapping is None:
            # pass
            if modifier == 'normal':
                pass
                # logger.info("Normal phenotypes not yet supported")
            else:
                logger.warn("Couldn't map ZP id to %s with modifier %s", "_".join(
                    (superterm1_id, subterm1_id, quality_id, superterm2_id, subterm2_id, mod_id)),
                            modifier)
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


    def _get_mapping_panel_info(self, panel):
        panel_hash = {
            'HS': {'id': 'ZDB-REFCROSS-000320-1', 'name': 'Heat Shock', 'type': 'meiotic', 'num_meioses': 42},
            'GAT': {'id': 'ZDB-REFCROSS-990308-7', 'name': 'Gates et al', 'type': 'meiotic', 'num_meioses' : 96},
            'LN54': {'id': 'ZDB-REFCROSS-990426-6', 'name': 'Loeb/NIH/5000/4000', 'dose': '4000 rads', 'type': 'Radiation Hybrid'},
            'MGH': {'id': 'ZDB-REFCROSS-980521-11', 'name': 'Boston MGH Cross', 'type': 'meiotic', 'num_meioses' : 40},
            'MOP': {'id': 'ZDB-REFCROSS-980526-5', 'name': 'Mother of Pearl', 'type': 'meiotic', 'num_meioses' : 96},
            'T51': {'id': 'ZDB-REFCROSS-990707-1', 'name': 'Goodfellow T51', 'dose': '3000 rads', 'type': 'Radiation Hybrid'},
        }
        p = None
        if panel in panel_hash:
            p = panel_hash[panel]

        return p

    def _make_variant_locus_id(self, gene_id, allele_id):
        """
        A convenience method to uniformly create variant loci.
        If we want to materialize these in the monarch space, then we wrap with the self.make_id function.
        :param gene_id:
        :param allele_id:
        :return:
        """

        i = '-'.join((gene_id, allele_id))
        i = '_' + re.sub('(ZFIN)?:', '', i)
        if self.nobnodes:
            i = ':' + i

        return i

    def _make_effective_genotype_id(self, intrinsic_id, extrinsic_id):
        effective_genotype_id = self.make_id('-'.join((intrinsic_id, extrinsic_id)))

        return effective_genotype_id

    def _remove_backslash_r(self, filename, encoding):
        """
        A helpful utility to remove '\r' from any file
        :param filename:
        :return:
        """

        f = open(filename, 'r', encoding=encoding, newline='\n')
        contents = f.read()
        f.close()
        contents = re.sub(r'\r', '', contents)
        with open(filename, "w") as f:
            f.truncate()
            f.write(contents)

        return

    def getTestSuite(self):
        import unittest
        from tests.test_zfin import ZFINTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(ZFINTestCase)

        return test_suite