import csv
import re
import logging
#from Bio.Seq import Seq

from dipper.utils import pysed
from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Genotype import Genotype
from dipper.models.OrthologyAssoc import OrthologyAssoc
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.GenomicFeature import Feature,makeChromID
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
            
        self.kd_reagent_hash = {'kd_reagent_id': {}, 'kd_reagent_label': {}, 'gene_label': {}}
        self.wildtype_hash = {'id': {}, 'symbol': {}}
        self.label_hash = {'gene_label': {}, 'allele_label': {}, 'construct_label': {}, 'genotype_label': {},
                           'background_label': {}, 'morpholino_label': {}}
        self.genotype_id_to_background_id_hash = {'genotype_id': {}}
        self.extrinsic_id_to_enviro_id_hash = {'extrinsic_id': {}}

        # These must be processed in a specific order
        self._process_wildtypes(limit)  # Must be processed before wildtype_expression
        self._process_genotype_backgrounds(limit)
        self._process_genotype_features(limit)
        self._process_morpholinos(limit)  # Process before talens/crisprs
        # NOTE: TALENs/CRISPRs are not currently added to the extrinsic genotype, but are added to the env label.
        self._process_talens(limit)  # Leaving TALENs out until further review.
        self._process_crisprs(limit)  # Leaving CRISPRs out until further review.
        self._process_pheno_enviro(limit)  #  Must be processed after morpholinos/talens/crisprs
        self._process_feature_affected_genes(limit)
        self._process_g2p(limit)

        #self._process_wildtype_expression(limit)
        self._process_gene_marker_relationships(limit)
        self._process_features(limit)
        self._process_genes(limit)
        self._process_genbank_ids(limit)
        self._process_uniprot_ids(limit)
        self._process_human_orthos(limit)
        #self._process_anatomy(limit)
        self._process_stages(limit)
        self._process_pubinfo(limit)
        self._process_pub2pubmed(limit)
        self._process_mappings(limit)

        logger.info("Finished parsing.")

        self.load_bindings()

        return

    def _process_genotype_features(self, limit=None):
        """
        We don't actually know the allele pairs, so we'll store some info in a hashmap for post-processing
        :param limit:
        :return:
        """
        raw = '/'.join((self.rawdir, self.files['geno']['file']))

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        gu = GraphUtils(curie_map.get())
        geno_hash = {}
        gvc_hash = {}
        vslc_label_hash = {}
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
                background_id = self.genotype_id_to_background_id_hash['genotype_id'].get(genotype_id)
                #if background_id is not None:
                    #genotype_name = genotype_name+' ['+self.label_hash['background_label'].get(background_id)+']'
                #else:
                    #genotype_name = genotype_name+' [n.s.]'
                #print(genotype_name)
                geno = Genotype(g)
                #gt = geno.addGenotype(genotype_id, genotype_name)
                gt = geno.addGenotype(genotype_id, None)
                if genotype_id not in geno_hash:
                    geno_hash[genotype_id] = {}
                genoparts = geno_hash[genotype_id]
                #if self.label_hash['genotype_label'].get(genotype_id) is None:
                    #self.label_hash['genotype_label'][genotype_id] = genotype_name

                # FIXME: This seems incorrect. Shouldn't the types available here be added to the sequence alteration?
                # reassign the allele_type to a proper GENO or SO class
                allele_type = self._map_allele_type_to_geno(allele_type)

                allele_id = 'ZFIN:' + allele_id.strip()
                geno.addAllele(allele_id, allele_name, allele_type)

                if gene_id is not None and gene_id.strip() != '':
                    gene_id = 'ZFIN:' + gene_id.strip()
                    geno.addGene(gene_id, gene_symbol)

                    if gene_id is not None and gene_id not in self.label_hash['gene_label']:
                        self.label_hash['gene_label'][gene_id] = gene_symbol
                    if allele_id is not None and allele_id not in self.label_hash['allele_label']:
                        self.label_hash['allele_label'][allele_id] = allele_name

                    # if it's a transgenic construct, then we'll have to get the other bits
                    if construct_id is not None and construct_id.strip() != '':
                        construct_id = 'ZFIN:' + construct_id.strip()
                        geno.addDerivesFrom(allele_id, construct_id)
                        if construct_id not in self.label_hash['construct_label']:
                            self.label_hash['construct_label'][construct_id] = construct_name

                    # allele to gene
                    geno.addAlleleOfGene(allele_id, gene_id)
                    if gene_id not in genoparts:
                        genoparts[gene_id] = [allele_id]
                    else:
                        genoparts[gene_id].append(allele_id)

                    if zygosity == 'homozygous':
                        genoparts[gene_id].append(allele_id)  # add the allele again
                    elif zygosity == 'unknown':
                        genoparts[gene_id].append('?')  # we'll just use this as a convention for unknown
                    #elif (zygosity == 'complex'):  # what are these?
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
                    if allele_id is not None and allele_id not in self.label_hash['allele_label']:
                        self.label_hash['allele_label'][allele_id] = allele_name

                if limit is not None and line_counter > limit:
                    break

                # end loop through file
            # now loop through the geno_hash, and build the vslcs

            for gt in geno_hash:
                if gt not in gvc_hash:
                    gvc_hash[gt] = {}
                gvcparts = gvc_hash[gt]
                vslc_counter = 0

                for gene_id in geno_hash.get(gt):
                    gene_label = self.label_hash['gene_label'].get(gene_id)
                    variant_locus_parts = geno_hash.get(gt).get(gene_id)
                    if gene_id in variant_locus_parts:
                        # reset the gene_id to none
                        gene_id = None

                    allele1_id = variant_locus_parts[0]
                    allele1_label = self.label_hash['allele_label'].get(allele1_id)
                    allele2_id = None
                    allele2_label = None
                    zygosity_id = None
                    # making the assumption that there are not more than 2 variant_locus_parts
                    if len(variant_locus_parts) > 1:
                        allele2_id = variant_locus_parts[1]
                        allele2_label = self.label_hash['allele_label'].get(allele2_id)
                    if allele2_id is not None:
                        if allele2_id == '?':
                            zygosity_id = geno.zygosity['indeterminate']
                            allele2_id = None
                            allele2_label = '?'
                        elif allele2_id == 'complex':
                            pass  # not sure what to assign here
                        elif allele1_id != allele2_id:
                            zygosity_id = geno.zygosity['heterozygous']
                            #print('heterozygous pair='+allele1_id+'_'+allele2_id)
                        elif allele1_id == allele2_id:
                            zygosity_id = geno.zygosity['homozygous']
                    else:
                        zygosity_id = geno.zygosity['indeterminate']

                    # create the vslc
                    if gene_id is None:
                        gn = ''
                    else:
                        gn = gene_id
                    if (allele2_id is None):
                        a2 = ''
                    else:
                        a2 = allele2_id

                    vslc_id = self.make_id(('-').join((gn, allele1_id, a2)))
                    #print(g+'_'+allele1_id+'_'+a2)
                    #print(gene_label)
                    #print(allele1_label)
                    #print(allele2_label)
                    #print(construct_label)
                    if gene_label is not None and allele1_label is not None and allele2_label is not None:
                        vslc_label = gene_label+'<'+allele1_label+'>/'+gene_label+'<'+allele2_label+'>'
                    elif gene_label is None and allele1_label is not None and allele2_label is not None:
                        vslc_label = '<'+allele1_label+'>/<'+allele2_label+'>'
                    elif gene_label is not None and allele1_label is not None and allele2_label is None:
                        vslc_label = gene_label+'<'+allele1_label+'>'
                    elif gene_label is None and allele1_label is not None and allele2_label is None:
                        vslc_label = '<'+allele1_label+'>'
                    else:
                        logger.error('No VSLC label created.')
                        vslc_label = ''
                    #print(vslc_label)

                    gu.addIndividualToGraph(g, vslc_id, vslc_label, geno.genoparts['variant_single_locus_complement'])
                    geno.addPartsToVSLC(vslc_id, allele1_id, allele2_id, zygosity_id)
                    # Remove this since I am now adding the VSLC to the GVC?
                    #geno.addVSLCtoParent(vslc_id,gt)

                    gt_vslc = gt+'vslc'

                    # FIXME: Refactor this. Don't like the counter approach.
                    #gvc_hash[vslc_id] = {};
                    if vslc_counter == 0:
                        gvcparts[gt_vslc] = [vslc_id]
                    elif vslc_id not in gvcparts:
                        gvcparts[gt_vslc].append(vslc_id)
                    vslc_counter += 1

                    if gt_vslc not in vslc_label_hash:
                        vslc_label_hash[gt_vslc] = [vslc_label]
                    elif vslc_id not in vslc_label_hash[gt_vslc]:
                        vslc_label_hash[gt_vslc].append(vslc_label)

                #print(gvcparts)

                    # end loop through geno_hash
            # now loop through the gvc_hash, and build the gvc
            # TODO: Possible to pass through VSLC label, assemble GVC label?
            for gt in gvc_hash:
                gvc_id = '<empty>'
                gvc_ids = []
                gvc_labels = []
                for vslc_id in gvc_hash.get(gt):
                    genomic_variation_complement_parts = gvc_hash.get(gt).get(vslc_id)
                    gvc_label_parts = vslc_label_hash[vslc_id]
                    #print(gvc_label_parts)
                    #print(genomic_variation_complement_parts)
                    # FIXME: Change to make_id after QA.
                    gvc_id = self.make_id('-'.join(genomic_variation_complement_parts))
                    #gvc_id = ('_split_').join(genomic_variation_complement_parts)
                    gvc_label = '; '.join(gvc_label_parts)
                    # Add the GVC
                    gu.addIndividualToGraph(g, gvc_id, gvc_label, geno.genoparts['genomic_variation_complement'])
                    #print(gvc_id)

                    # Add the VSLCs to the GVC
                    for i in genomic_variation_complement_parts:
                        #geno.addVSLCtoParent(i,gt)
                        gu.addTriple(g, gvc_id, geno.object_properties['has_alternate_part'], i)

                background_id = self.genotype_id_to_background_id_hash['genotype_id'].get(gt)
                if background_id is not None:
                    genotype_name = gvc_label+' ['+self.label_hash['background_label'].get(background_id)+']'
                else:
                    genotype_name = gvc_label+' [n.s.]'
                #print(genotype_name)
                geno = Genotype(g)
                geno.addGenotype(gt, genotype_name)
                #gt = geno.addGenotype(genotype_id, None)
                if gt not in geno_hash:
                    geno_hash[gt] = {}
                #genoparts = geno_hash[genotype_id]
                if self.label_hash['genotype_label'].get(gt) is None:
                    self.label_hash['genotype_label'][gt] = genotype_name
                    #print(genotype_name)

                # Add the GVC to the genotype
                #print(gt)
                gu.addTriple(g, gt, geno.object_properties['has_alternate_part'], gvc_id)

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

        logger.info("Processing genotype backgrounds")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir,self.files['backgrounds']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (genotype_id, genotype_name, background_id, empty) = row

                genotype_id = 'ZFIN:' + genotype_id.strip()
                background_id = 'ZFIN:' + background_id.strip()
                if self.genotype_id_to_background_id_hash['genotype_id'].get(genotype_id) is None:
                    self.genotype_id_to_background_id_hash['genotype_id'][genotype_id] = background_id

                # Add the taxon as a class
                taxon_id = 'NCBITaxon:7955'  # Danio rerio
                gu.addClassToGraph(self.graph, taxon_id, None)
                geno.addTaxon(taxon_id, background_id)

                # Add genotype
                # TODO: Need to adjust the genotype name to properly formatted labels
                # Need to break apart the single gene/dual allele notation
                # genotype_name = re.sub('</sup>','>', (re.sub('<sup>','<',genotype_name)))
                #print(genotype_name)
                geno.addGenotype(genotype_id, genotype_name)

                # Add background to the genotype
                geno.addGenomicBackgroundToGenotype(background_id, genotype_id)

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with genotype backgrounds")
        return


    def _process_wildtypes(self, limit=None):
        """
        This table provides the genotype IDs, name, and abbreviation of the wildtype genotypes.

        Triples created:
        <genotype id> a GENO:wildtype
        <genotype id> rdfs:label genotype_abbreviation
        <genotype id> dc:description genotype_name

        :param limit:
        :return:
        """

        logger.info("Processing wildtype genotypes")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['wild']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (genotype_id, genotype_name, genotype_abbreviation, empty) = row

                genotype_id = 'ZFIN:' + genotype_id.strip()

                # Add genotype to graph with label, type=wildtype, and description.
                geno.addGenotype(genotype_id, genotype_abbreviation, geno.genoparts['wildtype'], genotype_name)

                if self.label_hash['background_label'].get(genotype_id) is None:
                    self.label_hash['background_label'][genotype_id] = genotype_name

                if self.label_hash['genotype_label'].get(genotype_id) is None:
                    self.label_hash['genotype_label'][genotype_id] = '['+genotype_name+']'

                # Build the hash for the wild type genotypes.
                if self.wildtype_hash['id'].get(genotype_name) is None:
                    self.wildtype_hash['id'][genotype_name] = genotype_id
                    self.wildtype_hash['symbol'][genotype_name] = genotype_abbreviation

                #if self.kd_reagent_hash['kd_reagent_label'].get(morpholino_id) is None:
                    #self.kd_reagent_hash['kd_reagent_label'][morpholino_id] = morpholino_symbol
                #if self.kd_reagent_hash['gene_label'].get(gene_id) is None:
                    #self.kd_reagent_hash['gene_label'][gene_id] = gene_symbol

                if limit is not None and line_counter > limit:
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

                #genotype_id = 'ZFIN:' + genotype_id.strip()

                genotype_id = self.wildtype_hash['id'][genotype_name]
                genotype_id = 'ZFIN:' + genotype_id.strip()

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

        logger.info("Processing stages")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['stage']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (stage_id, stage_obo_id, stage_name, begin_hours, end_hours, empty) = row

                stage_id = 'ZFIN:' + stage_id.strip()

                # Make ID for the beginning hour, add to graph.
                #begin_hour_id = self.make_id(begin_hours)
                #gu.addIndividualToGraph(self.graph,begin_hour_id,begin_hours+' hours',gu.datatype_properties['hours'])

                # Make ID for the ending hour, add to graph.
                #end_hour_id = self.make_id(end_hours)
                #gu.addIndividualToGraph(self.graph,end_hour_id,end_hours+' hours',gu.datatype_properties['hours'])

                # Add the stage as an individual.
                gu.addIndividualToGraph(self.graph, stage_id, stage_name, stage_obo_id)

                # Add the beginning and end hours of the development stage.
                #gu.addTriple(self.graph,stage_id,gu.object_properties['existence_starts_at'],begin_hour_id)
                #gu.addTriple(self.graph,stage_id,gu.object_properties['existence_ends_at'],end_hour_id)

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with stages")
        return

    def _process_anatomy(self, limit=None):
        """
        This table provides mappings between ZFIN stage IDs and ZFA anatomy terms,
        indicating the starting and ending development stage for the anatomy term.

        Triples created:
        <anatomy_id> an individual
        <anatomy_id> rdfs:label anatomy_name

        <anatomy_id> uberon:existence_starts_at begin_hour_id
        <anatomy_id> uberon:existence_ends_at end_hour_id

        :param limit:
        :return:
        """
        # FIXME I don't think we need this [nlw]
        logger.info("Processing anatomy")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['anatomy']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (anatomy_id, anatomy_name, start_stage_id, end_stage_id, empty) = row

                start_stage_id = 'ZFIN:' + start_stage_id.strip()
                end_stage_id = 'ZFIN:' + end_stage_id.strip()

                #Is this correct? Should an anatomy term be declared as it's own type, or just an individual?
                gu.addIndividualToGraph(self.graph, anatomy_id, anatomy_name)

                gu.addTriple(self.graph, anatomy_id, gu.object_properties['existence_starts_at'], start_stage_id)
                gu.addTriple(self.graph, anatomy_id, gu.object_properties['existence_ends_at'], end_stage_id)

                if (limit is not None and line_counter > limit):
                    break

        logger.info("Done with anatomy")
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
                    geno.addGenotype(effective_genotype_id, effective_genotype_label, geno.genoparts['effective_genotype'])
                    if extrinsic_genotype_label is not None and extrinsic_genotype_label != '':
                        geno.addParts(extrinsic_geno_id, effective_genotype_id)
                    if intrinsic_genotype_label is not None and intrinsic_genotype_label != '':
                        geno.addParts(genotype_id, effective_genotype_id)

                phenotype_id = self._map_sextuple_to_phenotype(superterm1_id, subterm1_id, quality_id,
                                                               superterm2_id, subterm2_id, modifier)

                if (phenotype_id is None):
                    #add this phenotype to a set to report at the end
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
                    assoc.addAssociationNodeToGraph(g)
                    # FIXME: Hanging the environment and stages on the G2P association for now.
                    # Add environment to G2P association.
                    gu.addTriple(g, assoc, gu.object_properties['has_environment_qualifier'], env_id)
                    # Add starting stage to G2P association.
                    gu.addTriple(g, assoc, gu.object_properties['has_begin_stage_qualifier'], start_stage_id)
                    # Add ending stage to G2P association.
                    gu.addTriple(g, assoc, gu.object_properties['has_end_stage_qualifier'], end_stage_id)
                else:
                    # add normal phenotypes
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
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, ncbi_gene_id, empty) = row
                
                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue
                    
                geno = Genotype(g)

                gene_id = 'ZFIN:'+gene_id.strip()
                ncbi_gene_id = 'NCBIGene:'+ncbi_gene_id.strip()

                geno.addGene(gene_id, gene_symbol)
                gu.addEquivalentClass(g, gene_id, ncbi_gene_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with genes")
        return

    def _process_features(self, limit=None):
        """
        This module provides information for the intrinsic and extrinsic genotype features of zebrafish.

         sequence alteration ID, SO type, abbreviation, and relationship to
        the affected gene, with the gene's ID, symbol, and SO type (gene/pseudogene).

        Triples created:
        <gene id> a class:
        :param limit:
        :return:
        """

        logger.info("Processing features")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['features']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (genomic_feature_id, feature_so_id, genomic_feature_abbreviation, genomic_feature_name,
                genomic_feature_type, mutagen, mutagee, construct_id, construct_name, construct_so_id,empty) = row

                genomic_feature_id = 'ZFIN:' + genomic_feature_id.strip()

                gu.addIndividualToGraph(self.graph, genomic_feature_id, genomic_feature_name, feature_so_id)

                if construct_id is not None and construct_id != '':
                    construct_id = 'ZFIN:' + construct_id.strip()
                    geno.addConstruct(construct_id, construct_name, construct_so_id)
                    # FIXME: Need the appropriate relationship between the construct and the mutation/alteration.
                    # Derives from? Parent = construct, child = allele/feature?
                    geno.addDerivesFrom(genomic_feature_id, construct_id)

                # TODO: Have available a mutagen and mutagee (adult males, embryos, etc.)
                # How should this be modeled?
                # Mutagens: CRISPR, EMS, ENU, DNA, g-rays, not specified, spontaneous, TALEN, TMP, zinc finger nuclease
                # TODO: make a mapping function for mutagens, if needed.
                # Mutagees: adult females, adult males, embryos, not specified, sperm

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with features")
        return

    def _process_feature_affected_genes(self, limit=None):
        """
        This table provides the sequence alteration ID, SO type, abbreviation, and relationship to
        the affected gene, with the gene's ID, symbol, and SO type (gene/pseudogene).

        Triples created:
        <gene id> a class:
        <gene id> rdfs:label gene_symbol
        <gene id> subclass of type gene/pseudogene

        <variant locus id> is a GENO:allele
        <variant locus id> rdfs:label <variant_locus_label>
        <variant locus id> is an allele of <gene id>
        <variant locus id> has alternate part <sequence alteration id>

        <sequence alteration id> is an allele of <gene id>
        <sequence alteration id> rdf:type <sequence alteration type>

        :param limit:
        :return:
        """

        logger.info("Processing feature affected genes")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['feature_affected_gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                geno = Genotype(self.graph)
                (genomic_feature_id, feature_so_id, genomic_feature_abbreviation, gene_symbol, gene_id, gene_so_id,
                 genomic_feature_marker_relationship, empty) = row

                # Sequence alteration types present in file: SO:0000159 - deletion, SO:0000199 - translocation,
                # SO:0000667 - insertion, SO:0001059 - sequence_alteration,
                # SO:0001060 - sequence_variant, SO:0001218 - transgenic insertion, SO:1000005 - complex_substitution,
                # SO:1000008 - point_mutation, SO:1000029 - chromosomal_deletion, SO:1000032 - indel

                genomic_feature_id = 'ZFIN:' + genomic_feature_id.strip()
                # TODO: Can build the variant locus, sequence alteration, and sequence alteration type here.

                gene_id = 'ZFIN:' + gene_id.strip()
                # NOTE: There are a few pseudogenes in the file, using the gene_so_id will identify them.
                geno.addGene(gene_id, gene_symbol, gene_so_id)

                # Add variant locus (allele)
                # TODO: Confirm this is correct/holds true for all variant loci.
                # Should a different format be used for transgenic insertions?
                variant_locus_label = gene_symbol+'<'+genomic_feature_abbreviation+'>'
                #print(variant_locus)
                variant_locus_id = self.make_id(genomic_feature_id+gene_id)

                # FIXME:Do we want to filter out the translocations/deletions, or include them?
                geno.addAllele(variant_locus_id,variant_locus_label)

                #gu.addIndividualToGraph(self.graph,genomic_feature_id,genomic_feature_abbreviation,feature_so_id)

                # NOTE: Most feature_marker_relationship entries are 'is allele of' but there are some entries with
                # 'markers missing' corresponds to SO:1000029 - chromosomal_deletion) or
                # 'markers moved' (corresponds to SO:1000199 - translocation).
                # For now, only indicating as an allele of gene if the relationship is 'is allele of.'
                # TODO: Confirm that this conditional is the correct approach.
                # Doesn't a translocation or deletion count as an allele?
                if (genomic_feature_marker_relationship == 'is allele of'):
                    # FIXME: Should this be the variant locus ID or the genomic_feature_id? I think the VL ID.
                    # Add the gene to the allele.
                    geno.addAlleleOfGene(variant_locus_id, gene_id)
                # TODO: For the other relationships, is there a 'translocation_of' or 'deletion_of' that can be used?

                sequence_alteration_type = feature_so_id

                # Add the sequence alteration id, label, and type
                geno.addSequenceAlteration(genomic_feature_id, genomic_feature_abbreviation, sequence_alteration_type)

                # Add the sequence alteration to the variant locus
                geno.addSequenceAlterationToVariantLocus(genomic_feature_id, variant_locus_id)

                if limit is not None and line_counter > limit:
                    break

        logger.info("Done with feature affected genes")
        return

    def _process_gene_marker_relationships(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing gene marker relationships")
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = '/'.join((self.rawdir, self.files['gene_marker_rel']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol,
                 marker_id, marker_so_id, marker_symbol, relationship, empty) = row

                geno = Genotype(self.graph)

                gene_id = 'ZFIN:'+gene_id.strip()
                # File contains genes (SO:0000704) and psuedogenes (SO:0000336).
                if gene_so_id == 'SO:0000704':
                    geno.addGene(gene_id, gene_symbol)
                elif gene_so_id == 'SO:0000336':
                    gu.addIndividualToGraph(self.graph, gene_id, gene_symbol, gene_so_id)

                marker_id = 'ZFIN:'+marker_id.strip()
                gu.addIndividualToGraph(self.graph, marker_id, marker_symbol, marker_so_id)

                # TODO: Map gene-marker relationships.
                # Gene-marker relationships: clone contains gene, coding sequence of, contains polymorphism,
                # gene contains small segment, gene encodes small segment, gene has artifact,
                # gene hybridized by small segment, gene produces transcript, gene product recognized by antibody,
                # knockdown reagent targets gene, promoter of, transcript targets gene

                if limit is not None and line_counter > limit:
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
        process_pubinfo method, but this mapping includes additional internal pub to pubmed mappings.

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

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

        # TODO: We have the sequence information for each of the targeting reagents. How to model?
    def _process_morpholinos(self, limit=None):
        """
        This method processes the morpholino knockdown reagents, creating triples for the
        morpholinos and passing the morpholino data into a hash map for use in the pheno_enviro method.

        Triples created:
        <morpholino_id> is a gene_targeting_reagent
        <morpholino_id> rdfs:label <morpholino_symbol>
        <morpholino_id> has type <morpholino_so_id>
        <morpholino_id> has comment <note>

        <publication_id> is an individual
        <publication_id> mentions <morpholino_id>
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
        raw = '/'.join((self.rawdir, self.files['morph']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, morpholino_id, morpholino_so_id,
                 morpholino_symbol, morpholino_sequence, publication, note) = row

                if self.testMode and morpholino_id not in self.test_ids['morpholino']:
                    continue

                geno = Genotype(g)
                morpholino_id = 'ZFIN:'+morpholino_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()

                # TODO: map target sequence to morpholino.
                # FIXME: Is this correct? Leaving out for now.
                # Commenting out for now
                # Is the reverse complement of the morpholino sequence the target sequence or, like miRNAs, is there
                # a seed sequence that is the target sequence and it is not the full reverse complement of the sequence?
                # Also, does the morpholino require the exact sequence match or can there be mismatches?
                # Take the morpholino sequence and get the reverse complement as the target sequence.
                #seq = Seq(morpholino_sequence)
                #target_sequence = seq.reverse_complement()
                #print(seq)
                #print(target_sequence)
                #print(morpholino_id)

                # FIXME: This is incorrect, as it requires the concentration if available, and is related to the extrinsic genotype.
                # Should add the morpholino as a typed individual instead. Same for TALENs/CRISPRs.
                geno.addGeneTargetingReagent(morpholino_id,morpholino_symbol,morpholino_so_id)
                # Now adding the reagent targeted gene in the pheno_environment processing function.
                #geno.addReagentTargetedGene(morpholino_id,gene_id, gene_id)

                # Add publication
                if(publication != ''):
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g, pub_id, None)
                    gu.addTriple(g, pub_id, gu.properties['mentions'], morpholino_id)

                # Add comment?
                if(note != ''):
                    gu.addComment(g, morpholino_id, note)

                # Build the hash for the reagents and the gene targets
                if self.kd_reagent_hash['kd_reagent_id'].get(morpholino_id) is None:
                    reagent_target = []
                    reagent_target.append(gene_id)
                    self.kd_reagent_hash['kd_reagent_id'][morpholino_id] = reagent_target
                else:
                    reagent_target = self.kd_reagent_hash['kd_reagent_id'][morpholino_id]
                    reagent_target.append(gene_id)
                    self.kd_reagent_hash['kd_reagent_id'][morpholino_id] = reagent_target

                if self.kd_reagent_hash['kd_reagent_label'].get(morpholino_id) is None:
                    self.kd_reagent_hash['kd_reagent_label'][morpholino_id] = morpholino_symbol
                if self.kd_reagent_hash['gene_label'].get(gene_id) is None:
                    self.kd_reagent_hash['gene_label'][gene_id] = gene_symbol

                if self.label_hash['morpholino_label'].get(morpholino_id) is None:
                    self.label_hash['morpholino_label'][morpholino_id] = morpholino_symbol

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with Morpholinos")
        return


    def _process_talens(self, limit=None):
        """
        TALENs are artificial restriction enzymes that can be used for genome editing in situ.

        Triples created:
        <talen_id> is a gene_targeting_reagent
        <talen_id> rdfs:label <talen_symbol>
        <talen_id> has type <talen_so_id>
        <talen_id> has comment <note>

        <publication_id> is an individual
        <publication_id> mentions <talen_id>
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
                 talen_symbol, talen_target_sequence_1, talen_target_sequence_2, publication, note) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue

                geno = Genotype(g)
                talen_id = 'ZFIN:'+talen_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()

                geno.addGeneTargetingReagent(talen_id, talen_symbol, talen_so_id)
                # Now adding the reagent targeted gene in the pheno_environment processing function.
                #geno.addReagentTargetedGene(talen_id,gene_id,gene_id)

                # Add publication
                if(publication != ''):
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g,pub_id,None)
                    gu.addTriple(g, pub_id, gu.properties['mentions'], talen_id)

                # Add comment?
                if(note != ''):
                    gu.addComment(g, talen_id, note)

                # Build the hash for the reagents and the gene targets
                if self.kd_reagent_hash['kd_reagent_id'].get(talen_id) is None:
                    reagent_target = []
                    reagent_target.append(gene_id)
                    self.kd_reagent_hash['kd_reagent_id'][talen_id] = reagent_target
                else:
                    reagent_target = self.kd_reagent_hash['kd_reagent_id'][talen_id]
                    reagent_target.append(gene_id)
                    self.kd_reagent_hash['kd_reagent_id'][talen_id] = reagent_target

                if self.kd_reagent_hash['kd_reagent_label'].get(talen_id) is None:
                    self.kd_reagent_hash['kd_reagent_label'][talen_id] = talen_symbol
                if self.kd_reagent_hash['gene_label'].get(gene_id) is None:
                    self.kd_reagent_hash['gene_label'][gene_id] = gene_symbol

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        logger.info("Done with TALENS")
        return

    def _process_crisprs(self, limit=None):
        """
        CRISPRs are knockdown reagents, working similar to RNAi but at the transcriptional level instead of mRNA level.

        Triples created:
        <crispr_id> is a gene_targeting_reagent
        <crispr_id> rdfs:label <crispr_symbol>
        <crispr_id> has type <crispr_so_id>
        <crispr_id> has comment <note>

        <publication_id> is an individual
        <publication_id> mentions <crispr_id>
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
        raw = '/'.join((self.rawdir, self.files['crispr']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (gene_id, gene_so_id, gene_symbol, crispr_id, crispr_so_id,
                crispr_symbol, crispr_target_sequence, publication, note) = row

                if self.testMode and gene_id not in self.test_ids['gene']:
                    continue

                geno = Genotype(g)
                crispr_id = 'ZFIN:'+crispr_id.strip()
                gene_id = 'ZFIN:'+gene_id.strip()

                geno.addGeneTargetingReagent(crispr_id, crispr_symbol, crispr_so_id)
                # Now adding the reagent targeted gene in the pheno_environment processing function.
                #geno.addReagentTargetedGene(crispr_id, gene_id, gene_id)

                # Add publication
                if publication != '':
                    pub_id = 'ZFIN:'+publication.strip()
                    gu.addIndividualToGraph(g, pub_id, None)
                    gu.addTriple(g, pub_id, gu.properties['mentions'], crispr_id)

                # Add comment?
                if note != '':
                    gu.addComment(g, crispr_id, note)

                # Build the hash for the reagents and the gene targets
                if self.kd_reagent_hash['kd_reagent_id'].get(crispr_id) is None:
                    reagent_target = []
                    reagent_target.append(gene_id)
                    self.kd_reagent_hash['kd_reagent_id'][crispr_id] = reagent_target
                else:
                    reagent_target = self.kd_reagent_hash['kd_reagent_id'][crispr_id]
                    reagent_target.append(gene_id)
                    self.kd_reagent_hash['kd_reagent_id'][crispr_id] = reagent_target

                if self.kd_reagent_hash['kd_reagent_label'].get(crispr_id) is None:
                    self.kd_reagent_hash['kd_reagent_label'][crispr_id] = crispr_symbol
                if self.kd_reagent_hash['gene_label'].get(gene_id) is None:
                    self.kd_reagent_hash['gene_label'][gene_id] = gene_symbol

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        logger.info("Done with CRISPRs")
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

                # FIXME: Can the environment serve as the extrinsic genotype ID?
                # For now, making an ID from environment ID only, but may want to revisit.
                extrinsic_geno_id = self.make_id(environment_id)
                self.extrinsic_id_to_enviro_id_hash['extrinsic_id'][extrinsic_geno_id] = environment_id

                geno = Genotype(self.graph)

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
                        #FIXME: Better way to indicate > and < ?
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
                elif not re.match('ZDB.*',condition):
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
                gu.addIndividualToGraph(self.graph, environment_id, environment_label, gu.datatype_properties['environment'])

            # Now process through the extrinsic_part_hash to produce targeted_gene_subregion and targeted_gene_variant
            #print(extrinsic_part_hash)
            tgc_hash = {}
            for extrinsic_geno_id in extrinsic_part_hash:
                #print(extrinsic_part_hash[extrinsic_geno_id])

                geno = Genotype(self.graph)
                ex_geno = geno.addGenotype(extrinsic_geno_id, None, geno.genoparts['extrinsic_genotype'])
                for condition in extrinsic_part_hash[extrinsic_geno_id]:
                    kd_reagent_conc_id = kd_reagent_conc_hash[extrinsic_geno_id][condition]
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

    def getTestSuite(self):
        import unittest
        from tests.test_zfin import ZFINTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(ZFINTestCase)

        return test_suite
