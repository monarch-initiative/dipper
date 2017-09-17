import logging
import re
import csv
import gzip
import io
import hashlib
import os

from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.Model import Model
# from dipper.models.assoc.Association import Assoc
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.Environment import Environment
from dipper.utils.DipperUtil import DipperUtil
from dipper import config


logger = logging.getLogger(__name__)


class FlyBase(PostgreSQLSource):
    """
    This is the [Drosophila Genetics](http://www.flybase.org/) resource,
    from which we process genotype and phenotype data about fruitfly.
    Genotypes leverage the GENO genotype model.

    Here, we connect to their public database,
    and download a subset of tables/views to get specifically at
    the geno-pheno data, then iterate over the tables.
    We end up effectively performing joins when adding nodes to the graph.
    We connect using the
    [Direct Chado Access](http://gmod.org/wiki/Public_Chado_Databases#Direct_Chado_Access)

    When running the whole set,
    it performs best by dumping raw triples using the flag ```--format nt```.

    """

    tables = [
        'genotype',             # done
        'feature_genotype',
        'pub',                  # done
        'feature_pub',          # done
        'pub_dbxref',
        'feature_dbxref',
        'cvterm',               # done
        'stock_genotype',       # done
        'stock',                # done
        'organism',
        'organism_dbxref',
        'environment',          # done
        'phenotype',            # done
        'phenstatement',        # done
        'dbxref',               # done
        'phenotype_cvterm',     # done
        'phendesc',             # done
        'environment_cvterm',   # done
        'stockprop'
        # 'feature_cvterm'
        # TODO to get better feature types than (what) is in the feature table
        # itself  (watch out for is_not)
    ]

    resources = [
        {
            'query': '../../resources/sql/fb/feature_relationship.sql',
            'outfile': 'feature_relationship'
        }
    ]

    # columns = { WIP
        #'genotype': (
            #feature_genotype_id, feature_id, genotype_id, chromosome_id, rank,
            #cgroup, cvterm_id),
        #'feature_genotype': (
            #feature_genotype_id, feature_id, genotype_id, chromosome_id,
            #rank, cgroup, cvterm_id),
        #'pub': (
            #pub_id, title, volumetitle, volume, series_name, issue, pyear,
            #pages, miniref, type_id, is_obsolete, publisher, pubplace,
            #uniquename),
        #'feature_pub': (
            #feature_pub_id, feature_id, pub_id),
        #'pub_dbxref': (
            #pub_dbxref_id, pub_id, dbxref_id, is_current),
        #'feature_dbxref': (
            #feature_dbxref_id, feature_id, dbxref_id, is_current),
        #'feature_relationship': (
            #feature_relationship_id, subject_id, object_id, type_id, rank,
            #value),
        #'cvterm': (
            #cvterm_id, cv_id, definition, dbxref_id, is_obsolete,
            #is_relationshiptype, name),
        #'stock_genotype': (
            #stock_genotype_id, stock_id, genotype_id),
        #'stock': (
            #stock_id, dbxref_id, organism_id, name, uniquename, description,
            #type_id, is_obsolete),
        #'organism': (
            #organism_id, abbreviation, genus, species, common_name, comment),
        #'organism_dbxref': (
            #organism_dbxref_id, organism_id, dbxref_id, is_current),
        #'environment': (
            #environment_id, uniquename, description),
        #'phenotype': (
            #phenotype_id, uniquename, observable_id, attr_id, value, cvalue_id,
            #assay_id),
        #'phenstatement': (
            #phenstatement_id, genotype_id, environment_id, phenotype_id,
            #type_id, pub_id),
        #'dbxref': (
            #dbxref_id, db_id, accession, version, description, url),
        #'phenotype_cvterm': (
            #phenotype_cvterm_id, phenotype_id, cvterm_id, rank),
        #'phendesc':  (
            #phendesc_id, genotype_id, environment_id, description, type_id,
            #pub_id),
        #'environment_cvterm': (
            #environment_cvterm_id, environment_id, cvterm_id),
        #'stockprop': (stockprop_id, stock_id, type_id, value, rank)
    #}

    querys = {
        # WIP: start filtering closer to the source,
        # move towards joining there as well.
        # instead of pulling full tables accross the wire
        # then rejoining them here in python, take advantages of
        # the indexes and mature relational engine to the the
        # work on someone elses machine.
        # we can call it "in the cloud"
        # from Lilly: http://gmod.org/wiki/FlyBase_Field_Mapping_Tables

        'feature_dbxref_WIP': """  -- 17M rows in ~2 minutes
            SELECT
            feature.name feature_name, feature.uniquename feature_id,
            organism.abbreviation abbrev, organism.genus, organism.species,
            cvterm.name frature_type, db.name db, dbxref.accession
            FROM feature_dbxref
            JOIN dbxref ON  feature_dbxref.dbxref_id = dbxref.dbxref_id
            JOIN db ON  dbxref.db_id = db.db_id
            JOIN feature ON feature_dbxref.feature_id  = feature.feature_id
            JOIN organism ON feature.organism_id = organism.organism_id
            JOIN cvterm ON feature.type_id = cvterm.cvterm_id
            WHERE feature_dbxref.is_current = true
            AND feature.is_analysis = false
            AND feature.is_obsolete = false
            AND cvterm.is_obsolete = 0
            ;
        """,
        'feature': """
            SELECT feature_id, dbxref_id, organism_id, name, uniquename,
                null as residues, seqlen, md5checksum, type_id, is_analysis,
                timeaccessioned, timelastmodified
            FROM feature WHERE is_analysis = false and is_obsolete = 'f'
        """
    }

    files = {
        'disease_models': {
            'file': 'allele_human_disease_model_data.tsv.gz',
            'url': 'ftp://ftp.flybase.net/releases/current/precomputed_files/human_disease/allele_human_disease_model_data_fb_*.tsv.gz'
        }
    }

    test_keys = {
        'allele': [
            29677937, 23174110, 23230960, 23123654, 23124718, 23146222,
            29677936, 23174703, 11384915, 11397966, 53333044, 23189969,
            3206803, 29677937, 29677934, 23256689, 23213050, 23230614,
            23274987, 53323093, 40362726, 11380755, 11380754, 23121027,
            44425218, 28298666],
        'gene': [
            23220066, 10344219, 58107328, 3132660, 23193483, 3118401, 3128715,
            3128888, 23232298, 23294450, 3128626, 23255338, 8350351, 41994592,
            3128715, 3128432, 3128840, 3128650, 3128654, 3128602, 3165464,
            23235262, 3165510, 3153563, 23225695, 54564652, 3111381, 3111324],
        'annot': [
            437783, 437784, 437785, 437786, 437789, 437796, 459885, 436779,
            436780, 479826],
        'genotype': [
            267393, 267400, 130147, 168516, 111147, 200899, 46696, 328131,
            328132, 328134, 328136, 381024, 267411, 327436, 197293, 373125,
            361163, 403038],
        'feature': [
            11411407, 53361578, 53323094, 40377849, 40362727, 11379415,
            61115970, 11380753, 44425219, 44426878, 44425220],
        'pub': [
            359867, 327373, 153054, 153620, 370777, 154315, 345909, 365672,
            366057, 11380753],
        'strain': [8117, 3649, 64034, 213, 30131],
        'notes': [],
        'organism': [1, 226, 456]
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'flybase')
        logger.setLevel(logging.INFO)
        # to be used to store the version number to be acquired later
        self.version_num = None

        # update the dataset object with details about this resource
        self.dataset = Dataset(
            'flybase', 'FlyBase', 'http://www.flybase.org/', None, None,
            'http://flybase.org/wiki/FlyBase:FilesOverview')

        # source-specific warnings.  will be cleared when resolved.
        logger.warning("we are ignoring normal phenotypes for now")

        # so that we don't have to deal with BNodes, we will create
        # hash lookups for the internal identifiers the hash will
        # hold the type-specific-object-keys to FB public identifiers. Then,
        # subsequent views of the table will lookup the identifiers
        # in the hash.
        # This allows us to do the 'joining' on the fly
        self.idhash = {
            'allele': {}, 'gene': {}, 'publication': {}, 'stock': {},
            'genotype': {}, 'annot': {}, 'notes': {}, 'organism': {},
            'environment': {}, 'feature': {}, 'phenotype': {}, 'cvterm': {},
            'reagent': {}}
        self.dbxrefs = {}
        # to store if a marker is a class or indiv
        self.markers = {'classes': [], 'indiv': []}

        # use this to store internally generated labels for various features
        self.label_hash = {}
        # use this to store the genotype strain ids for genotype labels
        self.geno_bkgd = {}
        # mappings between internal phenotype db key and multiple cv terms
        self.phenocv = {}
        # keep this mapping so we can track fly things to be leaders
        self.feature_to_organism_hash = {}
        # store the feature types, they are needed for making some triples
        self.feature_types = {}
        # when we verify a tax id in eutils
        self.checked_organisms = set()
        self.deprecated_features = set()

        # check to see if there's any ids configured in the config;
        # otherwise, warn
        if 'test_ids' not in config.get_config() \
                or 'disease' not in config.get_config()['test_ids']:
            logger.warning("not configured with disease test ids.")
            self.test_ids = None
        else:
            # select ony those test ids that are omim's.
            self.test_ids = config.get_config()['test_ids']

        return

    def fetch(self, is_dl_forced=False):
        """
        :return:

        """

        # create the connection details for Flybase
        cxn = {
            'host': 'chado.flybase.org', 'database': 'flybase', 'port': 5432,
            'user': 'flybase', 'password': 'no password'}

        self.dataset.setFileAccessUrl(
            ''.join(('jdbc:postgresql://', cxn['host'], ':', str(cxn['port']),
                     '/', cxn['database'])), is_object_literal=True)

        # process the tables
        # self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables, cxn, None, is_dl_forced)

        for query_map in self.resources:
            query_fh = open(os.path.join(
                os.path.dirname(__file__), query_map['query']), 'r')
            query = query_fh.read()
            self.fetch_query_from_pgdb(
                query_map['outfile'], query, None, cxn)

        # we want to fetch the features,
        # but just a subset to reduce the processing time
        # query = \
        #    "SELECT " \
        #    " feature_id, dbxref_id, organism_id, name, uniquename, " \
        #    " null as residues, seqlen, md5checksum, type_id, is_analysis," \
        #    " timeaccessioned, timelastmodified, is_obsolete " \
        #    "FROM feature WHERE is_analysis = false"

        self.fetch_query_from_pgdb(
            'feature', self.querys['feature'], None, cxn, None, is_dl_forced)

        self._get_human_models_file()
        self.get_files(False)
        self.dataset.set_version_by_num(self.version_num)

        return

    def parse(self, limit=None):
        """
        We process each of the postgres tables in turn.
        The order of processing is important here, as we build up a hashmap of
        internal vs external identifers (unique keys by type to FB id).
        These include allele, marker (gene), publication, strain, genotype,
        annotation (association), and descriptive notes.
        :param limit: Only parse this many lines of each table
        :return:

        """
        if limit is not None:
            logger.info("Only parsing first %d rows of each file", limit)
        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        # the following will provide us the hash-lookups
        self._process_dbxref()
        self._process_cvterm()
        self._process_genotypes(limit)
        self._process_pubs(limit)

        # do this before environments to get the external ids
        self._process_environment_cvterm()
        self._process_environments()
        self._process_organisms(limit)  # must be done before features
        self._process_organism_dbxref(limit)
        self._process_features(limit)
        self._process_phenotype(limit)
        self._process_phenotype_cvterm()
        # gets external mappings for features (genes, variants, etc)
        self._process_feature_dbxref(limit)
        # do this after organisms to get the right taxonomy
        self._process_stocks(limit)
        # figures out types of some of the features
        self._get_derived_feature_types(limit)

        # These are the associations amongst the objects above
        self._process_stockprop(limit)
        self._process_pub_dbxref(limit)
        self._process_phendesc(limit)
        self._process_feature_genotype(limit)
        self._process_feature_pub(limit)
        self._process_stock_genotype(limit)
        self._process_phenstatement(limit)   # these are G2P associations

        self._process_feature_relationship(limit)

        self._process_disease_models(limit)
        # TODO add version info from file somehow
        # (in parser rather than during fetching)

        logger.info("Finished parsing.")
        logger.info("Loaded %d nodes", len(self.graph))
        return

    def _process_genotypes(self, limit):
        """
        Add the genotype internal id to flybase mapping to the idhashmap.
        Also, add them as individuals to the graph.

        Triples created:
        <genotype id> a GENO:intrinsic_genotype
        <genotype id> rdfs:label "<gvc> [bkgd]"

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0

        raw = '/'.join((self.rawdir, 'genotype'))
        logger.info("building labels for genotypes")
        geno = Genotype(g)
        fly_tax = 'NCBITaxon:7227'
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (genotype_num, uniquename, description, name) = line

                # if self.testMode is True:
                #     if int(object_key) not in self.test_keys.get('genotype'):
                #         continue

                # add the internal genotype to pub mapping
                genotype_id = 'MONARCH:FBgeno'+str(genotype_num)
                self.idhash['genotype'][genotype_num] = genotype_id

                if description == '':
                    description = None

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    pass
                else:
                    if self.testMode and \
                            int(genotype_num) not in \
                            self.test_keys['genotype']:
                        continue

                    model.addIndividualToGraph(
                        genotype_id, uniquename,
                        Genotype.genoparts['intrinsic_genotype'],
                        description)
                    # we know all genotypes are in flies
                    # FIXME we assume here they are in melanogaster,
                    # but that isn't necessarily true!!!
                    # TODO should the taxon be == genomic background?
                    geno.addTaxon(fly_tax, genotype_id)
                    genotype_iid = self._makeInternalIdentifier(
                        'genotype', genotype_num)
                    model.addComment(
                        genotype_id, genotype_iid)
                    if name.strip() != '':
                        model.addSynonym(genotype_id, name)

        return

    # todo moke singular
    def _process_stocks(self, limit):
        """
        Stock definitions.
        Here we instantiate them as instances of the given taxon.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0

        raw = '/'.join((self.rawdir, 'stock'))
        logger.info("building labels for stocks")

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (stock_id, dbxref_id, organism_id, name, uniquename,
                 description, type_id, is_obsolete) = line
# 2       12153979        1       2       FBst0000002     w[*]; betaTub60D[2] Kr[If-1]/CyO        10670

                stock_num = stock_id
                stock_id = 'FlyBase:'+uniquename
                self.idhash['stock'][stock_num] = stock_id
                stock_label = description

                organism_key = organism_id
                taxon = self.idhash['organism'][organism_key]

                # from what i can tell, the dbxrefs are just more FBst,
                # so no added information vs uniquename

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    pass
                else:
                    if self.testMode \
                            and int(stock_num) not in self.test_keys['strain']:
                        continue

                    # tax_label = self.label_hash[taxon]  # unused
                    # add the tax in case it hasn't been already
                    model.addClassToGraph(taxon)
                    model.addIndividualToGraph(stock_id, stock_label, taxon)
                    if is_obsolete == 't':
                        model.addDeprecatedIndividual(stock_id)

        return

    # todo make singular
    def _process_pubs(self, limit):
        """
        Flybase publications.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0

        raw = '/'.join((self.rawdir, 'pub'))
        logger.info("building labels for pubs")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                (pub_id, title, volumetitle, volume, series_name, issue, pyear,
                 pages, miniref, type_id, is_obsolete, publisher, pubplace,
                 uniquename) = line
# 2       12153979        1       2       FBst0000002     w[*]; betaTub60D[2] Kr[If-1]/CyO        10670
                # if self.testMode is True:
                #     if int(object_key) not in self.test_keys.get('genotype'):
                #         continue

                pub_num = pub_id
                pub_id = 'FlyBase:'+uniquename.strip()
                self.idhash['publication'][pub_num] = pub_id

                # TODO figure out the type of pub by type_id
                if not re.match(r'(FBrf|multi)', uniquename):
                    continue
                line_counter += 1

                reference = Reference(g, pub_id)
                if title != '':
                    reference.setTitle(title)
                if pyear != '':
                    reference.setYear(str(pyear))
                if miniref != '':
                    reference.setShortCitation(miniref)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    pass
                else:
                    if self.testMode \
                            and int(pub_num) not in self.test_keys['pub']:
                        continue

                    if is_obsolete == 't':
                        model.addDeprecatedIndividual(pub_id)
                    else:
                        reference.addRefToGraph()
        return

    # todo make singular
    def _process_environments(self):
        """
        There's only about 30 environments in which the phenotypes
        are recorded.
        There are no externally accessible identifiers for environments,
        so we make anonymous nodes for now.
        Some of the environments are comprised of >1 of the other environments;
        we do some simple parsing to match the strings of the environmental
        labels to the other atomic components.

        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        raw = '/'.join((self.rawdir, 'environment'))
        logger.info("building labels for environment")
        env_parts = {}
        label_map = {}
        env = Environment(g)
        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (environment_id, uniquename, description) = line
                # 22      heat sensitive | tetracycline conditional

                environment_num = environment_id
                environment_internal_id = self._makeInternalIdentifier(
                    'environment', environment_num)
                if environment_num not in self.idhash['environment']:
                    self.idhash['environment'][environment_num] = \
                        environment_internal_id

                environment_id = self.idhash['environment'][environment_num]
                environment_label = uniquename
                if environment_label == 'unspecified':
                    environment_label += ' environment'
                env.addEnvironment(environment_id, environment_label)
                self.label_hash[environment_id] = environment_label

                # split up the environment into parts
                # if there's parts, then add them to the hash;
                # we'll match the components in a second pass
                components = re.split(r'\|', uniquename)
                if len(components) > 1:
                    env_parts[environment_id] = components
                else:
                    label_map[environment_label] = environment_id

            # ### end loop through file

        # build the environmental components
        for eid in env_parts:
            eid = eid.strip()
            for e in env_parts[eid]:
                # search for the environmental component by label
                env_id = label_map.get(e.strip())
                env.addComponentToEnvironment(eid, env_id)

        return

    # todo make singular
    def _process_features(self, limit):
        """
        These are all of the genomic features genes, variations,
        transgenes, etc
        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'feature'))
        logger.info("building labels for features")

        line_counter = 0
        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (feature_id, dbxref_id, organism_id, name, uniquename,
                 residues, seqlen, md5checksum, type_id, is_analysis,
                 timeaccessioned, timelastmodified) = line

                feature_key = feature_id
                if re.search(r'[\|\s\[\]\{\}\\<\>]', uniquename):
                    # some uniquenames have pipes or other nasty chars!
                    # for example: FB||||FBrf0133242|Hugh-u1
                    feature_id = self._makeInternalIdentifier(
                        'feature', feature_key)
                else:
                    feature_id = 'FlyBase:'+uniquename
                self.idhash['feature'][feature_key] = feature_id
                self.feature_types[feature_key] = type_id
                self.label_hash[feature_id] = name

                if feature_key not in self.feature_to_organism_hash:
                    self.feature_to_organism_hash[feature_key] = set()
                self.feature_to_organism_hash[feature_key].add(organism_id)

                # HACK - FBgn are genes, and therefore classes,
                # all else be individuals
                is_gene = False
                if re.search(r'(FBgn|FBog)', feature_id):
                    self.idhash['gene'][feature_key] = feature_id
                    is_gene = True
                elif re.search(r'FBa[lb]', feature_id):
                    self.idhash['allele'][feature_key] = feature_id
                elif re.search(r'FBt[ip]', feature_id):
                    self.idhash['feature'][feature_key] = feature_id

                if self.testMode and \
                        int(feature_key) not in self.test_keys['gene'] + \
                        self.test_keys['allele'] + self.test_keys['feature']:
                    continue

                # now do something with it!
                # switch on type_id
                if name.strip() == '':
                    name = uniquename

                type_key = type_id
                type_id = self.idhash['cvterm'][type_key]

                # skip some features by type
                types_to_skip = [
                    'SO:0000316',  # CDS
                    'SO:0000696',  # oligos
                    'SO:0000358',  # polypeptide
                    'SO:0000234',  # transcripts
                    ]

                type_keys_to_skip = [
                    596,    # pcr_product
                    57096,  # mature peptide
                    57097,  # signal_peptide
                    57270,  # repeat masker
                    58210,  # alignment
                    59643,  # cDNA_clone
                    60006,  # uncharacterized_change_in_nucleotide_sequence
                    61351,  # oligo
                    61467,  # polypeptide_domain
                    257,    # exon
                    286,    # intron
                ]

                organisms_to_skip = [
                    2  # computational result
                ]

                if type_id in types_to_skip \
                        or int(type_key) in type_keys_to_skip\
                        or int(organism_id) in organisms_to_skip:
                    continue

                line_counter += 1

                if int(type_key) == 604:  # RNAi_reagent
                    # TODO add other reagents?
                    self.idhash['reagent'][feature_key] = feature_id

                # deal with the taxonomy
                # only get taxa for features that are actually used in our set
                tax_internal_id = self._makeInternalIdentifier(
                    'organism', organism_id)
                if organism_id not in self.checked_organisms:
                    # will get the NCBITax if necessary
                    tax_id = self._get_organism_id(organism_id)
                    self.checked_organisms.add(organism_id)
                else:
                    tax_id = self.idhash['organism'][organism_id]

                tax_label = self.label_hash.get(tax_id)
                if not re.search(r'FBog', feature_id) \
                        and re.search(r'Drosophila', tax_label):
                    # make only fly things leaders
                    model.makeLeader(feature_id)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    pass
                else:
                    if is_gene:
                        model.addClassToGraph(
                            feature_id, name, type_id)
                        g.addTriple(
                            feature_id, model.object_properties['in_taxon'],
                            tax_id)
                    else:
                        if re.search('FBa[lb]', feature_id):
                            type_id = Genotype.genoparts['allele']
                        model.addIndividualToGraph(feature_id, name, type_id)

                    # stop adding what we do not appreciate
                    # if is_obsolete == 't':
                    #    if is_gene:
                    #        model.addDeprecatedClass(feature_id)
                    #    else:
                    #        model.addDeprecatedIndividual(feature_id)
                    #    self.deprecated_features.add(feature_key)

                    model.addClassToGraph(tax_id)
                    if tax_id != tax_internal_id:
                        model.addEquivalentClass(tax_id, tax_internal_id)

                    model.addComment(
                        feature_id,
                        self._makeInternalIdentifier('feature', feature_key))

            # TODO save checked_organisms fbid to ncbitax mapping to
            # a local file to speed up subsequent searches

        return

    def _process_feature_genotype(self, limit):

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        raw = '/'.join((self.rawdir, 'feature_genotype'))
        logger.info("processing genotype features")
        geno = Genotype(g)
        line_counter = 0

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                (feature_genotype_id, feature_id, genotype_id, chromosome_id,
                 rank, cgroup, cvterm_id) = line
                # 1	23273518	2	23159230	0	0	60468

                feature_key = feature_id
                if feature_key not in self.idhash['feature']:
                    continue
                feature_id = self.idhash['feature'][feature_key]

                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]

                if self.testMode and not (
                        int(feature_key) in
                        self.test_keys['gene']+self.test_keys['allele'] and
                        int(genotype_key) in self.test_keys['genotype']):
                    continue

                # what is cvterm_id for in this context???
                # cgroup is the order of composition of things in
                # the genotype label (complementation group?).
                # rank is the order that they appear in the label
                # sometimes the same feature is listed twice;
                # not sure if this is a mistake, or zygosity, or?
                if feature_id is not None and genotype_id is not None:
                    geno.addParts(
                        feature_id, genotype_id,
                        geno.object_properties['has_alternate_part'])

                # TODO we will build up the genotypes here... lots to do

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def _process_phendesc(self, limit):
        """
        The description of the resulting phenotypes
        with the genotype+environment

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'phendesc'))
        logger.info("processing G2P")

        line_counter = 0
        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (phendesc_id, genotype_id, environment_id, description,
                 type_id, pub_id) = line
                # 1	2	1	Hemizygous males are wild type, homozygous males are sterile.	60466	209729

                line_counter += 1
                phendesc_key = phendesc_id
                phendesc_id = self._makeInternalIdentifier(
                    'phendesc', phendesc_key)

                # for now, just attach the description to the genotype
                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]
                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                environment_key = environment_id
                environment_id = self.idhash['environment'][environment_key]

                if self.testMode and\
                        int(genotype_key) not in self.test_keys['genotype']:
                    continue

                # TODO type id ==> ECO???

                # just make associations with abnormal phenotype
                phenotype_id = 'FBcv:0001347'
                assoc = G2PAssoc(g, self.name, genotype_id, phenotype_id)
                assoc.add_source(pub_id)
                assoc.set_description(description)
                assoc.set_environment(environment_id)
                assoc.add_association_to_graph()
                assoc_id = assoc.get_association_id()
                model.addComment(assoc_id, phendesc_id)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def _process_feature_pub(self, limit):
        """
        The description of the resulting phenotypes
        with the genotype+environment

        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'feature_pub'))
        logger.info("processing feature_pub")

        line_counter = 0

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (feature_pub_id, feature_id, pub_id) = line
                # 1440    3175682 62137
                # 2       3160606 99159

                feature_key = feature_id
                if self.testMode and not (
                        int(feature_key) in
                        self.test_keys['gene']+self.test_keys['allele']and
                        int(pub_id) in self.test_keys['pub']):
                    continue
                if feature_key not in self.idhash['feature']:
                    continue
                feature_id = self.idhash['feature'][feature_key]
                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                g.addTriple(
                    pub_id, model.object_properties['mentions'], feature_id)

                line_counter += 1

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def _process_stock_genotype(self, limit):
        """
        The genotypes of the stocks.

        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        # model = Model(g)  # unused
        raw = '/'.join((self.rawdir, 'stock_genotype'))
        logger.info("processing stock genotype")
        geno = Genotype(g)
        line_counter = 0

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (stock_genotype_id, stock_id, genotype_id) = line

                stock_key = stock_id
                stock_id = self.idhash['stock'][stock_key]
                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]

                if self.testMode and int(genotype_key) not in \
                        self.test_keys['genotype']:
                    continue

                g.addTriple(
                    stock_id, geno.object_properties['has_genotype'],
                    genotype_id)

                line_counter += 1

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def _process_pub_dbxref(self, limit):
        """
        Xrefs for publications (ie FBrf = PMID)
        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'pub_dbxref'))
        logger.info("processing pub_dbxref")

        line_counter = 0

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (pub_dbxref_id, pub_id, dbxref_id, is_current) = line
                # 49648	43222	395730	t

                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                if self.testMode and int(pub_key) not in self.test_keys['pub']:
                    continue

                # get any dbxrefs for pubs, including pmids and dois
                dbxref_key = dbxref_id
                if str(dbxref_key) in self.dbxrefs:
                    dbxrefs = self.dbxrefs[str(dbxref_key)]
                    # pub_dbs = [75, 51, 76, 95, 126]
                    pmid_ids = [50, 77, 275, 286, 347]
                    # flybase_ids = [4]  # TODO unused
                    isbn = [75, 51]
                    for d in dbxrefs:
                        dbxref_id = None
                        if int(d) in pmid_ids:
                            if re.match(r'^PMID', dbxrefs[d]):
                                dbxref_id = dbxrefs[d].strip()
                            else:
                                dbxref_id = 'PMID:'+dbxrefs[d].strip()
                            model.makeLeader(dbxref_id)
                        elif int(d) in isbn:
                            dbxref_id = 'ISBN:'+dbxrefs[d].strip()
                        elif int(d) == 161:
                            dbxref_id = 'DOI:'+dbxrefs[d].strip()
                        # elif int(d) == 4:
                        #     dbxref_id = 'FlyBase:'+dbxrefs[d].strip()

                        if dbxref_id is not None:
                            reference = Reference(
                                g, dbxref_id,
                                Reference.ref_types['publication'])
                            reference.addRefToGraph()
                            model.addSameIndividual(pub_id, dbxref_id)
                            line_counter += 1

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def _process_dbxref(self):
        """
        We bring in the dbxref identifiers and store them in a hashmap for
        lookup in other functions.
        Note that some dbxrefs aren't mapped to identifiers.
        For example, 5004018 is mapped to a string,
            "endosome & imaginal disc epithelial cell | somatic clone..."
        In those cases, there just isn't a dbxref that's used
        when referencing with a cvterm; it'll just use the internal key.

        :return:

        """

        raw = '/'.join((self.rawdir, 'dbxref'))
        logger.info("processing dbxrefs")
        line_counter = 0

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (dbxref_id, db_id, accession, version, description, url) = line
                # dbxref_id	db_id	accession	version	description	url
                # 1	2	SO:0000000	""

                db_ids = {          # the databases to fetch
                    50: 'PMID',     # pubmed
                    68: 'RO',       # obo-rel
                    71: 'FBdv',     # FBdv
                    74: 'FBbt',     # FBbt
                    # 28:,          # genbank
                    30: 'OMIM',     # MIM
                    # 38,           # ncbi
                    75: 'ISBN',     # ISBN
                    46: 'PMID',     # PUBMED
                    51: 'ISBN',     # isbn
                    52: 'SO',       # so
                    # 76,           # http
                    77: 'PMID',     # PMID
                    80: 'FBcv',     # FBcv
                    # 95,           # MEDLINE
                    98: 'REACT',    # Reactome
                    103: 'CHEBI',   # Chebi
                    102: 'MESH',    # MeSH
                    106: 'OMIM',    # OMIM
                    105: 'KEGG-path',  # KEGG pathway
                    107: 'DOI',     # doi
                    108: 'CL',      # CL
                    114: 'CHEBI',   # CHEBI
                    115: 'KEGG',    # KEGG
                    116: 'PubChem',  # PubChem
                    # 120,          # MA???
                    3: 'GO',        # GO
                    4: 'FlyBase',   # FlyBase
                    # 126,          # URL
                    128: 'PATO',    # PATO
                    # 131,          # IMG
                    2: 'SO',        # SO
                    136: 'MESH',    # MESH
                    139: 'CARO',    # CARO
                    140: 'NCBITaxon',  # NCBITaxon
                    # 151,          # MP  ???
                    161: 'DOI',     # doi
                    36: 'BDGP',     # BDGP
                    # 55,           # DGRC
                    # 54,           # DRSC
                    # 169,          # Transgenic RNAi project???
                    231: 'RO',      # RO ???
                    180: 'NCBIGene',  # entrezgene
                    # 192,          # Bloomington stock center
                    197: 'UBERON',  # Uberon
                    212: 'ENSEMBL',  # Ensembl
                    # 129,          # GenomeRNAi
                    275: 'PMID',    # PubMed
                    286: 'PMID',    # pmid
                    264: 'HGNC',
                    # 265: 'OMIM',  # OMIM_Gene
                    266: 'OMIM',    # OMIM_Phenotype
                    300: 'DOID',    # DOID
                    302: 'MESH',    # MSH
                    347: 'PMID',    # Pubmed
                }

                if accession.strip() != '' and int(db_id) in db_ids:
                    # scrub some identifiers here
                    m = re.match(
                        r'(doi|SO|GO|FBcv|FBbt_root|FBdv|FBgn|FBdv_root|FlyBase|FBbt):',
                        accession)
                    if m:
                        accession = re.sub(m.group(1)+r'\:', '', accession)
                    elif re.match(
                            r'(FlyBase miscellaneous CV|cell_lineprop|relationship type|FBgn$)',
                            accession):
                        continue
                    elif re.match(r'\:', accession):  # starts with a colon
                        accession = re.sub(r'\:', '', accession)
                    elif re.search(r'\s', accession):
                        # skip anything with a space
                        # logger.debug(
                        #   'dbxref %s accession has a space: %s',
                        #   dbxref_id, accession)
                        continue

                    if re.match(r'http', accession):
                        did = accession.strip()
                    else:
                        prefix = db_ids.get(int(db_id))
                        did = ':'.join((prefix, accession.strip()))
                        if re.search(r'\:', accession) and prefix != 'DOI':
                            logger.warning(
                                'id %s may be malformed; skipping', did)

                    self.dbxrefs[dbxref_id] = {db_id: did}

                elif url != '':
                    self.dbxrefs[dbxref_id] = {db_id: url.strip()}
                else:
                    continue

                # the following are some special cases that we scrub
                if int(db_id) == 2 \
                        and accession.strip() == 'transgenic_transposon':
                    # transgenic_transposable_element
                    self.dbxrefs[dbxref_id] = {db_id: 'SO:0000796'}

                line_counter += 1

        return

    def _process_phenotype(self, limit):
        """
        Get the phenotypes, and declare the classes.
        If the "observable" is "unspecified", then we assign the phenotype to
        the "cvalue" id; otherwise we convert the phenotype into a
        uberpheno-style identifier, simply based on the anatomical part that's
        affected...that is listed as the observable_id, concatenated with
        the literal "PHENOTYPE"

        Note that some of the phenotypes no not have a dbxref to a FBcv;
        for these cases it will make a node with an anonymous node with an
        internal id like,  "_fbcvtermkey100920PHENOTYPE".  This is awkward,
        but not sure how else to construct identifiers.
        Maybe they should be fed back into Upheno and then leveraged by FB?

        Note that assay_id is the same for all current items,
        so we do nothing with this.
        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'phenotype'))
        logger.info("processing phenotype")

        line_counter = 0

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (phenotype_id, uniquename, observable_id, attr_id, value,
                 cvalue_id, assay_id) = line

                # 8505	unspecified
                # 20142	mesothoracic leg disc | somatic clone 87719 60468  60468 60468
                # 8507	sex comb | ectopic 88877 60468  60468 60468
                # 8508	tarsal segment	83664 60468  60468 60468
                # 18404	oocyte | oogenesis stage S9	86769 60468  60468 60468
                # for now make these as phenotypic classes
                # will need to xref at some point
                phenotype_key = phenotype_id
                phenotype_id = None
                phenotype_internal_id = self._makeInternalIdentifier(
                    'phenotype', phenotype_key)
                phenotype_label = None
                self.label_hash[phenotype_internal_id] = uniquename
                cvterm_id = None
                if observable_id != '' \
                        and int(observable_id) == 60468:
                    # undefined - typically these are already phenotypes
                    if cvalue_id in self.idhash['cvterm']:
                        cvterm_id = self.idhash['cvterm'][cvalue_id]
                        phenotype_id = self.idhash['cvterm'][cvalue_id]
                elif observable_id in self.idhash['cvterm']:
                    # observations to anatomical classes
                    cvterm_id = self.idhash['cvterm'][observable_id]
                    phenotype_id = \
                        self.idhash['cvterm'][observable_id] + 'PHENOTYPE'
                    if cvterm_id is not None and cvterm_id in self.label_hash:
                        phenotype_label = self.label_hash[cvterm_id]
                        phenotype_label += ' phenotype'
                        self.label_hash[phenotype_id] = phenotype_label
                    else:
                        logger.info('cvtermid=%s not in label_hash', cvterm_id)

                else:
                    logger.info(
                        "No observable id or label for %s: %s",
                        phenotype_key, uniquename)

                # TODO store this composite phenotype in some way
                # as a proper class definition?
                self.idhash['phenotype'][phenotype_key] = phenotype_id

                # assay_id is currently only "undefined" key=60468

                if not self.testMode and\
                        limit is not None and line_counter > limit:
                    pass
                else:
                    if phenotype_id is not None:
                        # assume that these fit into the phenotypic uberpheno
                        # elsewhere
                        model.addClassToGraph(phenotype_id, phenotype_label)
                        line_counter += 1

        return

    def _process_phenstatement(self, limit):
        """
        The phenstatements are the genotype-to-phenotype associations,
        in the context of an environment.
        These are also curated to a publication. So we make oban associations,
        adding the pubs as a source.  We additionally add the internal key as
        a comment for tracking purposes.
        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'phenstatement'))
        logger.info("processing phenstatement")

        line_counter = 0

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (phenstatement_id, genotype_id, environment_id, phenotype_id,
                 type_id, pub_id) = line

                # 168549	166695	1	8507	60468	151256
                # 168550	166695	1	8508	60468	151256
                # 168551	166696	1	8509	60468	151256
                # 168552	166696	1	8510	60468	151256
                line_counter += 1
                phenstatement_key = phenstatement_id
                phenstatement_id = self._makeInternalIdentifier(
                    'phenstatement', phenstatement_key)
                genotype_key = genotype_id

                if self.testMode and \
                        int(genotype_key) not in self.test_keys['genotype']:
                    continue

                genotype_id = self.idhash['genotype'][genotype_key]
                environment_key = environment_id
                environment_id = self.idhash['environment'][environment_key]
                phenotype_key = phenotype_id
                phenotype_internal_id = self._makeInternalIdentifier(
                    'phenotype', phenotype_key)  # TEMP
                phenotype_internal_label = self.label_hash[
                    phenotype_internal_id]
                phenotype_id = self.idhash['phenotype'][phenotype_key]
                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                # figure out if there is a relevant stage

                assoc = G2PAssoc(g, self.name, genotype_id, phenotype_id)
                if phenotype_id in self.phenocv:
                    stages = set(
                        s for s in
                        self.phenocv[phenotype_id] if re.match(r'FBdv', s))
                    if len(stages) == 1:
                        s = stages.pop()
                        assoc.set_stage(s, s)
                    elif len(stages) > 1:
                        logger.warning(
                            "There's more than one stage specified per " +
                            "phenotype. I don't know what to do. %s",
                            str(stages))
                    non_stage_ids = self.phenocv[phenotype_id] - stages
                    logger.debug(
                        'Other non-stage bits: %s', str(non_stage_ids))
                    # TODO do something with the other parts
                    # of a pheno-cv relationship
                assoc.set_environment(environment_id)
                # TODO check remove unspecified environments?
                assoc.add_source(pub_id)
                assoc.add_association_to_graph()
                assoc_id = assoc.get_association_id()
                model.addComment(assoc_id, phenstatement_id)
                model.addDescription(assoc_id, phenotype_internal_label)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

        return

    def _process_phenotype_cvterm(self):
        """
        These are the qualifiers for the phenotype location itself.
        But are just the qualifiers.
        The actual "observable" part of the phenotype is only in
        the phenotype table. These get added to a lookup variable used to
        augment a phenotype association statement.
        :return:

        """

        line_counter = 0
        raw = '/'.join((self.rawdir, 'phenotype_cvterm'))
        logger.info("processing phenotype cvterm mappings")

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (phenotype_cvterm_id, phenotype_id, cvterm_id, rank) = line

                # 4532	8507	60793	0
                # 4533	8513	60830	0

                # add the internal genotype to pub mapping
                phenotype_key = phenotype_id
                cvterm_key = cvterm_id
                phenotype_id = self.idhash['phenotype'][phenotype_key]
                if cvterm_key in self.idhash['cvterm']:
                    cvterm_id = self.idhash['cvterm'][cvterm_key]
                    if phenotype_key not in self.phenocv:
                        self.phenocv[phenotype_id] = set()
                    self.phenocv[phenotype_id].add(cvterm_id)
                else:
                    logger.info(
                        "Not storing the cvterm info for %s", cvterm_key)

        return

    def _process_cvterm(self):
        """
        CVterms are the internal identifiers for any controlled vocab
        or ontology term.  Many are xrefd to actual ontologies.  The actual
        external id is stored in the dbxref table, which we place into
        the internal hashmap for lookup with the cvterm id.  The name of
        the external term is stored in the "name" element of this table, and
        we add that to the label hashmap for lookup elsewhere

        :return:

        """

        line_counter = 0
        raw = '/'.join((self.rawdir, 'cvterm'))
        logger.info("processing cvterms")

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (cvterm_id, cv_id, definition, dbxref_id, is_obsolete,
                 is_relationshiptype, name) = line

                # 316 6 1665919 0 0 rRNA_cleavage_snoRNA_primary_transcript
                # 28  5 1663309 0 0 synonym
                # 455 6 1665920 0 0 tmRNA

                # not sure the following is necessary
                # cv_prefixes = {
                #     6 : 'SO',
                #     20: 'FBcv',
                #     28: 'GO',
                #     29: 'GO',
                #     30: 'GO',
                #     31: 'FBcv',  # not actually FBcv - I think FBbt.
                #     32: 'FBdv',
                #     37: 'GO',   # these are relationships
                #     73: 'DOID'
                # }

                # if int(cv_id) not in cv_prefixes:
                #     continue
                cvterm_key = cvterm_id
                cvterm_id = self._makeInternalIdentifier('cvterm', cvterm_key)
                self.label_hash[cvterm_id] = name
                self.idhash['cvterm'][cvterm_key] = cvterm_id
                # look up the dbxref_id for the cvterm
                # hopefully it's one-to-one
                dbxrefs = self.dbxrefs.get(dbxref_id)
                if dbxrefs is not None:
                    if len(dbxrefs) > 1:
                        logger.info(
                            ">1 dbxref for this cvterm (%s: %s): %s",
                            str(cvterm_id), name, dbxrefs.values())
                    elif len(dbxrefs) == 1:
                        # replace the cvterm with
                        # the dbxref (external) identifier
                        did = dbxrefs.popitem()[1]
                        # get the value
                        self.idhash['cvterm'][cvterm_key] = did
                        # also add the label to the dbxref
                        self.label_hash[did] = name
        return

    def _process_environment_cvterm(self):
        """
        This is the mapping between the internal environment id
        and the external ones; here we map the internal environment id to
        the external one in the hashmap.
        :return:

        """

        line_counter = 0
        raw = '/'.join((self.rawdir, 'environment_cvterm'))
        logger.info("processing environment to cvterm mappings")

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (environment_cvterm_id, environment_id, cvterm_id) = line
                # 1	1	60468

                environment_key = environment_id

                cvterm_key = cvterm_id
                cvterm_id = self.idhash['cvterm'][cvterm_key]

                # look up the dbxref_id for the cvterm
                # hopefully it's one-to-one
                self.idhash['environment'][environment_key] = cvterm_id

        return

    def _process_feature_dbxref(self, limit):
        """
        This is the mapping between the flybase features and external
        repositories. Generally we want to leave the flybase feature id
        as the primary identifier. But we need to make the equivalences/sameAs.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'feature_dbxref'))
        logger.info("processing feature_dbxref mappings")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:

                (feature_dbxref_id, feature_id, dbxref_id, is_current) = line

                # 431890	3091292	596211	t
                # 2	9	55044	t
                # 3	9	55045	t
                # 437595	4551668	277309	t
                # 437596	4551662	277307	t

                if is_current == 'f':
                    # not sure what to do with it?
                    continue

                feature_key = feature_id

                if self.testMode \
                        and int(feature_key) not in \
                        self.test_keys['gene'] + self.test_keys['allele']:
                    continue

                if feature_key not in self.idhash['feature']:
                    # some features may not be found in the hash
                    # if they are "analysis features"
                    # logger.debug("Feature %s not found in hash", feature_key)
                    continue
                feature_id = self.idhash['feature'][feature_key]
                dbxref_key = dbxref_id
                dbxrefs = self.dbxrefs.get(dbxref_key)

                if dbxrefs is not None:
                    for d in dbxrefs:
                        # need to filter based on db ?
                        # TODO make other species' identifiers primary??
                        # instead of flybase?
                        did = dbxrefs.get(d)
                        if did.endswith('&class=protein'):
                            did = did[0:len(dbxrefs)-15]
                        # don't make something sameAs itself
                        if did == feature_id:
                            continue
                        dlabel = self.label_hash.get(did)
                        if re.search(r'FB(gn|og)', feature_id):
                            # only want to add equivalences for fly things
                            if not re.match(r'OMIM', did):
                                # these are only omim diseases, not genes;
                                # we shouldn't be adding these here anyway
                                # model.addClassToGraph(did, dlabel)
                                # model.addXref(feature_id, did)
                                True  # that
                        elif did is not None and dlabel is not None \
                                and feature_id is not None:
                            model.addIndividualToGraph(did, dlabel)
                            model.addXref(feature_id, did)
                        line_counter += 1

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break

                # FIXME - some flybase genes are xrefed to OMIM diseases!!!!!!
                # for example,
                # FBog0000375495 xref to omim 601181 (gene)
                # and 608033 (phenotype)

        return

    def _get_derived_feature_types(self, limit):
        """
        Make a pass through the feature table in order to properly type
        the FBal (allele) features, which are derived either from other
        sequence features (which can be things like RNAi products)
        or transgenic-transposons.  We'll save the allele type into a hasmap.

        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'feature_relationship'))
        logger.info("determining some feature types based on relationships")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                (feature_relationship_id, subject_id, object_id, name, rank,
                 value) = line

                if name == 'derived_tp_assoc_alleles':
                    # derived_tp_assoc_alleles
                    self.feature_types[subject_id] = \
                        Genotype.genoparts['transgenic_insertion']
                    sid = self.idhash['allele'].get(subject_id)
                    model.addType(sid, self.feature_types[subject_id])
                elif name == 'derived_sf_assoc_alleles':
                    # only take the derived_sf_assoc_alleles
                    # my subject is a reagent_targeted_gene
                    # my object is the dsRNA
                    self.feature_types[subject_id] = \
                        Genotype.genoparts['reagent_targeted_gene']
                    sid = self.idhash['allele'].get(subject_id)
                    model.addType(sid, self.feature_types[subject_id])

                else:
                    continue

        return

    def _process_feature_relationship(self, limit):
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        # ti_allele_map = {}  # TODO to be used when building genotypes

        line_counter = 0
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, 'feature_relationship'))
        logger.info("processing feature relationships")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                (feature_relationship_id, subject_id, object_id, name, rank,
                 value) = line
                # 7253191 11713123        3177614 27      0
                # 18513040        23683101        11507545        26      0
                # 7130199 9068909 11507822        26      0
                # 18513041        23683101        11507448        26      0
                # 7130197 9346315 11507821        26      0

                if self.testMode and not (
                        int(subject_id) in self.test_keys['gene'] +
                        self.test_keys['allele']+self.test_keys['feature'] and
                        int(object_id) in self.test_keys['gene'] +
                        self.test_keys['allele'] +
                        self.test_keys['feature']):
                    continue

                if subject_id in self.deprecated_features \
                        or object_id in self.deprecated_features:
                    logger.debug(
                        "Skipping deprecated feature_relationship %s %s",
                        str(subject_id), str(object_id))
                    continue

                # TODO move this out of the if later
                # allele of gene
                # in sql, we limited the
                # subject to type_id = 219,33 object type_id  219  #??? TEC
                # subject = variation
                # object = gene
                if name in [
                        'alleleof',
                        'molec_dups',
                        'molec_deletes',
                        'molec_partdeletes',
                        'molec_partdups',
                        'useful_Df_from_cyto',
                        'useful_Df_direct',
                        'useful_Dp_direct',
                        'useful_Dp_from_cyto',
                        'deletes',
                        'part_deletes',
                        'duplicates',
                        'part_duplicates']:
                    allele_id = None
                    gene_id = None
                    if subject_id in self.idhash['allele']:
                        allele_id = self.idhash['allele'][subject_id]
                    if object_id in self.idhash['gene']:
                        gene_id = self.idhash['gene'][object_id]
                    if gene_id is not None and gene_id in self.label_hash:
                        # TODO FAIL: KeyError: None   default?
                        logger.info("getting label for gene_id:\t%s", gene_id)
                        gene_label = self.label_hash[gene_id]
                    else:
                        if gene_id is None:
                            logger.error(
                                "The gene_id for object_id is None: %s \t %s",
                                str(subject_id), str(object_id))
                        if not gene_id not in self.label_hash:
                            logger.error(
                                "gene_id's label missing for: %s\t%s\t%s",
                                str(subject_id), str(object_id), str(object_id))
                        continue
                    # TODO move this out of the if later
                    line_counter += 1
                    if allele_id is not None and gene_id is not None:
                        if self.feature_types[subject_id] == \
                                Genotype.genoparts['reagent_targeted_gene']:
                            g.addTriple(
                                allele_id,
                                Genotype.object_properties[
                                    'is_targeted_expression_variant_of'],
                                gene_id)
                        elif self.feature_types[subject_id] == \
                                Genotype.genoparts['transgenic_insertion']:
                            geno.addSequenceDerivesFrom(allele_id, gene_id)
                        elif re.match(r'\w+\\', gene_label):
                            # the thing that this is the allele of,
                            # is in some other species,
                            # so we do not want to create
                            # the is_allele_of between then.
                            geno.addSequenceDerivesFrom(allele_id, gene_id)
                        else:
                            # assume that the gene is in the same species
                            geno.addAlleleOfGene(allele_id, gene_id)
                    else:
                        if allele_id is None \
                                and subject_id in self.idhash['feature']:
                            feature_id = self.idhash['feature'][subject_id]
                            logger.debug(
                                "this thing %s is not an allele", feature_id)
                        if gene_id is None \
                                and subject_id in self.idhash['feature']:
                            feature_id = self.idhash['feature'][subject_id]
                            logger.debug(
                                "this thing %s is not a gene", feature_id)
                elif name == 'associated_with':

                    allele_id = None
                    gene_id = None
                    reagent_id = None
                    ti_id = None

                    if object_id in self.idhash['allele']:
                        allele_id = self.idhash['allele'][object_id]
                    elif object_id in self.idhash['reagent']:
                        reagent_id = self.idhash['reagent'][object_id]
                    elif object_id in self.idhash['feature']:
                        of = self.idhash['feature'][object_id]
                        if re.search(r'FBt[ip]', of):
                            ti_id = of

                    if object_id in self.idhash['gene']:
                        gene_id = self.idhash['gene'][object_id]

                    if subject_id in self.idhash['gene']:
                        gene_id = self.idhash['gene'][subject_id]
                    elif subject_id in self.idhash['reagent']:
                        reagent_id = self.idhash['reagent'][subject_id]
                    elif subject_id in self.idhash['allele']:
                        allele_id = self.idhash['allele'][subject_id]

                    if allele_id is not None and gene_id is not None:
                        geno.addAlleleOfGene(allele_id, gene_id)
                    elif reagent_id is not None and gene_id is not None:
                        reagent_label = self.label_hash.get(reagent_id)
                        geno.addGeneTargetingReagent(
                            reagent_id, reagent_label, None, gene_id)
                        geno.addReagentTargetedGene(reagent_id, gene_id)
                    elif reagent_id is not None and allele_id is not None:
                        geno.addReagentTargetedGene(
                                reagent_id, None, allele_id)
                    # elif allele_id is not None and ti_id is not None:
                    #     # the FBti == transgenic insertion,
                    #     which is basically the sequence alteration
                    #     geno.addParts(
                    #       ti_id, allele_id,
                    #       geno.object_properties['has_alternate_part'])
                    elif reagent_id is not None and ti_id is not None:
                        g.addTriple(
                            ti_id, geno.object_properties['targeted_by'],
                            reagent_id)

                # derived_tp_assoc_alleles
                elif name == 'derived_tp_assoc_alleles':
                    # note that the internal type id changed
                    # from 129784 --> 133526 around 02.2016
                    # note that this relationship is only specified between
                    # an allele and a tp. therefore we know the FBal should be
                    # a transgenic_insertion
                    allele_id = None

                    if subject_id in self.idhash['allele']:
                        allele_id = self.idhash['allele'][subject_id]

                    if object_id not in self.idhash['feature']:
                        continue
                    tp_id = self.idhash['feature'][object_id]
                    # if allele_id is not None and tp_id is not None:
                    #     geno.addParts(
                    #       tp_id, allele_id,
                    #       geno.object_properties['has_alternate_part'])
                    model.addComment(allele_id, tp_id)

                # derived_sf_assoc_alleles
                elif name == 'derived_sf_assoc_alleles':
                    # note the internal type_id changed
                    # from 129791 --> 133533 around 02.2016
                    # the relationship between
                    #   a reagent-feature and the allele it targets
                    # the relationship between
                    #   a reagent-targeted-gene (FBal) and
                    #   the reagent that targetes it (FBsf)

                    allele_id = None
                    reagent_id = None

                    if subject_id in self.idhash['allele']:
                        allele_id = self.idhash['allele'][subject_id]

                    if object_id in self.idhash['reagent']:
                        reagent_id = self.idhash['reagent'][object_id]

                    if allele_id is not None and reagent_id is not None:
                        g.addTriple(
                            allele_id,
                            geno.object_properties['targeted_by'], reagent_id)

                # produced by
                elif name == 'producedby':
                    # i'm looking just for the relationships between
                    # ti and tp features... so doing a bit of a hack
                    ti_id = None
                    tp_id = None
                    if subject_id in self.idhash['feature']:
                        ti_id = self.idhash['feature'][subject_id]
                        if not re.search(r'FBti', ti_id):
                            ti_id = None
                    if object_id in self.idhash['feature']:
                        tp_id = self.idhash['feature'][object_id]
                        if not re.search(r'FBtp', tp_id):
                            tp_id = None
                    if ti_id is not None and tp_id is not None:
                        geno.addSequenceDerivesFrom(ti_id, tp_id, )

                # gets_expression_data_from
                elif name == 'gets_expression_data_from':
                    # FIXME i don't know if this is correct
                    if subject_id in self.idhash['allele']:
                        allele_id = self.idhash['allele'][subject_id]
                    if object_id not in self.idhash['feature']:
                        continue
                    tp_id = self.idhash['feature'][object_id]
                    if not re.search(r'FBtp', tp_id):
                        tp_id = None
                        # TODO there are FBmc features here;
                        # need to incorporate if necessary
                        # these are markers i think
                    if tp_id is not None and allele_id is not None:
                        geno.addSequenceDerivesFrom(allele_id, tp_id)

                if not self.testMode \
                        and limit is not None and line_counter > limit:
                    break
        return

    # todo make singular
    def _process_organisms(self, limit):
        """
        The internal identifiers for the organisms in flybase

        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'organism'))
        logger.info("processing organisms")

        line_counter = 0
        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (organism_id, abbreviation, genus, species, common_name,
                 comment) = line
                # 1	Dmel	Drosophila	melanogaster	fruit fly
                # 2	Comp	Computational	result

                line_counter += 1

                tax_internal_id = self._makeInternalIdentifier(
                    'organism', organism_id)
                tax_label = ' '.join((genus, species))
                tax_id = tax_internal_id

                self.idhash['organism'][organism_id] = tax_id
                self.label_hash[tax_id] = tax_label

                # we won't actually add the organism to the graph,
                # unless we actually use it therefore it is added outside of
                # this function

                if self.testMode and\
                        int(organism_id) not in self.test_keys['organism']:
                    continue

                if not self.testMode and\
                        limit is not None and line_counter > limit:
                    pass
                else:
                    model.addClassToGraph(tax_id)
                    for s in [common_name, abbreviation]:
                        if s is not None and s.strip() != '':
                            model.addSynonym(tax_id, s)
                            model.addComment(tax_id, tax_internal_id)

        return

    def _get_organism_id(self, organism_key):
        organism_id = None
        tax_num = None
        if organism_key in self.idhash['organism']:
            organism_id = self.idhash['organism'][organism_key]

            # check to see if NCBITaxon has been resolved; if not fetch it
            if not re.match(r'NCBITaxon', organism_id):
                # NCBITaxon is not available in the dbxref or cvterm tables.
                # so we look them up using NCBI eutils services
                tax_label = self.label_hash[organism_id]
                # FIXME comment this out to speed things up
                tax_num = DipperUtil.get_ncbi_taxon_num_by_label(tax_label)
                if tax_num is not None:
                    organism_id = ':'.join(('NCBITaxon', tax_num))
                    self.idhash['organism'][organism_key] = organism_id
                    self.label_hash[organism_id] = tax_label

        return organism_id

    def _process_organism_dbxref(self, limit):
        """
        This is the mapping between the flybase organisms and
        external identifier "FBsp". We will want to use the NCBITaxon as
        the primary, if possible, but will default to a blank node/internal id
        if that is all that is available
        But we need to make the equivalences/sameAs.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'organism_dbxref'))
        logger.info("processing organsim dbxref mappings")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:

                (organism_dbxref_id, organism_id, dbxref_id, is_current) = line

                if self.testMode and\
                        int(organism_id) not in self.test_keys['organism']:
                    continue

                organism_key = organism_id
                if organism_key not in self.idhash['organism']:
                    continue
                organism_id = self.idhash['organism'][organism_key]

                dbxref_key = dbxref_id
                dbxrefs = self.dbxrefs.get(dbxref_key)
                if dbxrefs is not None:
                    for d in dbxrefs:
                        did = dbxrefs.get(d)
                        # don't make something sameAs itself
                        if did == organism_id:
                            continue
                        dlabel = self.label_hash.get(did)
                        model.addXref(organism_id, did)
                        if re.match(r'NCBITaxon', did):
                            model.makeLeader(did)
                        else:
                            model.addIndividualToGraph(did, dlabel)
                        line_counter += 1

                if not self.testMode and\
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_disease_models(self, limit):
        """
        Here we make associations between a disease and the supplied "model".
        In this case it's an allele.
        FIXME consider changing this... are alleles really models?
        Perhaps map these alleles into actual animals/strains or genotypes?
        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, self.files['disease_models']['file']))
        logger.info("processing disease models")

        line_counter = 0
        geno = Genotype(g)
        fly_taxon = 'NCBITaxon:7227'

        with gzip.open(raw, 'rb') as f:
            filereader = csv.reader(
                io.TextIOWrapper(f, newline=""),
                delimiter='\t', quotechar='\"')
            for line in filereader:
                # skip comments
                if re.match(r'#', ''.join(line)) or ''.join(line) == '':
                    continue
                (allele_id, allele_symbol, qualifier, doid_label, doid_id,
                 evidence_or_interacting_allele, pub_id) = line
                line_counter += 1

                if self.testMode and\
                        self.test_ids['disease'] is not None and\
                        doid_id not in self.test_ids['disease']:
                    continue

                rel = None
                allele_id = 'FlyBase:'+allele_id
                if qualifier == 'model of':
                    rel = model.object_properties['model_of']
                else:
                    # TODO amelorates, exacerbates, and DOES NOT *
                    continue

                animal_id = geno.make_experimental_model_with_genotype(
                    allele_id, allele_symbol, fly_taxon, 'fly')

                assoc = G2PAssoc(g, self.name, animal_id, doid_id, rel)
                if pub_id != '':
                    pub_id = 'FlyBase:'+pub_id
                    assoc.add_source(pub_id)
                if evidence_or_interacting_allele == \
                        'inferred from mutant phenotype':
                    evidence_id = 'ECO:0000015'
                    assoc.add_evidence(evidence_id)
                else:
                    assoc.set_description(evidence_or_interacting_allele)

                assoc.add_association_to_graph()

                if not self.testMode and\
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_stockprop(self, limit):
        """
        This will add depiction association between a strain and
        images hosted at flybase.
        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        raw = '/'.join((self.rawdir, 'stockprop'))
        logger.info("processing stock-image depictions")

        line_counter = 0

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                # skip comments
                if re.match(r'#', ''.join(line)) or ''.join(line) == '':
                    continue
                (stockprop_id, stock_id, type_id, value, rank) = line

                line_counter += 1

                if self.testMode \
                        and self.test_keys['strain'] is not None \
                        and int(stock_id) not in self.test_keys['strain']:
                    continue

                sid = self.idhash['stock'].get(stock_id)
                # linked_image
                if int(type_id) == 136340 and re.match(r'FBim', value):
                    # FIXME make sure this image url is perm
                    image_url = \
                        'http://flybase.org/tmp-shared/reports/'+value+'.png'
                    if sid is not None:
                        model.addDepiction(sid, image_url)
                    # TODO should this be a Reference object?

                # TODO add the stockprop_pub table when there is data to pull

                if not self.testMode and\
                        limit is not None and line_counter > limit:
                    break

        return

    def _get_human_models_file(self):
        """
        This function uses ftp to probe the FTP site to get the name of
        the current human_models file, and sets it in the files object.
        :return:

        """

        base_url = 'ftp.flybase.net'
        human_disease_dir = 'releases/current/precomputed_files/human_disease'
        from ftplib import FTP
        ftp = FTP(base_url)     # connect to host
        ftp.login()
        ftp.cwd(human_disease_dir)
        l = ftp.nlst()          # get list of files
        ftp.quit()
        f = None
        f_list = [
            i for i, x in enumerate(l)
            if re.match(r'allele_human_disease_model', x)]
        if len(f_list) == 0:
            logger.error("Can't find the human_disease_model file")
        elif len(f_list) > 1:
            logger.error(
                "There's >1 human disease model file, " +
                "and I don't know which to choose: %s", str(l))
        else:
            f = l[f_list[0]]

        if f is not None:
            # cat the url together
            file_url = '/'.join(('ftp:/', base_url, human_disease_dir, f))
            self.files['disease_models']['url'] = file_url

            # while we're at it, set the version...
            m = re.match(
                r'allele_human_disease_model_data_fb_(\d+_\d+).tsv.gz', f)
            # allele_human_disease_model_data_fb_2015_03.tsv.gz
            if m:
                ver = 'FB' + m.group(1)
                self.version_num = ver

        return

    def _makeInternalIdentifier(self, prefix, key):
        """
        This is a special Flybase-to-MONARCH-ism.
        Flybase tables have unique keys that we use here, but don't want to
        re-distribute those internal identifiers.
        Therefore, we make them into keys in a consistent way here.

        :param prefix: the object type to prefix the key with,
                    since the numbers themselves are not unique across tables
                    nor guarenteed stable within a table over time

        TEC: blank nodes created that way would not allways be valid RDF
        Better is to use the input as a rdfs:label & maybe rdf:type
        then digest to poduce a known valid key string.
        Also other rdf tools DO rewrite blank node identifiers

        :param key: the number (unique key)
        :return:

        """

        return '_:' + hashlib.sha1(
            ('fb'+prefix+'key'+key).encode('utf-8')).hexdigest()[1:20]

    def getTestSuite(self):
        import unittest
        from tests.test_flybase import FlyBaseTestCase

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(FlyBaseTestCase)

        return test_suite
