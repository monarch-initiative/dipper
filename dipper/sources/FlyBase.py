import logging
import re
import csv
import gzip
import io

from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.Environment import Environment
from dipper import config
from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils


logger = logging.getLogger(__name__)


class FlyBase(Source):
    """
    This is the [Drosophila Genetics](http://www.flybase.org/) resource,
    from which we process genotype and phenotype data about fruitfly.
    Genotypes leverage the GENO genotype model.

    Here, we connect to their public database, and download a subset of tables/views to get specifically at the
    geno-pheno data, then iterate over the tables.  We end up effectively performing joins when adding nodes
    to the graph.
    We connect using the [Direct Chado Access](http://gmod.org/wiki/Public_Chado_Databases#Direct_Chado_Access)

    When running the whole set, it performs best by dumping raw triples using the flag ```--format nt```.

    """

    tables = [
        'genotype',  # done
        'feature_genotype',
        'pub',  # done
        'feature_pub',  # done
        'pub_dbxref',
        'feature_dbxref',
        # 'feature_relationship',
        'cvterm',  # done
        'stock_genotype',  # done
        'stock',  # done
        'organism',
        # 'organism_dbxref',  # TODO
        'environment',   # done
        'phenotype',  # done
        'phenstatement',  # done
        'dbxref',  # done
        'phenotype_cvterm',  # done
        'phendesc',  # done
        # 'strain',  # may not be necessary
        'environment_cvterm',  # done
        # 'feature_cvterm'  # TODO to get better feature types than is in the feature table itself  (watch out for is_not)
    ]

    files = {
        'disease_models': {
            'file': 'allele_human_disease_model_data.tsv.gz',
            'url': 'ftp://ftp.flybase.net/releases/current/precomputed_files/human_disease/allele_human_disease_model_data_fb_*.tsv.gz'
        }
    }


    def __init__(self):
        Source.__init__(self, 'flybase')
        self.version_num = None   # to be used to store the version number to be acquired later

        self.namespaces.update(curie_map.get())

        # update the dataset object with details about this resource
        self.dataset = Dataset('flybase', 'FlyBase', 'http://www.flybase.org/', None,
                               None, 'http://flybase.org/wiki/FlyBase:FilesOverview')


        # source-specific warnings.  will be cleared when resolved.
        logger.warn("we are ignoring normal phenotypes for now")

        # so that we don't have to deal with BNodes, we will create hash lookups for the internal identifiers
        # the hash will hold the type-specific-object-keys to FB public identifiers.  then, subsequent
        # views of the table will lookup the identifiers in the hash.  this allows us to do the 'joining' on the
        # fly
        self.idhash = {'allele': {}, 'marker': {}, 'publication': {}, 'stock': {},
                       'genotype': {}, 'annot': {}, 'notes': {}, 'organism': {},
                       'environment': {}, 'feature': {}, 'phenotype': {}, 'cvterm': {}}
        self.dbxrefs = {}
        self.markers = {'classes': [], 'indiv': []}  # to store if a marker is a class or indiv

        self.label_hash = {}  # use this to store internally generated labels for various features
        self.geno_bkgd = {}  # use this to store the genotype strain ids for building genotype labels
        self.phenocv = {}  # mappings between internal phenotype db key and multiple cv terms

        return

    def fetch(self, is_dl_forced=False):
        """
        For the MGI resource, we connect to the remote database, and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        """

        # create the connection details for MGI
        cxn = {'host': 'flybase.org', 'database': 'flybase', 'port': 5432, 'user': 'flybase',
               'password': 'no password'}

        self.dataset.setFileAccessUrl(''.join(('jdbc:postgresql://', cxn['host'], ':',
                                               str(cxn['port']), '/', cxn['database'])))

        # process the tables
        # self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
#         self.fetch_from_pgdb(self.tables, cxn, None, is_dl_forced)

        # we want to fetch the features, but just a subset to reduce the processing time
        query = "select feature_id, dbxref_id, organism_id, name, uniquename, null as residues,"\
                +"seqlen, md5checksum, type_id, is_analysis, timeaccessioned, timelastmodified, is_obsolete "\
                +"from feature where is_analysis = false"

#        self.fetch_query_from_pgdb('feature', query, None, cxn, None, is_dl_forced)

        self._get_human_models_file()
        self.get_files(False)
        self.dataset.set_version_by_num(self.version_num)

        return

    def parse(self, limit=None):
        """
        We process each of the postgres tables in turn.  The order of processing is important here, as we build
        up a hashmap of internal vs external identifers (unique keys by type to FB id).  These include
        allele, marker (gene), publication, strain, genotype, annotation (association), and descriptive notes.
        :param limit: Only parse this many lines of each table
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %d rows of each file", limit)
        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self.nobnodes = True

        # the following will provide us the hash-lookups
        self._process_dbxref()
        self._process_cvterm()
        self._process_genotypes(limit)
        self._process_stocks(limit)
        self._process_pubs(limit)
        self._process_environment_cvterm()  # do this before environments to get the external ids
        self._process_environments()
        self._process_organisms(limit)
        self._process_features(limit)
        self._process_phenotype(limit)
        self._process_phenotype_cvterm()
        self._process_feature_dbxref(limit)  # gets external mappings for features (genes, variants, etc)

        # These are the associations amongst the objects above
        self._process_pub_dbxref(limit)
        self._process_phendesc(limit)
        self._process_feature_genotype(limit)
        self._process_feature_pub(limit)
        self._process_stock_genotype(limit)
        self._process_phenstatement(limit)   # these are G2P associations

        self._process_disease_models(limit)
        # TODO add version info from file somehow (in parser rather than during fetching)

        logger.info("Finished parsing.")

        self.load_bindings()
        for g in [self.graph, self.testgraph]:
            Assoc(self.name).load_all_properties(g)
            gu = GraphUtils(curie_map.get())
            gu.loadAllProperties(g)

        logger.info("Loaded %d nodes", len(self.graph))
        return

    def _process_genotypes(self, limit):
        """
        Add the genotype internal id to flybase mapping to the idhashmap.  Also, add them as individuals to the graph.

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

        line_counter = 0

        raw = '/'.join((self.rawdir, 'genotype'))
        logger.info("building labels for genotypes")
        gu = GraphUtils(curie_map.get())
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
                genotype_id = self._makeInternalIdentifier('genotype', genotype_num)

                self.idhash['genotype'][genotype_num] = genotype_id

                if description == '':
                    description = None

                if not self.testMode and limit is not None and line_counter > limit:
                    pass
                else:
                    gu.addIndividualToGraph(g, genotype_id, uniquename, Genotype.genoparts['intrinsic_genotype'],
                                            description)
                    if name.strip() != '':
                        gu.addSynonym(g, genotype_id, name)

        return

    def _process_stocks(self, limit):
        """
        Stock definitions.  Here we instantiate the Instances.

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0

        raw = '/'.join((self.rawdir, 'stock'))
        logger.info("building labels for stocks")
        gu = GraphUtils(curie_map.get())
        taxon = 'NCBITaxon:7227'
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (stock_id, dbxref_id, organism_id, name, uniquename, description,
                 type_id, is_obsolete) = line
# 2       12153979        1       2       FBst0000002     w[*]; betaTub60D[2] Kr[If-1]/CyO        10670
                # if self.testMode is True:
                #     if int(object_key) not in self.test_keys.get('genotype'):
                #         continue

                stock_num = stock_id
                stock_id = 'FlyBase:'+uniquename
                self.idhash['stock'][stock_num] = stock_id
                stock_label = description

                # todo look up the species by organism_id in the hashmap
                organism_key = organism_id
                # taxon = self.idhash['organsim'][organism_key]

                # todo do something with the dbxref_id (external ids?)

                if not self.testMode and limit is not None and line_counter > limit:
                    pass
                else:
                    if is_obsolete == 't':
                        gu.addDeprecatedIndividual(g, stock_id)
                    else:
                        gu.addIndividualToGraph(g, stock_id, stock_label, taxon)

        return

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

        line_counter = 0

        raw = '/'.join((self.rawdir, 'pub'))
        logger.info("building labels for pubs")
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            for line in filereader:
                (pub_id, title, volumetitle, volume, series_name, issue, pyear, pages, miniref,
                 type_id, is_obsolete, publisher, pubplace, uniquename) = line
# 2       12153979        1       2       FBst0000002     w[*]; betaTub60D[2] Kr[If-1]/CyO        10670
                # if self.testMode is True:
                #     if int(object_key) not in self.test_keys.get('genotype'):
                #         continue

                pub_num = pub_id
                pub_id = 'FlyBase:'+uniquename.strip()
                self.idhash['publication'][pub_num] = pub_id

                # TODO figure out the type of pub by type_id
                if not re.match('(FBrf|multi)', uniquename):
                    continue
                line_counter += 1

                r = Reference(pub_id)
                if title != '':
                    r.setTitle(title)
                if pyear != '':
                    r.setYear(str(pyear))
                if miniref != '':
                    r.setShortCitation(miniref)

                if not self.testMode and limit is not None and line_counter > limit:
                    pass
                else:
                    if is_obsolete == 't':
                        gu.addDeprecatedIndividual(g, pub_id)
                    else:
                        r.addRefToGraph(g)

        return

    def _process_environments(self):
        """
        There's only a few (~30) environments in which the phenotypes are recorded...
        There are no externally accessible identifiers for environments, so we make anonymous nodes for now.
        Some of the environments are comprised of >1 of the other environments; we do some simple parsing to
        match the strings of the environmental labels to the other atomic components.

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
                environment_internal_id = self._makeInternalIdentifier('environment', environment_num)
                if environment_num not in self.idhash['environment']:
                    self.idhash['environment'][environment_num] = environment_internal_id

                environment_id = self.idhash['environment'][environment_num]
                environment_label = uniquename
                env.addEnvironment(environment_id, environment_label)
                self.label_hash[environment_id] = environment_label

                # split up the environment into parts
                # if there's parts, then add them to the hash; we'll match the components in a second pass
                components = re.split('\|', uniquename)
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

    def _process_features(self, limit):
        """
        These are all of the genomic features (genes, variations, transgenes, etc.).
        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        raw = '/'.join((self.rawdir, 'feature'))
        logger.info("building labels for features")

        line_counter = 0
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (feature_id, dbxref_id, organism_id, name, uniquename, residues, seqlen, md5checksum, type_id,
                 is_analysis, timeaccessioned, timelastmodified, is_obsolete) = line

                line_counter += 1

                feature_key = feature_id
                if re.search('[\|\s\[\]\{\}\\\\<\>]', uniquename):
                    # some uniquenames have pipes or other nasty chars!  for example: FB||||FBrf0133242|Hugh-u1
                    feature_id = self._makeInternalIdentifier('feature', feature_key)
                else:
                    feature_id = 'FlyBase:'+uniquename
                self.idhash['feature'][feature_key] = feature_id

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
                    ]

                # FIXME should we skip introns/exons?
                type_keys_to_skip = [
                    596,  # pcr_product
                    57096,
                    57097,
                    57270,
                    58210,
                    59643,
                    60006,
                    61351,
                    61467,
                    257,  # exon
                    286,  # intron
                ]

                if type_id in types_to_skip or int(type_key) in type_keys_to_skip:
                    continue

                tax_id = self._makeInternalIdentifier('organism', organism_id)

                # HACK - FBgn are genes, and therefore classes, all else be individuals
                is_gene = False
                if re.search('(FBgn|FBog)', feature_id):
                    is_gene = True

                if not self.testMode and limit is not None and line_counter > limit:
                    pass
                else:
                    if is_gene:
                        gu.addClassToGraph(g, feature_id, name, type_id)
                    else:
                        gu.addIndividualToGraph(g, feature_id, name, type_id)

                    if is_obsolete == 't':
                        if is_gene:
                            gu.addDeprecatedClass(g, feature_id)
                        else:
                            gu.addDeprecatedIndividual(g, feature_id)

                    gu.addTriple(g, feature_id, gu.object_properties['in_taxon'], tax_id)

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
                (feature_genotype_id, feature_id, genotype_id, chromosome_id, rank, cgroup, cvterm_id) = line
                # 1	23273518	2	23159230	0	0	60468
                feature_key = feature_id
                feature_id = self.idhash['feature'][feature_key]
                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]

                # what is cvterm_id for in this context???
                # cgroup is the order of composition of things in the genotype label.
                # sometimes the same feature is listed twice; not sure if this is a mistake, or zygosity, or?

                geno.addParts(feature_id, genotype_id, geno.object_properties['has_alternate_part'])

                # TODO we will build up the genotypes here... lots to do

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_phendesc(self, limit):
        """
        The description of the resulting phenotypes with the genotype+environment

        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        raw = '/'.join((self.rawdir, 'phendesc'))
        logger.info("processing G2P")

        line_counter = 0
        gu = GraphUtils(curie_map.get())
        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (phendesc_id, genotype_id, environment_id, description, type_id, pub_id) = line
                # 1	2	1	Hemizygous males are wild type, homozygous males are sterile.	60466	209729

                line_counter += 1
                phendesc_key = phendesc_id
                phendesc_id = self._makeInternalIdentifier('phendesc', phendesc_key)

                # for now, just attach the description to the genotype
                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]
                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                environment_key = environment_id
                environment_id = self.idhash['environment'][environment_key]

                # TODO type id ==> ECO???

                # just make associations with abnormal phenotype
                phenotype_id = 'FBcv:0001347'
                assoc = G2PAssoc(self.name, genotype_id, phenotype_id)
                assoc.add_source(pub_id)
                assoc.set_description(description)
                assoc.set_environment(environment_id)
                assoc.add_association_to_graph(g)
                assoc_id = assoc.get_association_id()
                gu.addComment(g, assoc_id, phendesc_id)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_feature_pub(self, limit):
        """
        The description of the resulting phenotypes with the genotype+environment

        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        raw = '/'.join((self.rawdir, 'feature_pub'))
        logger.info("processing feature_pub")

        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (feature_pub_id, feature_id, pub_id) = line
                # 1440    3175682 62137
                # 2       3160606 99159
                feature_key = feature_id
                if feature_key not in self.idhash['feature']:
                    continue
                feature_id = self.idhash['feature'][feature_key]
                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                gu.addTriple(g, pub_id, gu.object_properties['mentions'], feature_id)

                line_counter += 1

                if not self.testMode and limit is not None and line_counter > limit:
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

        raw = '/'.join((self.rawdir, 'stock_genotype'))
        logger.info("processing stock genotype")
        geno = Genotype(g)
        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (stock_genotype_id, stock_id, genotype_id) = line

                stock_key = stock_id
                stock_id = self.idhash['stock'][stock_key]
                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]

                gu.addTriple(g, stock_id, geno.object_properties['has_genotype'], genotype_id)

                line_counter += 1

                if not self.testMode and limit is not None and line_counter > limit:
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

        raw = '/'.join((self.rawdir, 'pub_dbxref'))
        logger.info("processing pub_dbxref")

        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (pub_dbxref_id, pub_id, dbxref_id, is_current) = line
                # 49648	43222	395730	t

                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                # get any dbxrefs for pubs, including pmids and dois
                dbxref_key = dbxref_id
                if str(dbxref_key) in self.dbxrefs:
                    dbxrefs = self.dbxrefs[str(dbxref_key)]
                    # pub_dbs = [75, 51, 76, 95, 126]
                    pmid_ids = [50, 77, 275, 286, 347]
                    flybase_ids = [4]
                    isbn = [75, 51]
                    for d in dbxrefs:
                        dbxref_id = None
                        if int(d) in pmid_ids:
                            dbxref_id = 'PMID:'+dbxrefs[d].strip()
                        elif int(d) in isbn:
                            dbxref_id = 'ISBN:'+dbxrefs[d].strip()
                        elif int(d) == 161:
                            dbxref_id = 'DOI:'+dbxrefs[d].strip()
                        # elif int(d) == 4:
                        #     dbxref_id = 'FlyBase:'+dbxrefs[d].strip()

                        if dbxref_id is not None:
                            r = Reference(dbxref_id, Reference.ref_types['publication'])
                            r.addRefToGraph(g)
                            gu.addSameIndividual(g, pub_id, dbxref_id)
                            line_counter += 1

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_dbxref(self):
        """
        We bring in the dbxref identifiers and store them in a hashmap for lookup in other functions

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

                db_ids = {50: 'PMID',  # pubmed
                          68: 'RO',  # obo-rel
                          71: 'FBdv',  # FBdv
                          74: 'FBbt',  # FBbt
                          # 28:,  # genbank
                          30: 'OMIM',  # MIM
                          # 38,  # ncbi
                          75: 'ISBN',  # ISBN
                          46: 'PMID',  # PUBMED
                          51: 'ISBN',  # isbn
                          52: 'SO',  # so
                          # 76,  # http
                          77: 'PMID',  # PMID
                          80: 'FBcv',  # FBcv
                          # 95,  # MEDLINE
                          98: 'REACT',  # Reactome
                          103: 'CHEBI', # Chebi
                          102: 'MESH', # MeSH
                          106: 'OMIM', # OMIM
                          105: 'KEGG-path', # KEGG pathway
                          107: 'DOI', # doi
                          108: 'CL', # CL
                          114: 'CHEBI', # CHEBI
                          115: 'KEGG', # KEGG
                          116: 'PubChem', # PubChem
                          # 120, # MA???
                          3: 'GO',   # GO
                          4: 'FlyBase',   # FlyBase
                          # 126, # URL
                          128: 'PATO', # PATO
                          # 131, # IMG
                          2: 'SO',   # SO
                          136: 'MESH', # MESH
                          139: 'CARO', # CARO
                          140: 'NCBITaxon', # NCBITaxon
                          # 151, # MP  ???
                          161: 'DOI', # doi
                          36: 'BDGP',  # BDGP
                          # 55,  # DGRC
                          # 54,  # DRSC
                          # 169, # Transgenic RNAi project???
                          231: 'RO', # RO ???
                          180: 'NCBIGene', # entrezgene
                          # 192, # Bloomington stock center
                          197: 'UBERON', # Uberon
                          212: 'ENSEMBL', # Ensembl
                          # 129, # GenomeRNAi
                          275: 'PMID', # PubMed
                          286: 'PMID', # pmid
                          265: 'OMIM', # OMIM_Gene
                          266: 'OMIM', # OMIM_Phenotype
                          300: 'DOID', # DOID
                          302: 'MESH', # MSH
                          347: 'PMID', # Pubmed
                }  # the databases to fetch

                if accession.strip() != '' and int(db_id) in db_ids:
                    # scrub some identifiers here
                    m = re.match('(doi|SO|GO|FBcv|FBbt_root|FBdv|FBgn|FBdv_root|FlyBase|FBbt):', accession)
                    if m:
                        accession = re.sub(m.group(1)+'\:', '', accession)
                    elif re.match('(FlyBase miscellaneous CV|cell_lineprop|relationship type|FBgn$)', accession):
                        continue
                    elif re.match('\:', accession):  # starts with a colon
                        accession = re.sub('\:', '', accession)
                    elif re.search('\s', accession):
                        # skip anything with a space
                        # logger.debug('dbxref %s accession has a space: %s', dbxref_id, accession)
                        continue

                    if re.match('http', accession):
                        did = accession.strip()
                    else:
                        prefix = db_ids.get(int(db_id))
                        did = ':'.join((prefix, accession.strip()))
                        if re.search('\:', accession) and prefix != 'DOI':
                            logger.warn('id %s may be malformed; skipping', did)

                    self.dbxrefs[dbxref_id] = {db_id: did}

                elif url != '':
                    self.dbxrefs[dbxref_id] = {db_id: url.strip()}
                else:
                    continue

                line_counter += 1

        return

    def _process_phenotype(self, limit):
        """
        Get the phenotypes, and declare the classes.
        If the "observable" is "unspecified", then we assign the phenotype to the "cvalue" id; otherwise
        we convert the phenotype into a uberpheno-style identifier, simply based on
        the anatomical part that's affected...that is listed as the observable_id,
        concatenated with the literal "PHENOTYPE"

        Note that assay_id is the same for all current items, so we do nothing with this.
        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        raw = '/'.join((self.rawdir, 'phenotype'))
        logger.info("processing phenotype")

        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (phenotype_id, uniquename, observable_id, attr_id, value, cvalue_id, assay_id) = line

                # 8505	unspecified
                # 20142	mesothoracic leg disc | somatic clone	87719	60468		60468	60468
                # 8507	sex comb | ectopic	88877	60468		60468	60468
                # 8508	tarsal segment	83664	60468		60468	60468

                # for now make these as phenotypic classes - will need to xref at some point
                phenotype_key = phenotype_id
                phenotype_id = None
                phenotype_internal_id = self._makeInternalIdentifier('phenotype', phenotype_key)
                phenotype_label = None
                self.label_hash[phenotype_internal_id] = uniquename
                cvterm_id = None
                if observable_id != '' and int(observable_id) == 60468:  # undefined - typically these are already phenotypes
                    if cvalue_id in self.idhash['cvterm']:
                        cvterm_id = self.idhash['cvterm'][cvalue_id]
                        phenotype_id = self.idhash['cvterm'][cvalue_id]
                elif observable_id in self.idhash['cvterm']:  # observations to anatomical classes
                    cvterm_id = self.idhash['cvterm'][observable_id]
                    phenotype_id = self.idhash['cvterm'][observable_id] + 'PHENOTYPE'
                    if cvterm_id is not None and cvterm_id in self.label_hash:
                        phenotype_label = self.label_hash[cvterm_id]
                        phenotype_label += ' phenotype'
                        self.label_hash[phenotype_id] = phenotype_label
                    else:
                        logger.info('cvtermid=%s not in label_hash', cvterm_id)

                else:
                    logger.info("No observable id or label for %s: %s", phenotype_key, uniquename)


                # TODO store this composite phenotype in some way as a proper class definition?
                self.idhash['phenotype'][phenotype_key] = phenotype_id

                # assay_id is currently only "undefined" key=60468

                if not self.testMode and limit is not None and line_counter > limit:
                    pass
                else:
                    if phenotype_id is not None:
                        # assume that these fit into the phenotypic uberpheno elsewhere
                        gu.addClassToGraph(g, phenotype_id, phenotype_label)
                        line_counter += 1

        return

    def _process_phenstatement(self, limit):
        """
        The phenstatements are the genotype-to-phenotype associations, in the context of an environment.
        These are also curated to a publication.
        So we make oban associations, adding the pubs as a source.  We additionally add the internal key
        as a comment for tracking purposes.
        :param limit:
        :return:
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        raw = '/'.join((self.rawdir, 'phenstatement'))
        logger.info("processing phenstatement")

        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (phenstatement_id, genotype_id, environment_id, phenotype_id, type_id, pub_id) = line

                # 168549	166695	1	8507	60468	151256
                # 168550	166695	1	8508	60468	151256
                # 168551	166696	1	8509	60468	151256
                # 168552	166696	1	8510	60468	151256
                line_counter += 1
                phenstatement_key = phenstatement_id
                phenstatement_id = self._makeInternalIdentifier('phenstatement', phenstatement_key)
                genotype_key = genotype_id
                genotype_id = self.idhash['genotype'][genotype_key]
                environment_key = environment_id
                environment_id = self.idhash['environment'][environment_key]
                phenotype_key = phenotype_id
                phenotype_internal_id = self._makeInternalIdentifier('phenotype', phenotype_key)  # TEMP
                phenotype_internal_label = self.label_hash[phenotype_internal_id]
                phenotype_id = self.idhash['phenotype'][phenotype_key]
                pub_key = pub_id
                pub_id = self.idhash['publication'][pub_key]

                assoc = G2PAssoc(self.name, genotype_id, phenotype_id)
                assoc.set_environment(environment_id)
                assoc.add_source(pub_id)
                assoc.add_association_to_graph(g)
                assoc_id = assoc.get_association_id()
                gu.addComment(g, assoc_id, phenstatement_id)
                gu.addDescription(g, assoc_id, phenotype_internal_label)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    def _process_phenotype_cvterm(self):
        """
        These are the qualifiers for the phenotype location itself.
        We don't really do anything with these yet.
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

                    self.idhash['phenotype'][phenotype_key] = cvterm_id
                    if phenotype_key not in self.phenocv:
                        self.phenocv[phenotype_id] = [cvterm_id]
                    else:
                        self.phenocv[phenotype_id] += [cvterm_id]
                else:
                    logger.info("Not storing the cvterm info for %s", cvterm_key)

        return

    def _process_cvterm(self):
        """
        CVterms are the internal identifiers for any controlled vocab or ontology term.  Many are
        xrefd to actual ontologies.  The actual external id is stored in the dbxref table, which we
        place into the internal hashmap for lookup with the cvterm id.  The name of the external term
        is stored in the "name" element of this table, and we add that to the label hashmap for lookup elsewhere

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

                (cvterm_id, cv_id, definition, dbxref_id, is_obsolete, is_relationshiptype, name) = line

                # 316	6		1665919	0	0	rRNA_cleavage_snoRNA_primary_transcript
                # 28	5		1663309	0	0	synonym
                # 455	6		1665920	0	0	tmRNA

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
                # look up the dbxref_id for the cvterm - hopefully it's one-to-one
                dbxrefs = self.dbxrefs.get(dbxref_id)
                if dbxrefs is not None:
                    if len(dbxrefs) > 1:
                        logger.info(">1 dbxref for this cvterm (%s: %s): ", str(cvterm_id), name, dbxrefs.values())
                    elif len(dbxrefs) == 1:
                        # replace the cvterm with the dbxref (external) identifier
                        did = dbxrefs.popitem()[1]
                        self.idhash['cvterm'][cvterm_key] = did  # get the value
                        # also add the label to the dbxref
                        self.label_hash[did] = name
        return

    def _process_environment_cvterm(self):
        """
        This is the mapping between the internal environment id and the external ones; here we
        map the internal environment id to the external one in the hashmap.
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

                # look up the dbxref_id for the cvterm - hopefully it's one-to-one
                self.idhash['environment'][environment_key] = cvterm_id

        return

    def _process_feature_dbxref(self, limit):
        """
        This is the mapping between the flybase features and external repositories.
        Generally we want to leave the flybase feature id as the primary identifier.
        But we need to make the equivalences/sameAs.

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        line_counter = 0
        raw = '/'.join((self.rawdir, 'feature_dbxref'))
        logger.info("processing feature dbxref mappings")
        gu = GraphUtils(curie_map.get())
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
                if feature_key not in self.idhash['feature']:
                    # some features may not be found in the hash if they are "analysis features"
                    # logger.debug("Feature %s not found in hash", feature_key)
                    continue
                feature_id = self.idhash['feature'][feature_key]

                dbxref_key = dbxref_id
                dbxrefs = self.dbxrefs.get(dbxref_key)
                if dbxrefs is not None:
                    for d in dbxrefs:
                        # need to filter based on db ?
                        # TODO should we make other species' identifiers primary, instead of flybase?
                        did = dbxrefs.get(d)
                        if did == feature_id:  # don't make something sameAs itself
                            continue
                        dlabel = self.label_hash.get(did)
                        gu.addIndividualToGraph(g, did, dlabel)
                        gu.addSameIndividual(g, feature_id, did)
                        line_counter += 1

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

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

        raw = '/'.join((self.rawdir, 'organism'))
        logger.info("processing organisms")

        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            filereader = csv.reader(f, delimiter='\t', quotechar='\"')
            f.readline()  # read the header row; skip
            for line in filereader:
                (organism_id, abbreviation, genus, species, common_name, comment) = line
                # 1	Dmel	Drosophila	melanogaster	fruit fly
                # 2	Comp	Computational	result

                line_counter += 1

                tax_id = self._makeInternalIdentifier('organism', organism_id)
                tax_label = ' '.join((genus, species))

                self.idhash['organism'][organism_id] = tax_id
                self.label_hash[tax_id] = tax_label

                if not self.testMode and limit is not None and line_counter > limit:
                    pass
                else:
                    gu.addClassToGraph(g, tax_id, tax_label)
                    for s in [common_name, abbreviation]:
                        if s is not None and s.strip() != '':
                            gu.addSynonym(g, tax_id, s)

        return

    def _process_disease_models(self, limit):

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        raw = '/'.join((self.rawdir, self.files['disease_models']['file']))
        logger.info("processing disease models")

        line_counter = 0
        gu = GraphUtils(curie_map.get())

        with gzip.open(raw, 'rb') as f:
            filereader = csv.reader(io.TextIOWrapper(f, newline=""), delimiter='\t', quotechar='\"')
            for line in filereader:
                if re.match('#', ''.join(line)) or ''.join(line) == '':  # skip comments
                    continue
                (allele_id, allele_symbol, qualifier, doid_label, doid_id, evidence_or_interacting_allele, pub_id) = line
                line_counter += 1
                rel = None
                allele_id = 'FlyBase:'+allele_id
                if qualifier == 'model of':
                    rel = gu.object_properties['model_of']
                else:
                    # TODO amelorates, exacerbates, and DOES NOT *
                    continue
                assoc = G2PAssoc(self.name, allele_id, doid_id, rel)
                if pub_id != '':
                    pub_id = 'FlyBase:'+pub_id
                    assoc.add_source(pub_id)
                if evidence_or_interacting_allele == 'inferred from mutant phenotype':
                    evidence_id = 'ECO:0000015'
                    assoc.add_evidence(evidence_id)
                else:
                    assoc.set_description(evidence_or_interacting_allele)

                assoc.add_association_to_graph(g)

        return

    def _get_human_models_file(self):
        """
        This function uses ftp to probe the FTP site to get the name of the current human_models file,
        and sets it in the files object.
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
        f_list = [i for i, x in enumerate(l) if re.match('allele_human_disease_model', x)]
        if len(f_list) == 0:
            logger.error("Can't find the human_disease_model file")
        elif len(f_list) > 1:
            logger.error("There's >1 human disease model file, and I don't know which to choose: %s", str(l))
        else:
            f = l[f_list[0]]

        if f is not None:
            # cat the url together
            file_url = '/'.join(('ftp:/', base_url, human_disease_dir, f))
            self.files['disease_models']['url'] = file_url

            # while we're at it, set the version...
            m = re.match('allele_human_disease_model_data_fb_(\d+_\d+).tsv.gz', f)
            # allele_human_disease_model_data_fb_2015_03.tsv.gz
            if m:
                ver = 'FB' + m.group(1)
                self.version_num = ver

        return


    def _makeInternalIdentifier(self, prefix, key):
        """
        This is a special Flybase-to-MONARCH-ism.  Flybase tables have unique keys that we use here, but don't want
        to necessarily re-distribute those internal identifiers.  Therefore, we make them into keys in a consistent
        way here.
        :param prefix: the object type to prefix the key with, since the numbers themselves are not unique across tables
        :param key: the number (unique key)
        :return:
        """
        iid = '_fb'+prefix+'key'+key
        if self.nobnodes:
            iid = ':'+ iid

        return iid

    # def getTestSuite(self):
    #     import unittest
    #     from tests.test_flybase import FlyBaseTestCase
    #
    #     test_suite = unittest.TestLoader().loadTestsFromTestCase(FlyBaseTestCase)
    #
    #     return test_suite