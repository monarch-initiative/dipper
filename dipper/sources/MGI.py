import csv
import os
from datetime import datetime
import logging
import re

from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.assoc.Association import Assoc
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper import config
from dipper.models.GenomicFeature import Feature, makeChromID


logger = logging.getLogger(__name__)


class MGI(PostgreSQLSource):
    """
    This is the
    [Mouse Genome Informatics](http://www.informatics.jax.org/) resource,
    from which we process genotype and phenotype data about laboratory mice.
    Genotypes leverage the GENO genotype model.

    Here, we connect to their public database, and download a subset of
    tables/views to get specifically at the geno-pheno data,
    then iterate over the tables.  We end up effectively performing joins
    when adding nodes to the graph.
    In order to use this parser, you will need to have user/password connection
    details in your conf.json file, like:
    dbauth : {'mgi' : {'user' : '<username>', 'password' : '<password>'}}
    You can request access by contacting mgi-help@jax.org

    """

    # CONSIDER IF WE NEED:
    # mgi_organism_acc_view:
    #    Consider using this for the taxon mapping instead of
    #   the hashmap encoded below
    # mgi_reference_allele_view:
    #   Don't believe this view is used in either
    #    the genotype of phenotype view
    # all_allele_cellline_view: When we want to start dealing with cell lines
    # mgi_note_strain_view: prose descriptions of strains.
    # prb_strain_summary_view:
    #   Don't believe this view is used in
    #    either the genotype of phenotype view
    # prb_strain_marker_view:
    #    eventually i think we want this because
    # it has other relevant markers that are affected

    resources = [
        {
          'query': '../../resources/sql/mgi/mgi_dbinfo.sql',
          'outfile': 'mgi_dbinfo',
          'Force': True
        },
        {
          'query': '../../resources/sql/mgi/gxd_genotype_view.sql',
          'outfile': 'gxd_genotype_view'
        },
        {
          'query': '../../resources/sql/mgi/gxd_genotype_summary_view.sql',
          'outfile': 'gxd_genotype_summary_view'
        },
        {
          'query': '../../resources/sql/mgi/gxd_allelepair_view.sql',
          'outfile': 'gxd_allelepair_view'
        },
        {
          'query': '../../resources/sql/mgi/all_summary_view.sql',
          'outfile': 'all_summary_view'
        },
        {
          'query': '../../resources/sql/mgi/all_allele_view.sql',
          'outfile': 'all_allele_view'
        },
        {
          'query': '../../resources/sql/mgi/all_allele_mutation_view.sql',
          'outfile': 'all_allele_mutation_view'
        },
        {
          'query': '../../resources/sql/mgi/mrk_marker_view.sql',
          'outfile': 'mrk_marker_view'
        },
        {
          'query': '../../resources/sql/mgi/voc_annot_view.sql',
          'outfile': 'voc_annot_view'
        },
        {
          'query': '../../resources/sql/mgi/voc_evidence_view.sql',
          'outfile': 'voc_evidence_view'
        },
        {
          'query': '../../resources/sql/mgi/bib_acc_view.sql',
          'outfile': 'bib_acc_view'
        },
        {
           'query': '../../resources/sql/mgi/prb_strain_view.sql',
           'outfile': 'prb_strain_view'
        },
        {
          'query': '../../resources/sql/mgi/mrk_summary_view.sql',
          'outfile': 'mrk_summary_view'
        },
        {
          'query': '../../resources/sql/mgi/mrk_acc_view.sql',
          'outfile': 'mrk_acc_view'
        },
        {
          'query': '../../resources/sql/mgi/prb_strain_acc_view.sql',
          'outfile': 'prb_strain_acc_view'
        },
        {
          'query': '../../resources/sql/mgi/prb_strain_genotype_view.sql',
          'outfile': 'prb_strain_genotype_view'
        },
        {
          'query': '../../resources/sql/mgi/mgi_note_vocevidence_view.sql',
          'outfile': 'mgi_note_vocevidence_view'
        },
        {
          'query': '../../resources/sql/mgi/mgi_note_allele_view.sql',
          'outfile': 'mgi_note_allele_view'
        },
        {
          'query': '../../resources/sql/mgi/mrk_location_cache.sql',
          'outfile': 'mrk_location_cache'  # gene locations
        }
    ]

    # for testing purposes, this is a list of internal db keys
    # to match and select only portions of the source
    test_keys = {
        'allele': [
            1612, 1609, 1303, 56760, 816699, 51074, 14595, 816707, 246, 38139,
            4334, 817387, 8567, 476, 42885, 3658, 1193, 6978, 6598, 16698,
            626329, 33649, 835532, 7861, 33649, 6308, 1285, 827608],
        'marker': [
            357, 38043, 305574, 444020, 34578, 9503, 38712, 17679, 445717,
            38415, 12944, 377, 77197, 18436, 30157, 14252, 412465, 38598,
            185833, 35408, 118781, 37270, 31169, 25040, 81079],
        'annot': [
            6778, 12035, 189442, 189443, 189444, 189445, 189446, 189447,
            189448, 189449, 189450, 189451, 189452, 318424, 717023, 717024,
            717025, 717026, 717027, 717028, 717029, 5123647, 928426, 5647502,
            6173775, 6173778, 6173780, 6173781, 6620086, 13487622, 13487623,
            13487624, 23241933, 23534428, 23535949, 23546035, 24722398,
            29645663, 29645664, 29645665, 29645666, 29645667, 29645682,
            43803707, 43804057, 43805682, 43815003, 43838073, 58485679,
            59357863, 59357864, 59357865, 59357866, 59357867, 60448185,
            60448186, 60448187, 62628962, 69611011, 69611253, 79642481,
            79655585, 80436328, 83942519, 84201418, 90942381, 90942382,
            90942384, 90942385, 90942386, 90942389, 90942390, 90942391,
            90942392, 92947717, 92947729, 92947735, 92947757, 92948169,
            92948441, 92948518, 92949200, 92949301, 93092368, 93092369,
            93092370, 93092371, 93092372, 93092373, 93092374, 93092375,
            93092376, 93092377, 93092378, 93092379, 93092380, 93092381,
            93092382, 93401080, 93419639, 93436973, 93436974, 93436975,
            93436976, 93436977, 93459094, 93459095, 93459096, 93459097,
            93484431, 93484432, 93491333, 93491334, 93491335, 93491336,
            93491337, 93510296, 93510297, 93510298, 93510299, 93510300,
            93548463, 93551440, 93552054, 93576058, 93579091, 93579870,
            93581813, 93581832, 93581841, 93581890, 93583073, 93583786,
            93584586, 93587213, 93604448, 93607816, 93613038, 93614265,
            93618579, 93620355, 93621390, 93624755, 93626409, 93626918,
            93636629, 93642680, 93643814, 93643825, 93647695, 93648755,
            93652704, 5123647, 71668107, 71668108, 71668109, 71668110,
            71668111, 71668112, 71668113, 71668114, 74136778, 107386012,
            58485691],
        'genotype': [
            81, 87, 142, 206, 281, 283, 286, 287, 341, 350, 384, 406, 407, 411,
            425, 457, 458, 461, 476, 485, 537, 546, 551, 553, 11702, 12910,
            13407, 13453, 14815, 26655, 28610, 37313, 38345, 59766, 60082,
            65406, 64235],
        'pub': [
            73197, 165659, 134151, 76922, 181903, 26681, 128938, 80054, 156949,
            159965, 53672, 170462, 206876, 87798, 100777, 176693, 139205,
            73199, 74017, 102010, 152095, 18062, 216614, 61933, 13385, 32366,
            114625, 182408, 140802],
        'strain': [
            30639, 33832, 33875, 33940, 36012, 59504, 34338, 34382, 47670,
            59802, 33946, 31421, 64, 40, 14, -2, 30639, 15975, 35077, 12610,
            -1, 28319, 27026, 141, 62299],
        'notes': [
            5114, 53310, 53311, 53312, 53313, 53314, 53315, 53316, 53317,
            53318, 53319, 53320, 71099, 501751, 501752, 501753, 501754, 501755,
            501756, 501757, 744108, 1055341, 6049949, 6621213, 6621216,
            6621218, 6621219, 7108498, 14590363, 14590364, 14590365, 25123358,
            25123360, 26688159, 32028545, 32028546, 32028547, 32028548,
            32028549, 32028564, 37833486, 47742903, 47743253, 47744878,
            47754199, 47777269, 65105483, 66144014, 66144015, 66144016,
            66144017, 66144018, 70046116, 78382808, 78383050, 103920312,
            103920318, 103920319, 103920320, 103920322, 103920323, 103920324,
            103920325, 103920326, 103920328, 103920330, 103920331, 103920332,
            103920333, 106390006, 106390018, 106390024, 106390046, 106390458,
            106390730, 106390807, 106391489, 106391590, 106579450, 106579451,
            106579452, 106579453, 106579454, 106579455, 106579456, 106579457,
            106579458, 106579459, 106579460, 106579461, 106579462, 106579463,
            106579464, 106949909, 106949910, 106969368, 106969369, 106996040,
            106996041, 106996042, 106996043, 106996044, 107022123, 107022124,
            107022125, 107022126, 107052057, 107052058, 107058959, 107058960,
            107058961, 107058962, 107058963, 107077922, 107077923, 107077924,
            107077925, 107077926, 107116089, 107119066, 107119680, 107154485,
            107155254, 107158128, 107159385, 107160435, 107163154, 107163183,
            107163196, 107163271, 107164877, 107165872, 107166942, 107168838,
            107170557, 107174867, 107194346, 107198590, 107205179, 107206725,
            107212120, 107214364, 107214911, 107215700, 107218519, 107218642,
            107219974, 107221415, 107222064, 107222717, 107235068, 107237686,
            107242709, 107244121, 107244139, 107248964, 107249091, 107250401,
            107251870, 107255383, 107256603]
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'mgi')

        # update the dataset object with details about this resource
        self.dataset = Dataset(
            'mgi', 'MGI', 'http://www.informatics.jax.org/', None,
            'http://www.informatics.jax.org/mgihome/other/copyright.shtml')

        # check if config exists; if it doesn't, error out and let user know
        if 'dbauth' not in config.get_config() and \
                'mgi' not in config.get_config()['dbauth']:
            logger.error("not configured with PG user/password.")

        # source-specific warnings.  will be cleared when resolved.
        logger.warning("we are ignoring normal phenotypes for now")

        # so that we don't have to deal with BNodes,
        # we will create hash lookups
        # for the internal identifiers the hash will hold
        # the type-specific-object-keys to MGI public identifiers.
        # then, subsequent views of the table will lookup the identifiers
        # in the hash.  this allows us to do the 'joining' on the fly
        self.idhash = {
            'allele': {}, 'marker': {}, 'publication': {}, 'strain': {},
            'genotype': {}, 'annot': {}, 'notes': {}, 'seqalt': {}}
        # to store if a marker is a class or indiv
        self.markers = {
            'classes': [], 'indiv': []}
        # use this to store internally generated labels for various features
        self.label_hash = {}
        # use this to store the genotype strain ids
        # for building genotype labels
        self.geno_bkgd = {}
        self.strain_to_genotype_map = {}

        self.wildtype_alleles = set()

        # also add the gene ids from the config
        # in order to capture transgenes of the test set
        if 'test_ids' not in config.get_config() or\
                'gene' not in config.get_config()['test_ids']:
            logger.warning("not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['gene']
        return

    def fetch(self, is_dl_forced=False):
        """
        For the MGI resource, we connect to the remote database,
        and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        """

        # create the connection details for MGI
        cxn = config.get_config()['dbauth']['mgi']
        cxn.update(
            {'host': 'mgi-adhoc.jax.org', 'database': 'mgd', 'port': 5432})

        self.dataset.setFileAccessUrl(
            ''.join(('jdbc:postgresql://', cxn['host'], ':', str(cxn['port']),
                    '/', cxn['database'])), is_object_literal=True)

        # process the tables
        # self.fetch_from_pgdb(self.tables, cxn, 100)  # for testing only
        # self.fetch_from_pgdb(self.tables, cxn, None, is_dl_forced)

        for query_map in self.resources:
            query_fh = open(os.path.join(
                os.path.dirname(__file__), query_map['query']), 'r')
            query = query_fh.read()
            force = False
            if 'Force' in query_map:
                force = query_map['Force']
            self.fetch_query_from_pgdb(
                query_map['outfile'], query, None, cxn, force=force)
        # always get this - it has the verion info
        self.fetch_transgene_genes_from_db(cxn)

        datestamp = ver = None
        # get the resource version information from
        # table mgi_dbinfo, already fetched above
        outfile = '/'.join((self.rawdir, 'mgi_dbinfo'))

        if os.path.exists(outfile):
            with open(outfile, 'r') as f:
                f.readline()  # read the header row; skip
                info = f.readline()
                cols = info.split('\t')
                ver = cols[0]  # col 0 is public_version
                ver = ver.replace('MGI ', '')  # MGI 5.20 --> 5.20
                # MGI has a datestamp for the data within the database;
                # use it instead of the download date
                # datestamp in the table: 2014-12-23 00:14:20[.12345]
                # modification date without micro seconds
                (d, ms) = cols[1].strip().split('.')
                datestamp = \
                    datetime.strptime(
                        d, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
                f.close()

        self.dataset.setVersion(datestamp, ver)

        return

    def parse(self, limit=None):
        """
        We process each of the postgres tables in turn.
        The order of processing is important here, as we build
        up a hashmap of internal vs external identifers
        (unique keys by type to MGI id).  These include allele, marker (gene),
        publication, strain, genotype, annotation (association),
        and descriptive notes.
        :param limit: Only parse this many lines of each table
        :return:

        """
        if limit is not None:
            logger.info("Only parsing first %d rows of each file", limit)
        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        # the following will provide us the hash-lookups
        # These must be processed in a specific order
        self._process_prb_strain_acc_view(limit)
        self._process_mrk_acc_view()
        self._process_all_summary_view(limit)
        self._process_bib_acc_view(limit)
        self._process_gxd_genotype_summary_view(limit)

        # The following will use the hash populated above
        # to lookup the ids when filling in the graph
        self._process_prb_strain_view(limit)
        # self._process_prb_strain_genotype_view(limit)
        self._process_gxd_genotype_view(limit)
        self._process_mrk_marker_view(limit)
        self._process_mrk_acc_view_for_equiv(limit)
        self._process_mrk_summary_view(limit)
        self._process_all_allele_view(limit)
        self._process_all_allele_mutation_view(limit)
        self._process_gxd_allele_pair_view(limit)
        self._process_voc_annot_view(limit)
        self._process_voc_evidence_view(limit)
        self._process_mgi_note_vocevidence_view(limit)
        self._process_mrk_location_cache(limit)
        self.process_mgi_relationship_transgene_genes(limit)
        self.process_mgi_note_allele_view(limit)

        logger.info("Finished parsing.")

        logger.info("Loaded %d nodes", len(self.graph))
        return

    def fetch_transgene_genes_from_db(self, cxn):
        """
        This is a custom query to fetch the non-mouse genes that
        are part of transgene alleles.

        :param cxn:
        :return:
        """

        query = "" + \
                "SELECT r._relationship_key as rel_key, " + \
                "r._object_key_1 as object_1, " + \
                "a.accid as allele_id, " + \
                "alabel.label as allele_label, " + \
                "rc._category_key as category_key, " + \
                "rc.name as category_name, " +\
                "t._term_key as property_key, " +\
                "t.term as property_name, " +\
                "rp.value as property_value " +\
                "FROM mgi_relationship r " +\
                "JOIN mgi_relationship_category rc " +\
                "ON r._category_key = rc._category_key " +\
                "JOIN acc_accession a " +\
                "ON r._object_key_1 = a._object_key " +\
                "AND rc._mgitype_key_1 = a._mgitype_key " +\
                "AND a._logicaldb_key = 1 " +\
                "JOIN all_label alabel " +\
                "ON a._object_key = alabel._allele_key " +\
                "AND alabel._label_status_key = 1 " +\
                "AND alabel.priority = 1 " +\
                "JOIN mgi_relationship_property rp " +\
                "ON r._relationship_key = rp._relationship_key " +\
                "AND rp._propertyname_key = 12948292  " +\
                "JOIN voc_term t " +\
                "ON rp._propertyname_key = t._term_key " +\
                "WHERE r._category_key = 1004  "

        self.fetch_query_from_pgdb(
            'mgi_relationship_transgene_genes', query, None, cxn)

        return

    def _process_gxd_genotype_view(self, limit=None):
        """
        This table indicates the relationship between a genotype
        and it's background strain.  It leverages the Genotype class methods
        to do this.

        Makes these triples:
        <MGI:genotypeid> GENO:has_reference_part <MGI:strainid>
        <MGI:strainid> a GENO:genomic_background

        If the genotype id isn't in the hashmap, it adds it here
        (but this shouldn't happen):
        <MGI:genotypeid> a GENO:genotype

        If the strain isn't in the hashmap, it also adds it here with a
        monarchized identifier using the unique key of the strain,
        formatted like:  :_mgistrainkey12345

        :param limit:
        :return:
        """

        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)
        model = Model(g)

        raw = '/'.join((self.rawdir, 'gxd_genotype_view'))
        logger.info("getting genotypes and their backgrounds")
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line = line.rstrip("\n")
                line_counter += 1
                (genotype_key, strain_key, strain, mgiid) = line.split('\t')

                if self.testMode is True:
                    if int(genotype_key) not in self.test_keys.get('genotype'):
                        continue

                if self.idhash['genotype'].get(genotype_key) is None:
                    # just in case we haven't seen it before,
                    # catch and add the id mapping here
                    self.idhash['genotype'][genotype_key] = mgiid
                    geno.addGenotype(mgiid, None)
                    # the label is elsewhere...
                    # need to add the MGI label as a synonym

                # if it's in the hash,
                # assume that the individual was created elsewhere
                strain_id = self.idhash['strain'].get(strain_key)
                background_type = geno.genoparts['genomic_background']
                if strain_id is None or int(strain_key) < 0:
                    if strain_id is None:
                        # some of the strains don't have public identifiers!
                        # so we make one up, and add it to the hash
                        strain_id = self._makeInternalIdentifier('strain',
                                                                 strain_key)
                        self.idhash['strain'].update({strain_key: strain_id})
                        model.addComment(strain_id, "strain_key:"+strain_key)
                    elif int(strain_key) < 0:
                        # these are ones that are unidentified/unknown.
                        # so add instances of each.
                        strain_id = \
                            self._makeInternalIdentifier(
                                'strain', re.sub(r':', '', str(strain_id)))
                        strain_id += re.sub(r':', '', str(mgiid))
                        strain_id = re.sub(r'^_', '_:', strain_id)
                        strain_id = re.sub(r'::', ':', strain_id)
                        model.addDescription(
                            strain_id,
                            "This genomic background is unknown.  " +
                            "This is a placeholder background for " +
                            mgiid + ".")
                        background_type = \
                            geno.genoparts['unspecified_genomic_background']

                    # add it back to the idhash
                    logger.info("adding background as internal id: %s %s: %s",
                                strain_key, strain, strain_id)

                geno.addGenomicBackgroundToGenotype(
                    strain_id, mgiid, background_type)

                self.label_hash[strain_id] = strain

                # add BG to a hash so we can build the genotype label later
                self.geno_bkgd[mgiid] = strain_id

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_gxd_genotype_summary_view(self, limit=None):
        """
        Add the genotype internal id to mgiid mapping to the idhashmap.
        Also, add them as individuals to the graph.
        We re-format the label to put the background strain in brackets
        after the gvc.

        We must pass through the file once to get the ids and
        aggregate the vslcs into a hashmap into the genotype

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
        geno_hash = {}
        raw = '/'.join((self.rawdir, 'gxd_genotype_summary_view'))
        logger.info("building labels for genotypes")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (object_key, preferred, mgiid, subtype,
                 short_description) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('genotype'):
                        continue

                # add the internal genotype to mgi mapping
                self.idhash['genotype'][object_key] = mgiid

                if preferred == '1':
                    d = re.sub(r'\,', '/', short_description.strip())
                    if mgiid not in geno_hash:
                        geno_hash[mgiid] = {'vslcs': [d], 'subtype': subtype,
                                            'key': object_key}
                    else:
                        vslcs = geno_hash[mgiid].get('vslcs')
                        vslcs.append(d)
                else:
                    pass
                    # TODO what to do with != preferred

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        # now, loop through the hash and add the genotypes as individuals
        # we add the mgi genotype as a synonym
        # (we generate our own label later)
        geno = Genotype(g)
        for gt in geno_hash:
            genotype = geno_hash.get(gt)
            gvc = sorted(genotype.get('vslcs'))
            label = '; '.join(gvc) + ' [' + genotype.get('subtype') + ']'
            geno.addGenotype(gt, None)
            model.addComment(gt, self._makeInternalIdentifier(
                'genotype', genotype.get('key')))
            model.addSynonym(gt, label.strip())

        return

    def _process_all_summary_view(self, limit):
        """
        Here, we get the allele definitions: id, label, description, type
        We also add the id to this source's global idhash for lookup later

        <alleleid> a OWL:NamedIndividual
            rdf:label "allele symbol"
            dc:description "long allele name"

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'all_summary_view'))
        logger.info(
            "alleles with labels and descriptions from all_summary_view")
        with open(raw, 'r') as f:
            col_count = f.readline().count('\t')  # read the header row; skip
            # head -1 workspace/build-mgi-ttl/dipper/raw/mgi/all_summary_view|\
            # tr '\t' '\n' | grep -n . | \
            # awk -F':' '{col=$1;$1="";print $0,",\t  #" col}'
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1
                cols = line.count('\t')
                # bail if the row is malformed
                if cols != col_count:
                    logger.warning('Expected ' + str(col_count) + ' columns.')
                    logger.warning('Received ' + str(cols) + ' columns.')
                    logger.warning(line.format())
                    continue
                # no stray tab in the description column
                (object_key, preferred, mgiid, description,
                 short_description) = line.split('\t')
                # NOTE: May want to filter alleles based on the preferred field
                # (preferred = 1) or will get duplicates
                # (24288, to be exact...
                # Reduced to 480 if filtered on preferred = 1)

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('allele'):
                        continue

                # we are setting the allele type to None,
                # so that we can add the type later
                # since we don't actually know
                # if it's a reference or altered allele
                altype = None  # temporary; we'll assign the type later

                # If we want to filter on preferred:
                if preferred == '1':
                    # add the allele key to the hash for later lookup
                    self.idhash['allele'][object_key] = mgiid
                    # TODO consider not adding the individuals in this one
                    model.addIndividualToGraph(
                        mgiid, short_description.strip(),
                        altype, description.strip())
                    self.label_hash[mgiid] = short_description.strip()

                # TODO deal with non-preferreds, are these deprecated?

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_all_allele_view(self, limit):
        """
        Add the allele as a variant locus (or reference locus if wild-type).
        If the marker is specified, we add the link to the marker.
        We assume that the MGI ids are available in the idhash,
        added in all_summary_view.
        We add the sequence alteration as a BNode here, if there is a marker.
        Otherwise, the allele itself is a sequence alteration.

        Triples:
        <MGI:allele_id> a GENO:variant_locus
            OR GENO:reference_locus
            OR GENO:sequence_alteration   IF no marker_id specified.

            [GENO:has_variant_part OR GENO:has_reference_part] <MGI:marker_id>
            GENO:derived_from <MGI:strain_id>
            GENO:has_variant_part <_seq_alt_id>
        <_seq_alt_id> a GENO:sequence_alteration
            derives_from <strain_id>

        :param limit:
        :return:

        """
        # transmission_key -> inheritance? Need to locate related table.
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        geno = Genotype(g)
        line_counter = 0
        logger.info(
            "adding alleles, mapping to markers, " +
            "extracting their sequence alterations " +
            "from all_allele_view")
        raw = '/'.join((self.rawdir, 'all_allele_view'))
        with open(raw, 'r') as f:
            col_count = f.readline().count('\t')  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1
                cols = line.count('\t')
                # bail if the row is malformed
                if cols != col_count:
                    logger.warning('Expected ' + str(col_count) + ' columns.')
                    logger.warning('Received ' + str(cols) + ' columns.')
                    logger.warning(line.format())
                    continue

                (allele_key, marker_key, strain_key, symbol,
                 name, iswildtype) = line.split('\t')

                # TODO update processing to use this view better
                # including jnums!

                if self.testMode is True:
                    if int(allele_key) not in self.test_keys.get('allele'):
                        continue

                allele_id = self.idhash['allele'].get(allele_key)
                if allele_id is None:
                    logger.error(
                        "what to do! can't find allele_id. skipping %s %s",
                        allele_key, symbol)
                    continue

                marker_id = None

                if marker_key is not None and marker_key != '':
                    # we make the assumption here that the markers
                    # have already been added to the table
                    marker_id = self.idhash['marker'].get(marker_key)
                    if marker_id is None:
                        logger.error(
                            "what to do! can't find marker_id. skipping %s %s",
                            marker_key, symbol)
                        continue

                iseqalt_id = self._makeInternalIdentifier('seqalt', allele_key)

                # for non-wild type alleles:
                if iswildtype == '0':
                    locus_type = geno.genoparts['variant_locus']
                    locus_rel = \
                        geno.properties['is_sequence_variant_instance_of']
                # for wild type alleles:
                elif iswildtype == '1':
                    locus_type = geno.genoparts['reference_locus']
                    locus_rel = geno.properties['is_reference_instance_of']
                    # add the allele to the wildtype set for lookup later
                    self.wildtype_alleles.add(allele_id)
                else:
                    locus_rel = None
                    locus_type = None

                model.addIndividualToGraph(allele_id, symbol, locus_type)
                model.makeLeader(allele_id)
                self.label_hash[allele_id] = symbol
                self.idhash['seqalt'][allele_key] = iseqalt_id

                # HACK - if the label of the allele == marker,
                # then make the thing a seq alt
                allele_label = self.label_hash.get(allele_id)
                marker_label = self.label_hash.get(marker_id)
                if allele_label is not None and allele_label == marker_label:
                    model.addSameIndividual(allele_id, marker_id)
                    self.idhash['seqalt'][allele_key] = allele_id
                    model.addComment(
                        allele_id,
                        self._makeInternalIdentifier('allele', allele_key))
                elif marker_id is not None:
                    # marker_id will be none if the allele
                    # is not linked to a marker
                    # (as in, it's not mapped to a locus)
                    geno.addAlleleOfGene(allele_id, marker_id, locus_rel)

                # sequence alteration in strain

                if iswildtype == '0':
                    sa_label = symbol
                    sa_id = iseqalt_id

                    if marker_key is not None and \
                            allele_label != marker_label and marker_key != '':
                        # sequence alteration has label reformatted(symbol)
                        if re.match(r".*<.*>.*", symbol):
                            sa_label = re.sub(r".*<", "<", symbol)
                        elif re.match(r"\+", symbol):
                            # TODO: Check to see if this is the proper handling
                            # as while symbol is just +,
                            # marker symbol has entries without any <+>.
                            sa_label = '<+>'
                        geno.addSequenceAlterationToVariantLocus(iseqalt_id,
                                                                 allele_id)
                    else:
                        # make the sequence alteration == allele
                        sa_id = allele_id

                    # else this will end up adding the non-located transgenes
                    # as sequence alterations also removing the < and > from sa
                    sa_label = re.sub(r'[\<\>]', '', sa_label)

                    # gu.addIndividualToGraph(g,sa_id,sa_label,None,name)
                    geno.addSequenceAlteration(sa_id, sa_label, None, name)
                    self.label_hash[sa_id] = sa_label

                    strain_id = self.idhash['strain'].get(strain_key)
                    # scrub out if the strain is "not specified"
                    if strain_id is not None and \
                            strain_id not in [
                                'MGI:4867032', 'MGI:5649511']:
                        geno.addSequenceDerivesFrom(allele_id, strain_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_gxd_allele_pair_view(self, limit):
        """
        This assumes that the genotype and alleles
        have already been added to the id hashmap.
        We use the Genotype methods to add all the parts we need.
        Triples added:
        <genotype_id> has_part <vslc>
        <vslc> has_part <allele1>
        <vslc> has_part <allele2>
        <vslc> has_zygosity <zygosity>

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        geno = Genotype(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'gxd_allelepair_view'))
        logger.info("processing allele pairs (VSLCs) for genotypes")
        geno_hash = {}
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (allelepair_key, genotype_key, allele_key_1, allele_key_2,
                 allele1, allele2, allelestate) = line.split('\t')
                # NOTE: symbol = gene/marker,
                # allele1 + allele2 = VSLC,
                # allele1/allele2 = variant locus,
                # allelestate = zygosity
                # FIXME Need to handle alleles not in the *<*> format,
                # incl gene traps, induced mut, & transgenics

                if self.testMode is True:
                    if int(genotype_key) not in self.test_keys.get('genotype'):
                        continue

                genotype_id = self.idhash['genotype'].get(genotype_key)
                if genotype_id not in geno_hash:
                    geno_hash[genotype_id] = set()
                if genotype_id is None:
                    logger.error("genotype_id not found for key %s; skipping",
                                 genotype_key)
                    continue

                allele1_id = self.idhash['allele'].get(allele_key_1)
                allele2_id = self.idhash['allele'].get(allele_key_2)

                # Need to map the allelestate to a zygosity term
                zygosity_id = self._map_zygosity(allelestate)
                ivslc_id = self._makeInternalIdentifier('vslc', allelepair_key)

                geno_hash[genotype_id].add(ivslc_id)
                # TODO: VSLC label likely needs processing similar to
                # the processing in the all_allele_view
                # FIXME: handle null alleles
                vslc_label = allele1+'/'
                if allele2_id is None:
                    if zygosity_id in [
                            geno.zygosity['hemizygous'],
                            geno.zygosity['hemizygous-x'],
                            geno.zygosity['hemizygous-y']]:
                        vslc_label += '0'
                    elif zygosity_id == geno.zygosity['heterozygous']:
                        vslc_label += '+'
                    elif zygosity_id == geno.zygosity['indeterminate']:
                        vslc_label += '?'
                    elif zygosity_id == geno.zygosity['homozygous']:
                        # we shouldn't get here, but for testing this is handy
                        vslc_label += allele1
                    else:
                        logger.info(
                            "A different kind of zygosity is found: %s",
                            zygosity_id)
                        vslc_label += '?'
                else:
                    vslc_label += allele2

                model.addIndividualToGraph(
                    ivslc_id, vslc_label,
                    geno.genoparts['variant_single_locus_complement'])
                self.label_hash[ivslc_id] = vslc_label
                rel1 = rel2 = geno.object_properties['has_alternate_part']
                if allele1_id in self.wildtype_alleles:
                    rel1 = geno.object_properties['has_reference_part']
                if allele2_id in self.wildtype_alleles:
                    rel2 = geno.object_properties['has_reference_part']
                geno.addPartsToVSLC(
                    ivslc_id, allele1_id, allele2_id, zygosity_id, rel1, rel2)

                # if genotype_id not in geno_hash:
                #     geno_hash[genotype_id] = [vslc_label]
                # else:
                #     geno_hash[genotype_id] += [vslc_label]

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        # build the gvc and the genotype label
        for gt in geno_hash.keys():
            if gt is None:   # not sure why, but sometimes this is the case
                continue
            vslcs = sorted(list(geno_hash[gt]))
            gvc_label = None
            if len(vslcs) > 1:
                gvc_id = re.sub(r'_', '', ('-'.join(vslcs)))
                gvc_id = re.sub(r':', '', gvc_id)
                gvc_id = '_:'+gvc_id
                vslc_labels = []
                for v in vslcs:
                    vslc_labels.append(self.label_hash[v])
                gvc_label = '; '.join(vslc_labels)

                model.addIndividualToGraph(
                    gvc_id, gvc_label,
                    geno.genoparts['genomic_variation_complement'])
                self.label_hash[gvc_id] = gvc_label
                for v in vslcs:
                    geno.addParts(
                        v, gvc_id,
                        geno.object_properties['has_alternate_part'])
                    geno.addVSLCtoParent(v, gvc_id)
                geno.addParts(
                    gvc_id, gt,
                    geno.object_properties['has_alternate_part'])
            elif len(vslcs) == 1:
                gvc_id = vslcs[0]
                gvc_label = self.label_hash[gvc_id]
                # type the VSLC as also a GVC
                model.addIndividualToGraph(
                    gvc_id, gvc_label,
                    geno.genoparts['genomic_variation_complement'])
                geno.addVSLCtoParent(gvc_id, gt)
            else:
                logger.info("No VSLCs for %s", gt)

            # make the genotype label = gvc + background
            bkgd_id = self.geno_bkgd.get(gt)
            if bkgd_id is not None:
                bkgd_label = self.label_hash.get(bkgd_id)
                if bkgd_label is None:
                    bkgd_label = bkgd_id  # just in case
            else:
                bkgd_label = 'n.s.'
            if gvc_label is not None:
                genotype_label = gvc_label + ' ['+bkgd_label+']'
            else:
                genotype_label = '['+bkgd_label+']'

            model.addIndividualToGraph(gt, genotype_label)
            self.label_hash[gt] = genotype_label

        return

    def _process_all_allele_mutation_view(self, limit):
        """
        This fetches the mutation type for the alleles,
        and maps them to the sequence alteration.
        Note that we create a BNode for the sequence alteration because
        it isn't publicly identified.
        <sequence alteration id> a <SO:mutation_type>

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'all_allele_mutation_view'))
        logger.info("getting mutation types for sequence alterations")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (allele_key, mutation) = line.split('\t')
                iseqalt_id = self.idhash['seqalt'].get(allele_key)
                if iseqalt_id is None:
                    iseqalt_id = \
                        self._makeInternalIdentifier(
                            'seqalt', allele_key)

                if self.testMode and \
                        int(allele_key) not in \
                        self.test_keys.get('allele'):
                    continue

                # TODO we might need to map the seq alteration to the MGI id
                # for unlocated things; need to use hashmap
                # map the sequence_alteration_type
                seq_alt_type_id = self._map_seq_alt_type(mutation)
                # HACK - if the seq alteration is a transgene,
                # then make sure it is a transgenic insertion
                allele_id = self.idhash['allele'].get(allele_key)
                if allele_id is not None:
                    allele_label = self.label_hash.get(allele_id)
                    if allele_label is not None and re.search(r'Tg\(',
                                                              allele_label):
                        logger.info("Found a transgenic insertion for %s",
                                    allele_label)
                        # transgenic_insertion, instead of plain old insertion
                        seq_alt_type_id = 'SO:0001218'

                model.addIndividualToGraph(iseqalt_id, None, seq_alt_type_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_voc_annot_view(self, limit):
        """
        This MGI table represents associations between things.
        We now filter this table on abnormal Genotype-Phenotype associations,
        but may be expanded in the future.

        We add the internal annotation id to the idhashmap.
        It is expected that the genotypes have already been added to the idhash

        :param limit:
        :return:

        """

        # TODO also get Strain/Attributes (annottypekey = 1000)
        # TODO what is Phenotype (Derived) vs
        # non-derived?  (annottypekey = 1015)
        # TODO is evidence in this table?  what is the evidence vocab key?
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        logger.info("getting G2P associations")
        raw = '/'.join((self.rawdir, 'voc_annot_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")

                (annot_key, annot_type, object_key, term_key,
                 qualifier_key, qualifier, term, accid) = line.split('\t')

                if self.testMode is True:
                    if int(annot_key) not in self.test_keys.get('annot'):
                        continue

                # iassoc_id = self._makeInternalIdentifier('annot', annot_key)
                # assoc_id = self.make_id(iassoc_id)

                assoc_id = None
                # Mammalian Phenotype/Genotype are curated G2P assoc
                if annot_type == 'Mammalian Phenotype/Genotype':
                    line_counter += 1

                    # TODO add NOT annotations
                    # skip 'normal'
                    if qualifier == 'norm':
                        logger.info("found normal phenotype: %s", term)
                        continue

                    # We expect the label for the phenotype
                    # to be taken care of elsewhere
                    model.addClassToGraph(accid, None)

                    genotype_id = self.idhash['genotype'].get(object_key)
                    if genotype_id is None:
                        logger.error("can't find genotype id for %s",
                                     object_key)
                    else:
                        # add the association
                        assoc = G2PAssoc(g, self.name, genotype_id, accid)
                        assoc.add_association_to_graph()
                        assoc_id = assoc.get_association_id()
                # OMIM/Genotype are disease-models
                elif annot_type == 'DO/Genotype':
                    # skip NOT annotations for now FIXME
                    if qualifier_key == '1614157':
                        continue
                    genotype_id = self.idhash['genotype'].get(object_key)
                    if genotype_id is None:
                        logger.error("can't find genotype id for %s",
                                     object_key)
                    else:
                        # add the association
                        assoc = Assoc(g, self.name)
                        # TODO PYLINT
                        # Redefinition of assoc type from
                        # dipper.models.assoc.G2PAssoc.G2PAssoc to
                        # dipper.models.assoc.Association.Assoc
                        assoc.set_subject(genotype_id)
                        assoc.set_object(accid)
                        assoc.set_relationship(
                            model.object_properties['model_of'])
                        assoc.add_association_to_graph()
                        assoc_id = assoc.get_association_id()
                elif annot_type == 'MCV/Marker':
                    # marker category == type
                    marker_id = self.idhash['marker'].get(object_key)
                    term_id = self._map_marker_category(str(term_key))
                    # note that the accid here is an internal mouse cv term,
                    # and we don't use it.
                    if term_id is not None and marker_id is not None:
                        # do something special for transgenics -
                        # make sure these are transgenic insertions
                        model.addType(marker_id, term_id)
                elif annot_type == 'DO/Allele':  # allele/Disease
                    allele_id = self.idhash['allele'].get(object_key)
                    if allele_id is None:
                        logger.error("can't find genotype id for %s",
                                     object_key)
                    else:
                        # add the association
                        assoc = Assoc(g, self.name)
                        assoc.set_subject(allele_id)
                        assoc.set_object(accid)
                        assoc.set_relationship(
                            model.object_properties['model_of'])
                        assoc.add_association_to_graph()
                        assoc_id = assoc.get_association_id()

                if assoc_id is not None:
                    # add the assoc to the hashmap (using the monarch id)
                    self.idhash['annot'][annot_key] = assoc_id
                    model.addComment(assoc_id, "annot_key:"+annot_key)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_voc_evidence_view(self, limit):
        """
        Here we fetch the evidence (code and publication) for the associations.
        The evidence codes are mapped from the standard GO codes to ECO.
        J numbers are added for publications.
        We will only add the evidence if the annotation is in our idhash.

        Triples:
        <annot_id> dc:evidence <evidence_id>
        <pub_id> a owl:NamedIndividual
        <annot_id> dc:source <pub_id>

        :param limit:
        :return:

        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        logger.info("getting evidence and pubs for annotations")
        raw = '/'.join((self.rawdir, 'voc_evidence_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (annot_evidence_key, annot_key,
                 evidence_code, jnumid) = line.split('\t')

                if self.testMode is True:
                    if int(annot_key) not in self.test_keys.get('annot'):
                        continue

                # add the association id to map to the evidence key
                # (to attach the right note to the right assn)
                self.idhash['notes'][annot_evidence_key] = annot_key

                assoc_id = self.idhash['annot'].get(annot_key)

                if assoc_id is None:
                    # assume that we only want to add the evidence/source
                    # for annots that we have in our db
                    continue

                evidence_id = self._map_evidence_id(evidence_code)

                reference = Reference(g, jnumid)
                reference.addRefToGraph()

                # add the ECO and citation information to the annot
                g.addTriple(assoc_id,
                            Assoc.object_properties['has_evidence'],
                            evidence_id)
                g.addTriple(assoc_id,
                            Assoc.object_properties['has_source'],
                            jnumid)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_bib_acc_view(self, limit):
        """
        This traverses the table twice:
        once to look up the internal key to J number mapping
        for the id hashmap then again to make the equivalences.
        All internal keys have both a J and MGI identifier.
        This will make equivalences between the different pub ids
        Triples:
            <pub_id> a owl:NamedIndividual
            <other_pub_id> a owl:NamedIndividual
            <pub_id> owl:sameAs <other_pub_id>
        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        # firstpass, get the J number mapping, and add to the global hash
        line_counter = 1
        logger.info('populating pub id hash')
        raw = '/'.join((self.rawdir, 'bib_acc_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            header = next(filereader)
            if len(header) != 6:
                logger.error('bib_acc_view expected 6 columns got: %s', header)
            for line in filereader:
                line_counter += 1
                (accid, prefixpart, numericpart, object_key,
                 logical_db, logicaldb_key) = line

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('pub'):
                        continue

                # we use the J number here because
                # it is the externally-accessible identifier
                if prefixpart != 'J:':
                    continue
                self.idhash['publication'][object_key] = accid
                reference = Reference(g, accid)
                reference.addRefToGraph()

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        # 2nd pass, look up the MGI identifier in the hash
        logger.info("getting pub equivalent ids")
        line_counter = 1
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            header = next(filereader)

            for line in filereader:
                line_counter += 1
                (accid, prefixpart, numericpart, object_key,
                 logical_db, logicaldb_key) = line

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('pub'):
                        continue

                logical_db = logical_db.strip()

                jid = self.idhash['publication'].get(object_key)

                pub_id = None
                if logicaldb_key == '29':  # pubmed
                    pub_id = 'PMID:'+accid
                elif logicaldb_key == '1' and re.match(r'MGI:', prefixpart):
                    # don't get the J numbers,
                    # because we dont' need to make the equiv to itself.
                    pub_id = accid
                elif logical_db == 'Journal Link':
                    # some DOIs seem to have spaces
                    # FIXME MGI needs to FIX THESE UPSTREAM!!!!
                    # we'll scrub them here for the time being
                    accid = re.sub(r'\s+', '', accid)
                    # some DOIs have un-urlencoded brackets <>
                    accid = re.sub(r'<', '%3C', accid)
                    accid = re.sub(r'>', '%3E', accid)
                    pub_id = 'DOI:' + accid

                elif logicaldb_key == '1' and re.match(r'J:', prefixpart):
                    # we can skip the J numbers
                    continue

                if pub_id is not None:
                    # only add these to the graph if
                    # it's mapped to something we understand
                    reference = Reference(g, pub_id)

                    # make the assumption that if it is a PMID, it is a journal
                    if re.match(r'PMID', pub_id):
                        reference.setType(
                            Reference.ref_types['journal_article'])
                        model.makeLeader(pub_id)
                    reference.addRefToGraph()

                    model.addSameIndividual(jid, pub_id)
                else:
                    logger.warning("Publication from (%s) not mapped for %s",
                                   logical_db, object_key)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_prb_strain_view(self, limit):
        """
        Process a table to get strains (with internal ids), and their labels.
        These strains are created as instances of the species that they are.
        Triples:
            <strain id> a GENO:intrinsic_genotype
                rdf:label "strain label"
                RO:in_taxon <NCBI taxon id>

        :param limit:
        :return:

        """
        # Only 9 strain types if we want to map them
        #   recombinant congenci,   inbred strain,  NA,
        #   congenic,               consomic,       coisogenic,
        #   recombinant inbred,     NS,             conplastic

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        line_counter = 0
        geno = Genotype(g)
        raw = '/'.join((self.rawdir, 'prb_strain_view'))
        logger.info("getting strains and adding their taxa")
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue
                (strain_key, strain, species) = line

                if self.testMode is True:
                    if int(strain_key) not in self.test_keys.get('strain'):
                        continue

                strain_id = self.idhash['strain'].get(strain_key)

                if strain_id is not None:
                    self.label_hash[strain_id] = strain

                    # add the species to the graph as a class
                    sp = self._map_strain_species(species)
                    if sp is not None:
                        model.addClassToGraph(sp, None)
                        geno.addTaxon(sp, strain_id)
                    model.addIndividualToGraph(strain_id, strain, sp)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_mrk_marker_view(self, limit):
        """
        This is the definition of markers
        (as in genes, but other genomic loci types as well).
        It looks up the identifiers in the hashmap
        This includes their labels, specific class, and identifiers
        TODO should we use the mrk_mouse_view instead?

        Triples:
        <marker_id> a owl:Class OR owl:NamedIndividual
            GENO:marker_type
            rdf:label <symbol>
            RO:in_taxon <NCBITaxon_id>

        :param limit:
        :return:
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        geno = Genotype(g)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'mrk_marker_view'))
        logger.info("getting markers and assigning types")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (marker_key, organism_key, marker_status_key,
                 symbol, name, latin_name, marker_type) = line.split('\t')

                if self.testMode is True:
                    if int(marker_key) not in self.test_keys.get('marker'):
                        continue

                # use only non-withdrawn markers
                if marker_status_key != '2':
                    marker_id = self.idhash['marker'].get(marker_key)

                    # only pull info for mouse genes for now
                    # other species should come from other dbs
                    if organism_key != '1':
                        continue

                    if marker_id is None:
                        logger.error("can't find %s %s in the id hash",
                                     marker_key, symbol)

                    mapped_marker_type = self._map_marker_type(marker_type)

                    # if it's unlocated, or is not a gene,
                    # then don't add it as a class because
                    # it's not added as a gene.
                    # everything except for genes are modeled as individuals

                    # it's a gene or pseudogene
                    if mapped_marker_type in ['SO:0000704', 'SO:0000336']:
                        model.addClassToGraph(
                            marker_id, symbol, mapped_marker_type, name)
                        model.addSynonym(
                            marker_id,
                            name, Assoc.properties['hasExactSynonym'])
                        self.markers['classes'].append(marker_id)
                    else:
                        model.addIndividualToGraph(
                            marker_id, symbol, mapped_marker_type, name)
                        model.addSynonym(
                            marker_id, name,
                            Assoc.properties['hasExactSynonym'])
                        self.markers['indiv'].append(marker_id)

                    self.label_hash[marker_id] = symbol
                    # add the taxon
                    taxon_id = self._map_taxon(latin_name)
                    geno.addTaxon(taxon_id, marker_id)

                    # make MGI the leader for mouse genes.
                    if taxon_id == 'NCBITaxon:10090':
                        model.makeLeader(marker_id)

                    if not self.testMode and \
                            limit is not None and line_counter > limit:
                        break

        return

    def _process_mrk_summary_view(self, limit):
        """
        Here we pull the mgiid of the features, and make equivalent (or sameAs)
        associations to referenced ids.
        Only adding the ENSEMBL genes and NCBI gene ids.
        Will wait on other ids later.

        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("getting markers and equivalent ids from mrk_summary_view")
        line_counter = 0
        raw = '/'.join((self.rawdir, 'mrk_summary_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (accid, logicaldb_key, object_key, preferred,
                 mgiid, subtype, short_description) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                if preferred == '1':

                    if self.idhash['marker'].get(object_key) is None:
                        # can't find the marker in the hash; add it here:
                        self.idhash['marker'][object_key] = mgiid
                        logger.error(
                            "this marker hasn't been seen before %s %s",
                            mgiid, short_description)

                    if accid == mgiid:
                        # don't need to make equivalences to itself
                        continue

                    mapped_id = None
                    if logicaldb_key == '60':
                        mapped_id = 'ENSEMBL:'+accid
                    elif logicaldb_key == '1':
                        # don't need to add the equivalence to itself.
                        continue
                    elif logicaldb_key == '55':
                        mapped_id = 'NCBIGene:'+accid

                    if mapped_id is not None:
                        if mgiid in self.markers['classes'] or \
                                subtype in ['Gene', 'Pseudogene']:
                            model.addClassToGraph(mapped_id, None)
                            model.addEquivalentClass(mgiid, mapped_id)
                        elif mgiid in self.markers['indiv']:
                            model.addIndividualToGraph(mapped_id, None)
                            model.addSameIndividual(mgiid, mapped_id)

                    # could parse the "subtype" string
                    # to get the kind of thing the marker is

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_mrk_acc_view(self):
        """
        Use this table to create the idmap between the internal marker id and
        the public mgiid.
        No triples are produced in this process
        :return:

        """

        # make a pass through the table first,
        # to create the mapping between the external and internal identifiers
        line_counter = 0
        logger.info("mapping markers to internal identifiers")
        raw = '/'.join((self.rawdir, 'mrk_acc_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip('\n')
                line_counter += 1
                (accid, prefix_part, logicaldb_key, object_key,
                 preferred, organism_key) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                # get the hashmap of the identifiers
                if logicaldb_key == '1' and \
                        prefix_part == 'MGI:' and preferred == '1':
                    self.idhash['marker'][object_key] = accid

        return

    def _process_mrk_acc_view_for_equiv(self, limit):
        """
        Add the equivalences, either sameAs or equivalentClass,
        depending on the nature of the marker.
        We only process the ENSEMBL genes and NCBI gene ids.
        :param limit:
        :return:

        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        # pass through the file again,
        # and make the equivalence statements to a subset of the idspaces.
        # TODO verify the difference between what the
        # mrk_acc_view vs mrk_summary_view buys us here.
        # if nothing, then we should remove one or the other.
        logger.info("mapping marker equivalent identifiers in mrk_acc_view")
        line_counter = 0
        with open('/'.join((self.rawdir, 'mrk_acc_view')), 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1
                (accid, prefix_part, logicaldb_key, object_key,
                 preferred, organism_key) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                # right now not caring about other organisms
                if organism_key != 1:
                    continue

                mgiid = self.idhash['marker'].get(object_key)
                if mgiid is None:
                    # presumably we've already added the relevant MGI ids,
                    # so skip those that we can't find
                    logger.debug("can't find mgiid for %s", object_key)
                    continue
                marker_id = None
                if preferred == '1':  # TODO what does it mean if it's 0?
                    if logicaldb_key == '55':  # entrez/ncbi
                        marker_id = 'NCBIGene:'+accid
                    elif logicaldb_key == '1' and prefix_part != 'MGI:':  # mgi
                        marker_id = accid
                    elif logicaldb_key == '60':
                        marker_id = 'ENSEMBL:'+accid
                    # TODO get non-preferred ids==deprecated?

                if marker_id is not None:
                    if mgiid in self.markers['classes']:
                        model.addClassToGraph(marker_id, None)
                        model.addEquivalentClass(mgiid, marker_id)
                    elif mgiid in self.markers['indiv']:
                        model.addIndividualToGraph(marker_id, None)
                        model.addSameIndividual(mgiid, marker_id)
                    else:
                        logger.error("mgiid not in class or indiv hash %s",
                                     mgiid)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_prb_strain_acc_view(self, limit):
        """
        Use this table to create the idmap between
        the internal marker id and the public mgiid.
        Also, add the equivalence statements between strains for MGI and JAX
        Triples:
        <strain_id> a GENO:intrinsic_genotype
        <other_strain_id> a GENO:intrinsic_genotype
        <strain_id> owl:sameAs <other_strain_id>

        :param limit:
        :return:

        """

        # make a pass through the table first,
        # to create the mapping between the external and internal identifiers
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("mapping strains to internal identifiers")
        raw = '/'.join((self.rawdir, 'prb_strain_acc_view'))

        tax_id = 'NCBITaxon:10090'  # hardcode mouse

        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1
                (accid, prefixpart, logicaldb_key,
                 object_key, preferred) = line.split('\t')
                # scrub out the backticks from accids
                # TODO notify the source upstream
                accid = re.sub(r'`', '', accid).strip()
                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('strain'):
                        continue

                # get the hashmap of the identifiers
                if logicaldb_key == '1' and \
                        prefixpart == 'MGI:' and preferred == '1':
                    self.idhash['strain'][object_key] = accid
                    model.addIndividualToGraph(accid, None, tax_id)

        # The following are the stock centers for the strains
        # (asterisk indicates complete)
        # *1	MGI	    Mouse Genome Informatics
        # *22	JAX Registry	(null)
        # *37	EMMA	European Mutant Mouse Archive
        # *38	MMRRC	Mutant Mouse Regional Resource Center
        # 39	Harwell	Mammalian Genome Unit Stock List
        # *40	ORNL	Oak Ridge National Lab mutant resource
        # *54	NCIMR	NCI Mouse Repository
        # *56	NMICE	Neuromice.org, a consortium of three NIH-sponsored
        #                   mutagenesis projects designed to search for
        #                   neurological mutations
        # 57	CARD	Center for Animal Resources and Development @ Kumamoto U
        # *70	RIKEN BRC	RIKEN BioResource Center
        # *71	CMMR	Canadian Mouse Mutant Resource
        # 84	JPGA	The Center for New Mouse Models of
        #                   Heart, Lung, BLood and Sleep Disorders,
        #                   JAX-PGA at The Jackson Laboratory
        # *87	MUGEN	Network of Excellence in Integrated Functional Genomics
        #                   in Mutant Mouse Models as Tools to Investigate the
        #                   Complexity of Human Immunological Disease
        # *90	APB	    Australian Phenomics Bank
        # ? 91	EMS	    Elizabeth M. Simpson
        # ? 93	NIG 	National Institute of Genetics,
        #                   Mammalian Genetics Laboratory, Japan
        # 94	TAC	    Taconic
        # 154	OBS 	Oriental BioService , Inc.
        # 161	RMRC-NLAC	National Applied Research Laboratories,Taiwan, R.O.C.

        # pass through the file again,
        # and make the equivalence statements to a subset of the idspaces
        logger.info("mapping strain equivalent identifiers")
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1
                (accid, prefixpart, logicaldb_key,
                 object_key, preferred) = line.split('\t')
                # scrub out the backticks from accids
                # TODO notify the source upstream
                accid = re.sub(r'`', '', accid).strip()
                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('strain'):
                        continue
                mgiid = self.idhash['strain'].get(object_key)
                if mgiid is None:
                    # presumably we've already added the relevant MGI ids,
                    # so skip those that we can't find
                    # logger.info("can't find mgiid for %s",object_key)
                    continue
                strain_id = None
                deprecated = False
                comment = None
                if preferred == '1':  # what does it mean if it's 0?
                    if logicaldb_key == '22':  # JAX
                        # scrub out the backticks from accids
                        # TODO notify the source upstream
                        accid = re.sub(r'`', '', accid).strip()
                        strain_id = 'JAX:' + accid
                    elif logicaldb_key == '38':  # MMRRC
                        strain_id = accid
                        if not re.match(r'MMRRC:', strain_id):
                            strain_id = 'MMRRC:'+strain_id
                    elif logicaldb_key == '37':  # EMMA
                        strain_id = re.sub(r'EM:', 'EMMA:', accid)
                    elif logicaldb_key == '90':  # APB
                        strain_id = 'APB:' + accid  # Check
                    elif logicaldb_key == '40':  # ORNL
                        # ORNL is not in existence any more.
                        # these are deprecated, and we will prefix with JAX
                        strain_id = 'JAX:' + accid
                        comment = "Originally from ORNL."
                        deprecated = True
                        # add these as synonyms of the MGI mouse
                        model.addSynonym(mgiid, accid)

                    elif logicaldb_key == '54':  # NCIMR
                        strain_id = 'NCIMR:'+accid
                    # CMMR  not great - doesn't resolve well
                    # elif logicaldb_key == '71':
                    #     strain_id = 'CMMR:'+accid
                    elif logicaldb_key == '56':  # neuromice
                        # neuromice.org doesn't exist any more.
                        # but all these are actually MGI ids
                        strain_id = accid
                    elif logicaldb_key == '70':  # RIKEN
                        # like
                        # http://www2.brc.riken.jp/lab/animal/detail.php?brc_no=RBRC00160
                        strain_id = 'RBRC:' + accid
                    elif logicaldb_key == '87':
                        strain_id = 'MUGEN:' + accid
                        # I can't figure out how to get to some of the strains

                    # TODO get non-preferred ids==deprecated?

                # TODO make these strains, rather than instance of taxon?
                if strain_id is not None:
                    model.addIndividualToGraph(strain_id, None, tax_id)
                    if deprecated:
                        model.addDeprecatedIndividual(strain_id, [mgiid])
                        model.addSynonym(mgiid, accid)
                    else:
                        model.addSameIndividual(mgiid, strain_id)
                    if re.match(r'MMRRC', strain_id):
                        model.makeLeader(strain_id)
                    if comment is not None:
                        model.addComment(strain_id, comment)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_mgi_note_vocevidence_view(self, limit):
        """
        Here we fetch the free text descriptions of the phenotype associations.
        Triples:
        <annot_id> dc:description "description text"
        :param limit:
        :return:

        """

        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("getting free text descriptions for annotations")
        raw = '/'.join((self.rawdir, 'mgi_note_vocevidence_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (object_key, note) = line

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('notes'):
                        continue
                # object_key == evidence._annotevidence_key
                annotkey = self.idhash['notes'].get(object_key)
                annot_id = self.idhash['annot'].get(annotkey)
                # only add the description for the annotations
                # we have captured through processing

                if annot_id is not None:
                    model.addDescription(annot_id, note.strip())

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def _process_mrk_location_cache(self, limit):
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        # model = Model(g)  # unused
        logger.info("getting marker locations")
        raw = '/'.join((self.rawdir, 'mrk_location_cache'))
        geno = Genotype(g)

        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (marker_key, organism_key, chromosome, startcoordinate,
                 endcoordinate, strand, version) = line

                # only get the location information for mouse
                if str(organism_key) != '1' or str(chromosome) == 'UN':
                    continue

                if self.testMode is True:
                    if int(marker_key) not in self.test_keys.get('marker'):
                        continue

                # make the chromsomome, and the build-instance
                chrom_id = makeChromID(chromosome, 'NCBITaxon:10090', 'CHR')
                if version is not None and \
                        version != '' and version != '(null)':

                    # switch on maptype or mapkey
                    assembly = version
                    build_id = 'NCBIGenome:'+assembly
                    geno.addChromosomeInstance(chromosome, build_id, assembly,
                                               chrom_id)
                    chrom_id = makeChromID(chromosome, build_id, 'MONARCH')

                if marker_key in self.idhash['marker']:
                    gene_id = self.idhash['marker'][marker_key]
                    feature = Feature(g, gene_id, None, None)
                    if strand == '(null)' or strand == '':
                        strand = None
                    if startcoordinate == '(null)' or startcoordinate == '':
                        startcoordinate = None
                    if endcoordinate == '(null)' or endcoordinate == '':
                        endcoordinate = None

                    if startcoordinate is not None:
                        feature.addFeatureStartLocation(
                            int(float(startcoordinate)), chrom_id, strand)
                    else:
                        feature.addFeatureStartLocation(
                            startcoordinate, chrom_id, strand,
                            [Feature.types['FuzzyPosition']])
                    if endcoordinate is not None:
                        feature.addFeatureEndLocation(
                            int(float(endcoordinate)), chrom_id, strand)
                    # note we don't add the uncertain end coordinate,
                    # because we don't know what it is.
                    add_as_class = False
                    if gene_id in self.markers['classes']:
                        add_as_class = True
                    feature.addFeatureToGraph(True, None, add_as_class)

                else:
                    logger.warning('marker key %s not in idhash',
                                   str(marker_key))

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def process_mgi_relationship_transgene_genes(self, limit=None):
        """
        Here, we have the relationship between MGI transgene alleles,
        and the non-mouse gene ids that are part of them.
        We augment the allele with the transgene parts.

        :param limit:
        :return:

        """
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        logger.info("getting transgene genes")
        raw = '/'.join((self.rawdir, 'mgi_relationship_transgene_genes'))
        geno = Genotype(g)

        # gu = GraphUtils(curie_map.get())  # TODO unused
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (rel_key, allele_key, allele_id, allele_label, category_key,
                 category_name, property_key, property_name, gene_num) = line

                if self.testMode is True:
                    if int(allele_key) not in self.test_keys.get('allele')\
                            and int(gene_num) not in self.test_ids:
                        continue

                gene_id = 'NCBIGene:'+gene_num

                # geno.addParts(gene_id, allele_id,
                #               geno.object_properties['has_alternate_part'])
                seqalt_id = self.idhash['seqalt'].get(allele_key)
                if seqalt_id is None:
                    seqalt_id = allele_id
                geno.addSequenceDerivesFrom(seqalt_id, gene_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    def process_mgi_note_allele_view(self, limit=None):
        """
        These are the descriptive notes about the alleles.
        Note that these notes have embedded HTML -
        should we do anything about that?
        :param limit:
        :return:

        """

        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        logger.info("Assembling notes on alleles")
        raw = '/'.join((self.rawdir, 'mgi_note_allele_view'))
        # geno = Genotype(g)  # TODO unused

        notehash = {}
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (object_key, notetype, note, sequencenum) = line

                # read all the notes into a hash to concatenate
                if object_key not in notehash:
                    notehash[object_key] = {}
                if notetype not in notehash[object_key]:
                    notehash[object_key][notetype] = []
                if len(notehash[object_key][notetype]) < int(sequencenum):
                    for i in range(
                            len(notehash[object_key][notetype]),
                            int(sequencenum)):
                        notehash[object_key][notetype].append('')
                notehash[object_key][notetype][int(sequencenum)-1] = \
                    note.strip()

            # finish iteration over notes

        line_counter = 0
        for allele_key in notehash:
            if self.testMode is True:
                if int(allele_key) not in self.test_keys.get('allele'):
                    continue
            line_counter += 1
            allele_id = self.idhash['allele'].get(allele_key)
            if allele_id is None:
                continue
            for n in notehash[allele_key]:
                logger.info(
                    "found %d %s notes for %s",
                    len(notehash[allele_key]), n, allele_id)
                notes = ''.join(notehash[allele_key][n])
                notes += ' ['+n+']'
                model.addDescription(allele_id, notes)

            if not self.testMode and \
                    limit is not None and line_counter > limit:
                break

        return

    def _process_prb_strain_genotype_view(self, limit=None):
        """
        Here we fetch the free text descriptions of the phenotype associations.
        Triples:
        <annot_id> dc:description "description text"
        :param limit:

        :return:
        """

        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        # model = Model(g)  # unused
        logger.info("Getting genotypes for strains")
        raw = '/'.join((self.rawdir, 'prb_strain_genotype_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (strain_key, genotype_key) = line

                if self.testMode is True:
                    if int(genotype_key) not in \
                            self.test_keys.get('genotype') and \
                            int(strain_key) not in \
                            self.test_keys.get('strain'):
                        continue

                strain_id = self.idhash['strain'].get(strain_key)
                if strain_id is None:
                    strain_id = self._makeInternalIdentifier(
                        'strain', strain_key)
                genotype_id = self.idhash['genotype'].get(genotype_key)
                if genotype_id is None:
                    genotype_id = self._makeInternalIdentifier(
                        'genotype', genotype_key)

                if strain_id is not None and genotype_id is not None:
                    self.strain_to_genotype_map[strain_id] = genotype_id

                g.addTriple(
                    strain_id,
                    Genotype.object_properties['has_genotype'],
                    genotype_id)
                # TODO
                # verify if this should be contingent on the exactness or not
                # if qualifier == 'Exact':
                #     gu.addTriple(
                #       g, strain_id,
                #       Genotype.object_properties['has_genotype'],
                #       genotype_id)
                # else:
                #     gu.addXref(g, strain_id, genotype_id)

                if not self.testMode and \
                        limit is not None and line_counter > limit:
                    break

        return

    # TODO: Finish identifying SO/GENO terms
    # for mappings for those found in MGI
    @staticmethod
    def _map_seq_alt_type(sequence_alteration_type):
        seqalttype = 'SO:0001059'  # default to sequence_alteration
        type_map = {
            'Deletion': 'SO:0000159',  # deletion
            # insertion - correct?
            'Disruption caused by insertion of vector': 'SO:0000667',
            'Duplication': 'SO:1000035',  # duplication
            'Insertion': 'SO:0000667',  # insertion
            # transgenic insertion - correct?
            # TODO gene_trap_construct: SO:0001477
            'Insertion of gene trap vector': 'SO:0001218',
            # deletion  # TODO return a list? SO:0001628 intergenic_variant
            'Intergenic deletion': 'SO:0000159',
            # deletion  # TODO return a list?  SO:0001564 gene_variant
            'Intragenic deletion': 'SO:0000159',
            'Inversion': 'SO:1000036',  # inversion
            'Not Applicable': 'SO:0001059',
            'Not Specified': 'SO:0001059',
            # tandem duplication  # TODO ask for another term
            'Nucleotide repeat expansion': 'SO:1000039',
            # multiple nucleotide variant
            'Nucleotide substitutions': 'SO:0002007',
            'Other': 'SO:0001059',
            'Single point mutation': 'SO:1000008',  # point_mutation
            'Translocation': 'SO:0000199',  # translocation
            'Transposon insertion': 'SO:0001837',  # mobile element insertion
            'Undefined': 'SO:0001059',
            # novel sequence insertion (no viral version)
            'Viral insertion': 'SO:0001838',
            'wild type': 'SO:0000817'  # wild type
        }
        if sequence_alteration_type.strip() in type_map:
            seqalttype = type_map.get(sequence_alteration_type.strip())
        else:
            logger.error(
                "Sequence Alteration Type (%s) not mapped; " +
                "defaulting to sequence_alteration",
                sequence_alteration_type)

        return seqalttype

    @staticmethod
    def _map_marker_category(marker_type_key):
        """
        These map the internal "mouse CV" terms for marker categories
        to SO classes.
        There remains one open ticket to satisfy the lncRNA terms
            that don't have a specific type:
        https://sourceforge.net/p/song/term-tracker/436/
        :param marker_type_key:
        :return:

        """

        so_id = None
        marker_category_to_so_map = {
            '6238160': 'SO:0000704',  # gene

            '15406202': 'SO:0001263',  # lncRNA gene --> ncRNA gene FIXME
            '6238161': 'SO:0001217',  # protein coding gene
            '6238162': 'SO:0001263',  # non-coding RNA gene
            # antisense lncRNA gene --> ncRNA gene FIXME
            '15406203': 'SO:0001263',
            '6238163': 'SO:0001637',  # rRNA gene
            '6238164': 'SO:0001272',  # tRNA gene
            '6238165': 'SO:0001268',  # snRNA gene
            '6238166': 'SO:0001267',  # snoRNA gene
            '6238167': 'SO:0001265',  # miRNA gene
            '6238168': 'SO:0001266',  # scRNA gene
            '6238169': 'SO:0001641',  # lincRNA gene
            '6238170': 'SO:0001500',  # heritable phenotypic marker
            '6238171': 'SO:3000000',  # gene segment
            '7313348': 'SO:0000336',  # pseudogene
            '6238180': 'SO:0001269',  # SRP RNA gene
            '6238181': 'SO:0001639',  # RNase P RNA gene
            '6238182': 'SO:0001640',  # RNase MRP RNA gene
            '6238183': 'SO:0001643',  # telomerase RNA gene
            '6238184': 'SO:0000704',  # unclassified gene  --> gene
            '6238186': 'SO:0001263',  # unclassified non-coding RNA gene
            '6967235': 'SO:0001741',  # pseudogenic gene segment
            '7196768': 'SO:1000029',  # chromosomal deletion
            '7196769': 'SO:0000667',  # insertion
            '7196770': 'SO:1000030',  # chromosomal inversion
            '7196771': 'SO:1000043',  # Robertsonian fusion
            '7196772': 'SO:1000048',  # reciprocal chromosomal translocation
            '7196773': 'SO:1000044',  # chromosomal translocation
            '7196774': 'SO:1000037',  # chromosomal duplication
            '7196775': 'SO:0000453',  # chromosomal transposition
            # unclassified cytogenetic marker --> chromosome_part
            '7222413': 'SO:0000830',
            '7288448': 'SO:0000341',  # pseudogenic region
            '7288449': 'SO:0001841',  # polymorphic pseudogene
            '7648966': 'SO:0000180',  # retrotransposon
            '7648967': 'SO:0000624',  # telomere
            '7648968': 'SO:0000643',  # minisatellite
            # unclassified other genome feature --> sequence feature
            '7648969': 'SO:0000110',
            # endogenous retroviral region --> endogenous retroviral sequence
            '9272146': 'SO:0000903',
            # mutation defined region --> sequence variant
            '11928467': 'SO:0001060',
            # intronic lncRNA gene  --> ncRNA gene   #FIXME
            '15406204': 'SO:0001263',
            '15406205': 'SO:0000307',  # CpG island
            '15406206': 'SO:0000374',  # ribozyme gene  --> ribozyme
            '15406207': 'SO:0000167',  # promoter
        }
        if marker_type_key.strip() in marker_category_to_so_map:
            so_id = marker_category_to_so_map.get(marker_type_key)
        else:
            logger.error("Marker Category (%s) not mapped", marker_type_key)

        return so_id

    @staticmethod
    def _map_zygosity(zygosity):
        zygtype = None
        type_map = {
            'Heterozygous': 'GENO:0000135',
            'Heteroplasmic': 'GENO:0000603',
            'Homozygous': 'GENO:0000136',
            'Homoplasmic': 'GENO:0000602',
            'Hemizygous Insertion': 'GENO:0000606',
            'Hemizygous Deletion': 'GENO:0000134',  # hemizygous
            'Hemizygous X-linked': 'GENO:0000605',
            'Hemizygous Y-linked': 'GENO:0000604',
            'Indeterminate': 'GENO:0000137'
        }
        if zygosity.strip() in type_map:
            zygtype = type_map.get(zygosity)
        else:
            logger.error("Zygosity (%s) not mapped", zygosity)

        return zygtype

    @staticmethod
    def _map_marker_type(marker_type):
        marktype = None
        type_map = {
            'Complex/Cluster/Region': 'SO:0000110',  # sequence feature
            # transgene
            # WE redefine these as transgenic insertion features instead
            'Transgene': 'SO:0001218',
            'Gene': 'SO:0000704',  # gene
            'QTL': 'SO:0000771',  # QTL
            # sequence_feature. sequence_motif=SO:0001683? region=SO:0000001
            'DNA Segment': 'SO:0000110',
            'Pseudogene': 'SO:0000336',  # pseudogene
            'Cytogenetic Marker': 'SO:0001645',  # genetic_marker?   # fixme
            # sequence_feature. Or sequence_motif=SO:0001683?
            'Other Genome Feature': 'SO:0000110',
            # BAC_end: SO:0000999, YAC_end: SO:00011498; using parent term
            'BAC/YAC end': 'SO:0000150'
        }
        if marker_type.strip() in type_map:
            marktype = type_map.get(marker_type)
        else:
            logger.error("Marker Type (%s) not mapped", marker_type)

        return marktype

    @staticmethod
    def _map_allele_type(allele_type):
        """
        This makes the assumption that all things are variant_loci
        (including Not Specified)
        :param allele_type:
        :return:

        """
        # assume it's a variant locus
        altype = Genotype.genoparts['variant_locus']
        type_map = {
            'Not Applicable': Genotype.genoparts['reference_locus'],
            # should QTLs be something else?  or SO:QTL?
            'QTL': Genotype.genoparts['reference_locus'],
        }
        if allele_type.strip() in type_map:
            altype = type_map.get(allele_type)

        return altype

    @staticmethod
    def _map_taxon(taxon_name):
        taxtype = None
        type_map = {
            'Bos taurus': 'NCBITaxon:9913',
            'Canis familiaris': 'NCBITaxon:9615',
            'Capra hircus': 'NCBITaxon:9925',
            'Cavia porcellus': 'NCBITaxon:10141',
            'Cricetulus griseus': 'NCBITaxon:10029',
            'Danio rerio': 'NCBITaxon:7955',
            'Equus caballus': 'NCBITaxon:9796',
            'Felis catus': 'NCBITaxon:9685',
            'Gallus gallus': 'NCBITaxon:9031',
            'Gorilla gorilla': 'NCBITaxon:9593',
            'Homo sapiens': 'NCBITaxon:9606',
            'Macaca mulatta': 'NCBITaxon:9544',
            'Macropus eugenii': 'NCBITaxon:9315',
            'Mesocricetus auratus': 'NCBITaxon:10036',
            'Microcebus murinus': 'NCBITaxon:30608',
            #  10090=Mus musculus, 10092=Mus musculus domesticus
            'Mus musculus/domesticus': 'NCBITaxon:10090',
            'Ornithorhynchus anatinus': 'NCBITaxon:9258',
            'Oryctolagus cuniculus': 'NCBITaxon:9986',
            'Ovis aries': 'NCBITaxon:9940',
            'Pan troglodytes': 'NCBITaxon:9598',
            'Pongo pygmaeus': 'NCBITaxon:9600',
            'Rattus norvegicus': 'NCBITaxon:10116',
            # 9823=Sus scrofa, 9825=Sus scrofa domestica
            'Sus scrofa domestica L.': 'NCBITaxon:9823',
            'Xenopus (Silurana) tropicalis': 'NCBITaxon:8364',
        }
        if taxon_name.strip() in type_map:
            taxtype = type_map.get(taxon_name)
        else:
            logger.error("Taxon Name (%s) not mapped", taxon_name)

        return taxtype

    @staticmethod
    def _map_strain_species(species):
        # make the assumption that it is a Mus genus, unless if specified
        tax = '10088'
        id_map = {
            'laboratory mouse and M. m. domesticus (brevirostris)': '116058',
            'laboratory mouse and M. m. bactrianus': '35531',
            # does not include lab mouse
            'laboratory mouse and M. m. castaneus and M. m. musculus': '477816',
            'M. m. domesticus and M. m. molossinus and M. m. castaneus': '1266728',
            'M. setulosus': '10102',
            'laboratory mouse and wild-derived': '10088',  # Mus genus
            'laboratory mouse and M. m. musculus (Prague)': '39442',
            'M. m. castaneus and M. m. musculus': '477816',
            'M. m. domesticus (Canada)': '10092',
            # FIXME
            'M. m. domesticus and M. m. domesticus poschiavinus': '10092',
            'M. m. musculus and M. spretus or M. m. domesticus': '186842',
            'M. chypre': '862507',  # unclassified mus
            'M. cookii (Southeast Asia)': '10098',
            'M. cervicolor': '10097',
            'M. m. gentilulus': '80274',
            'M. m. domesticus poschiavinus (Tirano, Italy)': '10092',  # FIXME
            'M. m. castaneus (Phillipines)': '10091',
            'M. m. domesticus (brevirostris) (France)': '116058',
            'M. m. domesticus (Lake Casitas, CA)': '10092',
            'laboratory mouse or wild': '862507',
            'M. spretus (Morocco)': '10096',
            'M. m. bactrianus (Russia)': '35531',
            'M. m. musculus (Pakistan)': '39442',
            'M. tenellus': '397330',
            'M. baoulei': '544437',
            'laboratory mouse and M. spretus or M. m. castaneus': '862507',
            'M. haussa': '273922',
            'M. m. domesticus (Montpellier, France)': '10092',
            'M. m. musculus': '39442',
            'M. m. domesticus and Not specified': '10090',
            'laboratory mouse and M. spretus or M. m. musculus': '862507',
            'Rat': '10114',
            'M. m. musculus (Czech)': '39442',
            'M. spretus or laboratory mouse and M. m. musculus': '862507',
            'M. m. domesticus (Skokholm, UK)': '10092',
            'M. shanghai': '862507',
            'M. m. musculus (Korea)': '39442',
            'laboratory mouse and M. m. domesticus (Peru)': '39442',
            'M. m. molossinus and laboratory mouse': '57486',
            'laboratory mouse and M. m. praetextus x M. m. domesticus': '10088',
            'M. cervicolor popaeus': '135828',
            'M. m. domesticus or M. m. musculus': '10090',
            'M. abbotti': '10108',
            'M. nitidulus': '390848',
            'M. gradsko': '10088',
            'M. m. domesticus poschiavinus (Zalende)': '10092',
            'M. m. domesticus (MD)': '10092',
            'M. triton': '473865',
            'M. m. domesticus (Tunisia)': '10092',
            'Wild and outbred Swiss and M. m. domesticus and M. m. molossinus': '862507',
            'M. m. molossinus and M. spretus': '186842',
            'M. spretus (France)': '10096',
            'M. m. praetextus x M. m. domesticus and M. m. domesticus': '10090',
            'M. m. bactrianus': '35531',
            'M. m. domesticus (Bulgaria)': '10092',
            'Not Specified and M. m. domesticus': '10090',
            'laboratory mouse and M. m. molossinus': '10090',
            'M. spretus or M. m. domesticus': '862507',
            'M. m. molossinus (Anjo, Japan)': '57486',
            'M. m. castaneus or M. m. musculus': '10090',
            'laboratory mouse and M. m. castaneus and M. m. domesticus poschiavinus': '10090',
            'M. spretus (Spain)': '10096',
            'M. m. molossinus (Japan)': '57486',
            'M. lepidoides': '390847',
            'M. m. macedonicus': '10090',
            'M. m. musculus (Pohnpei-Micronesia)': '39442',
            'M. m. domesticus (Morocco)': '10092',
            'laboratory mouse and M. m. domesticus poschiavinus': '10090',
            'M. (Nannomys) (Zakouma NP, Chad)': '862510',
            'laboratory mouse and M. m. musculus (Tubingen, S Germany)': '10090',
            'M. cypriacus': '468371',
            'laboratory mouse and M. m. domesticus': '10090',
            'M. m. musculus and M. m. domesticus': '477815',
            'M. crociduroides': '41269',
            'laboratory mouse and M. m. domesticus (Israel)': '10090',
            'M. m. domesticus (Ann Arbor, MI)': '10092',
            'M. caroli': '10089',
            'M. emesi': '887131',
            'M. gratus': '229288',
            'laboratory mouse or M. spretus': '10090',
            'M. m. domesticus (PA)': '10092',
            'M. m. domesticus (DE)': '10092',
            'laboratory mouse and M. m. castaneus': '10090',
            'M. hortulanus (Austria)': '10103',  # synonymous with Mus spicilegus
            'M. m. musculus and M. hortulanus': '862507',
            'Not Applicable': None,
            'M. mattheyi': '41270',
            'M. m. macedonicus or M. m. musculus': '10090',
            'M. m. musculus (Georgia)': '39442',
            'M. famulus': '83773',
            'M. m. bactrianus (Iran)': '35531',
            'M. m. musculus (Toshevo, Bulgaria)': '39442',
            'M. m. musculus (Prague)': '39442',
            'M. platythrix (India)': '10101',
            'M. m. domesticus': '10092',
            'M. m. castaneus and laboratory mouse': '10090',
            'M. m. domesticus (Australia)': '10092',
            'laboratory mouse and M. m. domesticus and M. m. molossinus': '10090',
            'M. m. domesticus (brevirostris) (Morocco)': '10092',
            'laboratory mouse and M. spretus and M. m. musculus': '10090',
            'M. m. molossinus': '57486',
            'M. hortulanus': '10103',  # synonymous with Mus spicilegus
            'M. dunni': '254704',  # synonymous with Mus terricolor
            'M. m. castaneus (Pathumthani, Thailand)': '10091',
            'M. terricolor': '254704',
            'M. m. domesticus (China)': '10092',
            'M. fragilicauda': '186193',
            'Not Resolved': '10088',  # assume Mus
            'M. m. musculus or M. spretus and laboratory mouse': '10088',
            'M. macedonicus macedonicus': '270353',
            'laboratory mouse': '10090',
            'laboratory mouse and M. abbotti': '10088',
            'M. cervicolor cervicolor': '135827',
            'laboratory mouse and M. spretus': '10088',
            'M. m. domesticus (praetextus)': '10092',
            'M. spretus (London)': '10096',
            'M. m. musculus (Prague) and laboratory mouse': '10090',
            'laboratory mouse and M. m. musculus (Central Jutland, Denmark)': '10090',
            'M. m. musculus (Northern Jutland, Denmark)': '39442',
            'M. m. castaneus (Taiwan)': '10091',
            'M. brevirostris': '116058',
            'M. spretus': '10096',
            'M. m. domesticus poschiavinus': '10092',
            'M. pahari': '10093',
            'M. m. molossinus (Hakozaki, Japan)': '57486',
            'M. m. musculus (Hokkaido)': '39442',
            'M. m. musculus or M. spretus': '862507',
            'laboratory mouse and M. m. musculus': '10090',
            'Not Specified and M. m. musculus': '10090',
            'M. m. castaneus (Masinagudi, South India)': '10091',
            'M. indutus': '273921',
            'M. saxicola': '10094',
            'M. m. praetextus x M. m. musculus': '10090',
            'M. spretus and laboratory mouse': '862507',
            'laboratory mouse and M. m. castaneus and M. m. domesticus': '10090',
            'M. m. musculus (Bulgaria)': '39442',
            'M. booduga': '27681',
            'Not Specified': '10090',  # OK?
            'Not Specified and M. m. molossinus': '10090',
            'M. m. musculus and M. spretus': '862507',
            'M. minutoides': '10105',
            'M. spretus (Tunisia)': '10096',
            'M. spicilegus': '10103',
            'Peru Coppock': '10088',  # unclassified mus
            'M. m. domesticus (CA)': '10092',
            'M. macedonicus spretoides': '270352',
            'M. m. musculus and M. m. castaneus': '477816',
            'M. m. praetextus x M. m. domesticus': '10090',
            'M. m. domesticus and M. spretus': '862507',
            'M. m. molossinus (Mishima)': '57486',
            'M. hortulanus and M. m. macedonicus': '10088',
            'Not Specified and M. spretus': '10088',
            'M. m. domesticus (Peru)': '10092',
            'M. m. domesticus (OH)': '10092',
            'M. m. bactrianus and laboratory mouse': '10090',
            'laboratory mouse and M. m. domesticus (OH)': '10090',
            'M. m. gansuensis': '1385377',
            'laboratory mouse and M. m. castaneus and M. spretus': '10088',
            'M. m. castaneus': '10091',
            'Wild and M. m. domesticus': '10088',
            'M. m. molossinus (Nagoya)': '57486'
        }
        if species.strip() in id_map:
            tax = id_map.get(species.strip())
        else:
            logger.warning("Species (%s) not mapped; defaulting to Mus genus.",
                           species)

        if tax is not None:
            tax = 'NCBITaxon:'+tax
        return tax

    @staticmethod
    def _map_evidence_id(evidence_code):

        ecotype = None
        type_map = {
            'EXP': 'ECO:0000006',
            'IBA': 'ECO:0000318',
            'IC': 'ECO:0000001',
            'IDA': 'ECO:0000314',
            'IEA': 'ECO:0000501',
            'IEP': 'ECO:0000008',
            'IGI': 'ECO:0000316',
            'IKR': 'ECO:0000320',
            'IMP': 'ECO:0000315',
            'IPI': 'ECO:0000353',
            'ISA': 'ECO:0000200',
            'ISM': 'ECO:0000202',
            'ISO': 'ECO:0000201',
            'ISS': 'ECO:0000250',
            'NAS': 'ECO:0000303',
            'ND': 'ECO:0000035',
            'RCA': 'ECO:0000245',
            'TAS': 'ECO:0000304'
        }
        if evidence_code.strip() in type_map:
            ecotype = type_map.get(evidence_code)
        else:
            logger.error("Evidence code (%s) not mapped", evidence_code)

        return ecotype

    @staticmethod
    def _makeInternalIdentifier(prefix, key):
        """
        This is a special MGI-to-MONARCH-ism.
        MGI tables have unique keys that we use here, but don't want to
        necessarily re-distribute those internal identifiers.
        Therefore, we make them into keys in a consistent way here.
        :param prefix: the object type to prefix the key with,
        since the numbers themselves are not unique across tables
        :param key: the number
        :return:

        """
        iid = '_:mgi'+prefix+'key'+key

        return iid

    # def _querysparql(self):
    #
    #     #load the graph
    #     vg = Graph()
    #     vg.parse(self.outfile, format="turtle")
    #
    #     qres = g.query(
    #         """SELECT DISTINCT ?aname ?bname
    #             WHERE {
    #                 ?a foaf:knows ?b .
    #                 ?a foaf:name ?aname .
    #                 ?b foaf:name ?bname .
    #             }""")
    #
    #     for row in qres:
    #         print("%s knows %s" % row)
    #
    #     return

    def getTestSuite(self):
        import unittest
        from tests.test_mgi import MGITestCase
        # TODO test genotypes

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MGITestCase)

        return test_suite
