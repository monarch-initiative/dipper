import csv
import os
from datetime import datetime
import logging
import re
from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.assoc.Association import Assoc
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper import config
from dipper.models.GenomicFeature import Feature, makeChromID

LOG = logging.getLogger(__name__)


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
    details in your conf.yaml file, like:
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

    resources = {
        'query_map': [
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
                'query': '../../resources/sql/mgi/evidence.sql',
                'outfile': 'evidence_view'
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
        ],
        'test_keys': '../../resources/mgi_test_keys.yaml'
    }

    # for testing purposes, this is a list of internal db keys
    # to match and select only portions of the source

    def __init__(
            self,
            graph_type,
            are_bnodes_skolemized
    ):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            name='mgi',
            ingest_title='Mouse Genome Informatics',
            ingest_url='http://www.informatics.jax.org/',
            license_url='http://www.informatics.jax.org/mgihome/other/copyright.shtml',
            data_rights=None,
            file_handle=None)

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

        # also add the gene ids from the test_ids
        # in order to capture transgenes of the test set

        if 'gene' in self.all_test_ids:
            self.test_ids = self.all_test_ids['gene']
        else:
            LOG.warning("not configured with gene test ids.")
            self.test_ids = []

        self.test_keys = self.open_and_parse_yaml(self.resources['test_keys'])

        return

    def fetch(self, is_dl_forced=False):
        """
        For the MGI resource, we connect to the remote database,
        and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        """
        # check if config exists; if it doesn't, error out and let user know
        if 'dbauth' not in config.get_config() and 'mgi' \
                not in config.get_config()['dbauth']:
            LOG.error("not configured with PG user/password.")

        # create the connection details for MGI
        cxn = config.get_config()['dbauth']['mgi']

        self.dataset.setFileAccessUrl(''.join((
            'jdbc:postgresql://', cxn['host'], ':', str(cxn['port']), '/',
            cxn['database'])), is_object_literal=True)

        # process the tables
        # self.fetch_from_pgdb(self.tables, cxn, 100)  # for testing only
        # self.fetch_from_pgdb(self.tables, cxn, None, is_dl_forced)

        for query_map in self.resources['query_map']:
            query_fh = open(os.path.join(
                os.path.dirname(__file__), query_map['query']), 'r')
            query = query_fh.read()
            force = False
            if 'Force' in query_map:
                force = query_map['Force']
            self.fetch_query_from_pgdb(
                query_map['outfile'], query, None, cxn)
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
                dat = cols[1].strip().split('.')[0]
                datestamp = datetime.strptime(
                    dat, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
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
            LOG.info("Only parsing first %d rows of each file", limit)
        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

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
        self._process_evidence_view(limit)
        self._process_mgi_note_vocevidence_view(limit)
        self._process_mrk_location_cache(limit)
        self.process_mgi_relationship_transgene_genes(limit)
        self.process_mgi_note_allele_view(limit)

        LOG.info("Finished parsing.")

        LOG.info("Loaded %d nodes", len(self.graph))
        return

    def fetch_transgene_genes_from_db(self, cxn):
        """
        This is a custom query to fetch the non-mouse genes that
        are part of transgene alleles.

        :param cxn:
        :return:
        """

        query = '''
SELECT  r._relationship_key as rel_key,
        r._object_key_1 as object_1,
        a.accid as allele_id,
        alabel.label as allele_label,
        rc._category_key as category_key,
        rc.name as category_name,
        t._term_key as property_key,
        t.term as property_name,
        rp.value as property_value
    FROM mgi_relationship r
    JOIN mgi_relationship_category rc ON r._category_key = rc._category_key
    JOIN acc_accession a  ON r._object_key_1 = a._object_key
        AND rc._mgitype_key_1 = a._mgitype_key
        AND a._logicaldb_key = 1
    JOIN all_label alabel ON a._object_key = alabel._allele_key
        AND alabel._label_status_key = 1
        AND alabel.priority = 1
    JOIN mgi_relationship_property rp ON r._relationship_key = rp._relationship_key
        AND rp._propertyname_key = 12948292
        JOIN voc_term t ON rp._propertyname_key = t._term_key
    WHERE r._category_key = 1004
        '''

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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        geno = Genotype(graph)
        model = Model(graph)

        raw = '/'.join((self.rawdir, 'gxd_genotype_view'))
        LOG.info("getting genotypes and their backgrounds")
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line = line.rstrip("\n")
                line_counter += 1
                (genotype_key, strain_key, strain, mgiid) = line.split('\t')

                if self.test_mode is True:
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
                background_type = self.globaltt['genomic_background']
                if strain_id is None or int(strain_key) < 0:
                    if strain_id is None:
                        # some of the strains don't have public identifiers!
                        # so we make one up, and add it to the hash
                        strain_id = self._makeInternalIdentifier('strain', strain_key)
                        self.idhash['strain'].update({strain_key: strain_id})
                        model.addComment(strain_id, "strain_key:" + strain_key,
                                         subject_category=
                                         blv.PopulationOfIndividualOrganisms.value)
                    elif int(strain_key) < 0:
                        # these are ones that are unidentified/unknown.
                        # so add instances of each.
                        strain_id = self._makeInternalIdentifier(
                            'strain', re.sub(r':', '', str(strain_id)))
                        strain_id += re.sub(r':', '', str(mgiid))
                        strain_id = re.sub(r'^_', '_:', strain_id)
                        strain_id = re.sub(r'::', ':', strain_id)
                        model.addDescription(
                            strain_id,
                            "This genomic background is unknown.  " +
                            "This is a placeholder background for " +
                            mgiid + ".",
                            subject_category=
                            blv.PopulationOfIndividualOrganisms.value)
                        background_type = self.globaltt[
                            'unspecified_genomic_background']

                    # add it back to the idhash
                    LOG.info(
                        "adding background as internal id: %s %s: %s",
                        strain_key, strain, strain_id)

                geno.addGenomicBackgroundToGenotype(
                    strain_id, mgiid, background_type)

                self.label_hash[strain_id] = strain

                # add BG to a hash so we can build the genotype label later
                self.geno_bkgd[mgiid] = strain_id

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        model = Model(graph)
        line_counter = 0
        geno_hash = {}
        raw = '/'.join((self.rawdir, 'gxd_genotype_summary_view'))
        LOG.info("building labels for genotypes")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (object_key, preferred, mgiid, subtype,
                 short_description) = line.split('\t')

                if self.test_mode is True:
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

                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        # now, loop through the hash and add the genotypes as individuals
        # we add the mgi genotype as a synonym
        # (we generate our own label later)
        geno = Genotype(graph)
        for gt in geno_hash:
            genotype = geno_hash.get(gt)
            gvc = sorted(genotype.get('vslcs'))
            label = '; '.join(gvc) + ' [' + genotype.get('subtype') + ']'
            geno.addGenotype(gt, None)
            model.addComment(gt, self._makeInternalIdentifier(
                'genotype', genotype.get('key')), subject_category=blv.Genotype.value)
            model.addSynonym(gt, label.strip(), subject_category=blv.Genotype.value)

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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'all_summary_view'))
        LOG.info(
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
                    LOG.warning('Expected ' + str(col_count) + ' columns.')
                    LOG.warning('Received ' + str(cols) + ' columns.')
                    LOG.warning(line.format())
                    continue
                # no stray tab in the description column
                (object_key, preferred, mgiid, description,
                 short_description) = line.split('\t')
                # NOTE: May want to filter alleles based on the preferred field
                # (preferred = 1) or will get duplicates
                # (24288, to be exact...
                # Reduced to 480 if filtered on preferred = 1)

                if self.test_mode is True:
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
                        altype, description.strip(),
                        blv.SequenceVariant.value)
                    self.label_hash[mgiid] = short_description.strip()

                # TODO deal with non-preferreds, are these deprecated?

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        line_counter = 0
        LOG.info(
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
                    LOG.warning('Expected ' + str(col_count) + ' columns.')
                    LOG.warning('Received ' + str(cols) + ' columns.')
                    LOG.warning(line.format())
                    continue

                (allele_key, marker_key, strain_key, symbol,
                 name, iswildtype) = line.split('\t')

                # TODO update processing to use this view better
                # including jnums!

                if self.test_mode is True:
                    if int(allele_key) not in self.test_keys.get('allele'):
                        continue

                allele_id = self.idhash['allele'].get(allele_key)
                if allele_id is None:
                    LOG.error(
                        "what to do! can't find allele_id. skipping %s %s",
                        allele_key, symbol)
                    continue

                marker_id = None

                if marker_key is not None and marker_key != '':
                    # we make the assumption here that the markers
                    # have already been added to the table
                    marker_id = self.idhash['marker'].get(marker_key)
                    if marker_id is None:
                        LOG.error(
                            "what to do! can't find marker_id. skipping %s %s",
                            marker_key, symbol)
                        continue

                iseqalt_id = self._makeInternalIdentifier('seqalt', allele_key)

                # for non-wild type alleles:
                if iswildtype == '0':
                    locus_type = self.globaltt['variant_locus']
                    locus_rel = self.globaltt['is_allele_of']
                # for wild type alleles:
                elif iswildtype == '1':
                    locus_type = self.globaltt['reference_locus']
                    locus_rel = self.globaltt['is_reference_allele_of']
                    # add the allele to the wildtype set for lookup later
                    self.wildtype_alleles.add(allele_id)
                else:
                    locus_rel = None
                    locus_type = None

                model.addIndividualToGraph(allele_id, symbol, locus_type,
                                           blv.SequenceVariant.value)
                model.makeLeader(allele_id, blv.SequenceVariant.value)
                self.label_hash[allele_id] = symbol
                self.idhash['seqalt'][allele_key] = iseqalt_id

                # HACK - if the label of the allele == marker,
                # then make the thing a seq alt
                allele_label = self.label_hash.get(allele_id)
                marker_label = self.label_hash.get(marker_id)
                if allele_label is not None and allele_label == marker_label:
                    model.addSameIndividual(allele_id, marker_id,
                                            subject_category=
                                            blv.SequenceVariant.value,
                                            object_category=
                                            blv.InformationContentEntity.value)
                    self.idhash['seqalt'][allele_key] = allele_id
                    model.addComment(
                        allele_id,
                        self._makeInternalIdentifier('allele', allele_key),
                        blv.SequenceVariant.value)
                elif marker_id is not None:
                    # marker_id will be none if the allele
                    # is not linked to a marker
                    # (as in, it's not mapped to a locus)
                    geno.addAlleleOfGene(allele_id, marker_id, locus_rel)

                # sequence alteration in strain

                if iswildtype == '0':
                    sa_label = symbol
                    sa_id = iseqalt_id

                    if marker_key is not None \
                            and allele_label != marker_label and marker_key != '':
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

                    geno.addSequenceAlteration(sa_id, sa_label, None, name)
                    self.label_hash[sa_id] = sa_label

                    strain_id = self.idhash['strain'].get(strain_key)
                    # scrub out if the strain is "not specified"
                    if strain_id is not None and \
                            strain_id not in ['MGI:4867032', 'MGI:5649511']:
                        geno.addSequenceDerivesFrom(allele_id, strain_id,
                                                    subject_category=
                                                    blv.SequenceVariant.value,
                                                    object_category=
                                                    blv.PopulationOfIndividualOrganisms.value)

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'gxd_allelepair_view'))
        LOG.info("processing allele pairs (VSLCs) for genotypes")
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

                if self.test_mode is True:
                    if int(genotype_key) not in self.test_keys.get('genotype'):
                        continue

                genotype_id = self.idhash['genotype'].get(genotype_key)
                if genotype_id not in geno_hash:
                    geno_hash[genotype_id] = set()
                if genotype_id is None:
                    LOG.error(
                        "genotype_id not found for key %s; skipping", genotype_key)
                    continue

                allele1_id = self.idhash['allele'].get(allele_key_1)
                allele2_id = self.idhash['allele'].get(allele_key_2)

                # Need to map the allelestate to a zygosity term
                zygosity_id = self.resolve(allelestate.strip())
                ivslc_id = self._makeInternalIdentifier('vslc', allelepair_key)

                geno_hash[genotype_id].add(ivslc_id)
                # TODO: VSLC label likely needs processing similar to
                # the processing in the all_allele_view
                # FIXME: handle null alleles
                vslc_label = allele1+'/'
                if allele2_id is None:
                    if zygosity_id in [
                            self.globaltt['hemizygous insertion-linked'],
                            self.globaltt['hemizygous-x'],
                            self.globaltt['hemizygous-y'],
                            self.globaltt['hemizygous']]:
                        vslc_label += '0'
                    elif zygosity_id == self.globaltt['heterozygous']:
                        vslc_label += '+'
                    elif zygosity_id == self.globaltt['indeterminate']:
                        vslc_label += '?'
                    elif zygosity_id == self.globaltt['homozygous']:
                        # we shouldn't get here, but for testing this is handy
                        vslc_label += allele1
                    else:  # heteroplasmic,  homoplasmic,  FIXME add these if possible
                        LOG.info(
                            "A different kind of zygosity found is: %s",
                            self.globaltcid[zygosity_id])
                        vslc_label += '?'
                else:
                    vslc_label += allele2

                model.addIndividualToGraph(
                    ivslc_id, vslc_label,
                    self.globaltt['variant single locus complement'],
                    ind_category=blv.SequenceVariant.value)
                self.label_hash[ivslc_id] = vslc_label
                rel1 = rel2 = self.globaltt['has_variant_part']
                if allele1_id in self.wildtype_alleles:
                    rel1 = self.globaltt['has_reference_part']
                if allele2_id in self.wildtype_alleles:
                    rel2 = self.globaltt['has_reference_part']
                geno.addPartsToVSLC(
                    ivslc_id, allele1_id, allele2_id, zygosity_id, rel1, rel2)

                # if genotype_id not in geno_hash:
                #     geno_hash[genotype_id] = [vslc_label]
                # else:
                #     geno_hash[genotype_id] += [vslc_label]

                if not self.test_mode and limit is not None and line_counter > limit:
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
                    gvc_id, gvc_label, self.globaltt['genomic_variation_complement'],
                    ind_category=blv.SequenceVariant.value)
                self.label_hash[gvc_id] = gvc_label
                for v in vslcs:
                    geno.addParts(v, gvc_id, self.globaltt['has_variant_part'],
                                  part_category=blv.SequenceVariant.value,
                                  parent_category=blv.Genotype.value)
                    geno.addVSLCtoParent(v, gvc_id)
                geno.addParts(gvc_id, gt, self.globaltt['has_variant_part'],
                              part_category=blv.SequenceVariant.value,
                              parent_category=blv.Genotype.value)
            elif len(vslcs) == 1:
                gvc_id = vslcs[0]
                gvc_label = self.label_hash[gvc_id]
                # type the VSLC as also a GVC
                model.addIndividualToGraph(
                    gvc_id, gvc_label, self.globaltt['genomic_variation_complement'],
                    ind_category=blv.SequenceVariant.value)
                geno.addVSLCtoParent(gvc_id, gt, parent_category=blv.Genotype.value)
            else:
                LOG.info("No VSLCs for %s", gt)

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

            model.addIndividualToGraph(gt, genotype_label,
                                       ind_category=blv.Genotype.value)
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'all_allele_mutation_view'))
        LOG.info("getting mutation types for sequence alterations")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (allele_key, mutation) = line.split('\t')
                iseqalt_id = self.idhash['seqalt'].get(allele_key)
                if iseqalt_id is None:
                    iseqalt_id = self._makeInternalIdentifier('seqalt', allele_key)

                if self.test_mode and int(allele_key) \
                        not in self.test_keys.get('allele'):
                    continue

                # TODO we might need to map the seq alteration to the MGI id
                # for unlocated things; need to use hashmap
                # map the sequence_alteration_type

                seq_alt_type_id = self.resolve(mutation, False)
                if seq_alt_type_id == mutation:
                    LOG.error("No mappjng found for seq alt '%s'", mutation)
                    LOG.info("Defaulting to 'sequence_alteration'")
                    seq_alt_type_id = self.globaltt['sequence_alteration']

                # HACK - if the seq alteration is a transgene,
                # then make sure it is a transgenic insertion
                allele_id = self.idhash['allele'].get(allele_key)
                if allele_id is not None:
                    allele_label = self.label_hash.get(allele_id)
                    if allele_label is not None and re.search(r'Tg\(', allele_label):
                        LOG.info(
                            "Found a transgenic insertion for %s", allele_label)
                        # transgenic_insertion, instead of plain old insertion
                        seq_alt_type_id = self.globaltt["transgenic_insertion"]

                model.addIndividualToGraph(iseqalt_id, None, seq_alt_type_id,
                                           ind_category=blv.SequenceVariant.value)

                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        return

    def _process_voc_annot_view(self, limit):
        """
        This MGI table represents associations between things.

        We add the internal annotation id to the idhashmap.
        It is expected that the genotypes have already been added to the idhash

        :param limit:
        :return:

        """

        # TODO also get Strain/Attributes (annottypekey = 1000)
        # TODO what is Phenotype (Derived) vs
        # non-derived?  (annottypekey = 1015)
        # TODO is evidence in this table?  what is the evidence vocab key?
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        line_counter = 0
        LOG.info("getting G2P associations")
        raw = '/'.join((self.rawdir, 'voc_annot_view'))
        col = [
            'annot_key', 'annot_type', 'object_key', 'term_key', 'qualifier_key',
            'qualifier', 'term', 'accid']
        with open(raw, 'r') as f:
            header = f.readline()  # read the header row; skip
            if header != col:
                LOG.error("\nExpected header: %s\nReceived header: %s", col, header)

            for line in f:
                row = line.rstrip('\n').split('\t')

                annot_key = row[col.index('annot_key')]
                annot_type = row[col.index('annot_type')]
                object_key = row[col.index('object_key')]
                term_key = row[col.index('term_key')]
                qualifier_key = row[col.index('qualifier_key')]
                # qualifier,
                # term,
                accid = row[col.index('accid')]

                if self.test_mode is True:
                    if int(annot_key) not in self.test_keys.get('annot'):
                        continue

                # iassoc_id = self._makeInternalIdentifier('annot', annot_key)
                # assoc_id = self.make_id(iassoc_id)

                assoc_id = None
                # Mammalian Phenotype/Genotype are curated G2P assoc
                if annot_type == 'Mammalian Phenotype/Genotype':
                    line_counter += 1

                    # We expect the label for the phenotype
                    # to be taken care of elsewhere
                    model.addClassToGraph(accid, None,
                                          class_category=blv.PhenotypicFeature.value)

                    genotype_id = self.idhash['genotype'].get(object_key)
                    if genotype_id is None:
                        LOG.error(
                            "can't find genotype id for %s", object_key)
                    else:
                        # add the association
                        assoc = G2PAssoc(graph, self.name, genotype_id, accid)
                        assoc.add_association_to_graph()
                        assoc_id = assoc.get_association_id()
                # OMIM/Genotype are disease-models
                elif annot_type == 'DO/Genotype':
                    # skip NOT annotations for now FIXME
                    if qualifier_key == '1614157':
                        continue
                    genotype_id = self.idhash['genotype'].get(object_key)
                    if genotype_id is None:
                        LOG.error("can't find genotype id for %s", object_key)
                    else:
                        # add the association
                        assoc = Assoc(graph, self.name)
                        # TODO PYLINT
                        # Redefinition of assoc type from
                        # dipper.models.assoc.G2PAssoc.G2PAssoc to
                        # dipper.models.assoc.Association.Assoc
                        assoc.set_subject(genotype_id)
                        assoc.set_object(accid)
                        assoc.set_relationship(self.globaltt['is model of'])
                        assoc.add_association_to_graph()
                        assoc_id = assoc.get_association_id()
                elif annot_type == 'MCV/Marker':
                    # marker category == type
                    marker_id = self.idhash['marker'].get(object_key)
                    if str(term_key).strip() in self.localtt:
                        # check "Not Applicable": "reference_locus"
                        term_id = self.resolve(str(term_key).strip())
                    else:
                        term_id = None
                        logging.warning('No type mapping for: %s', term_key)
                    # note that the accid here is an internal mouse cv term,
                    # and we don't use it.
                    if term_id is not None and marker_id is not None:
                        # do something special for transgenics -
                        # make sure these are transgenic insertions
                        model.addType(marker_id, term_id)
                elif annot_type == 'DO/Allele':  # allele/Disease
                    allele_id = self.idhash['allele'].get(object_key)
                    if allele_id is None:
                        LOG.error("can't find genotype id for %s", object_key)
                    else:
                        # add the association
                        assoc = Assoc(graph, self.name)
                        assoc.set_subject(allele_id)
                        assoc.set_object(accid)
                        assoc.set_relationship(self.globaltt['is model of'])
                        assoc.add_association_to_graph()
                        assoc_id = assoc.get_association_id()

                if assoc_id is not None:
                    # add the assoc to the hashmap (using the monarch id)
                    self.idhash['annot'][annot_key] = assoc_id
                    model.addComment(assoc_id, "annot_key:" + annot_key)

                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        return

    def _process_evidence_view(self, limit):
        """
        Here we fetch the evidence (code and publication) for the associations.
        The evidence codes are mapped from the standard GO codes to ECO.
        J numbers are added for publications.
        We will only add the evidence if the annotation is in our idhash.

        We also pull in evidence qualifiers, as of June 2018 they are
        Data Interpretation Center (eg IMPC)
        external ref (eg UniProtKB:Q9JHI2-3 for Proteoform/Marker assoc)
        Phenotyping Center (eg WTSI)
        Resource Name (eg MGP)
        MP-Sex-Specificity (eg NA, M, F)

        Triples:
        <annot_id> dc:evidence <evidence_id>
        <pub_id> a owl:NamedIndividual
        <annot_id> dc:source <pub_id>

        :param limit:
        :return:

        """

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        model = Model(graph)
        line_counter = 0
        LOG.info("getting evidence and pubs for annotations")
        raw = '/'.join((self.rawdir, 'evidence_view'))
        col = [
            'annot_evidence_key', 'annot_key', 'evidence_code', 'jnumid', 'qualifier',
            'qualifier_value', 'annotation_type']
        with open(raw, 'r') as reader:
            reader.readline()  # read the header row; skip
            for line in reader:
                line = line.rstrip("\n")
                line_counter += 1
                row = line.split('\t')

                annot_evidence_key = row[col.index('annot_evidence_key')]
                annot_key = row[col.index('annot_key')]
                evidence_code = row[col.index('evidence_code')]
                jnumid = row[col.index('jnumid')]
                qualifier = row[col.index('qualifier')]
                qualifier_value = row[col.index('qualifier_value')]
                # annotation_type = row[col.index('annotation_type')]

                if self.test_mode and annot_key not in self.test_keys.get('annot'):
                    continue

                # add the association id to map to the evidence key
                # (to attach the right note to the right assn)
                self.idhash['notes'][annot_evidence_key] = annot_key

                assoc_id = self.idhash['annot'].get(annot_key)

                if assoc_id is None:
                    # assume that we only want to add the evidence/source
                    # for annots that we have in our db
                    continue

                evidence_id = self.resolve(evidence_code)

                reference = Reference(graph, jnumid)
                reference.addRefToGraph()

                # add the ECO and citation information to the annot
                model.addTriple(assoc_id, self.globaltt['has evidence'], evidence_id)
                model.addTriple(assoc_id, self.globaltt['source'], jnumid)

                # For Mammalian Phenotype/Genotype annotation types
                # MGI adds sex specificity qualifiers here
                if qualifier == 'MP-Sex-Specificity' and qualifier_value in ('M', 'F'):
                    model._addSexSpecificity(assoc_id, self.resolve(qualifier_value))

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        # firstpass, get the J number mapping, and add to the global hash

        LOG.info('populating pub id hash')
        raw = '/'.join((self.rawdir, 'bib_acc_view'))
        col = [
            'accid', 'prefixpart', 'numericpart', 'object_key', 'logical_db',
            'logicaldb_key']
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            header = next(filereader)
            if header != col:
                LOG.error('bib_acc_view expected:\n%s\n\tBut got:\n%s', col, header)
            for row in filereader:

                accid = row[col.index('accid')]
                prefixpart = row[col.index('prefixpart')]
                # 'numericpart'
                object_key = int(row[col.index('object_key')])
                # logical_db = row[col.index('logical_db')]
                # logicaldb_key = row[col.index('logicaldb_key')]

                if self.test_mode and object_key not in self.test_keys.get('pub'):
                    continue

                # we use the J number here because
                # it is the externally-accessible identifier
                if prefixpart != 'J:':
                    continue
                self.idhash['publication'][object_key] = accid
                reference = Reference(graph, accid)
                reference.addRefToGraph()

                if not self.test_mode and limit is not None and \
                        filereader.line_num > limit:
                    break

        # 2nd pass, look up the MGI identifier in the hash
        LOG.info("getting pub equivalent ids")
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            header = next(filereader)

            for row in filereader:
                accid = row[col.index('accid')]
                prefixpart = row[col.index('prefixpart')]
                # 'numericpart'
                object_key = int(row[col.index('object_key')])
                logical_db = row[col.index('logical_db')]
                logicaldb_key = row[col.index('logicaldb_key')]

                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('pub'):
                        continue
                logical_db = logical_db.strip()
                jid = self.idhash['publication'].get(object_key)
                pub_id = None
                if logicaldb_key == '29':  # pubmed
                    pub_id = 'PMID:' + accid
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
                    reference = Reference(graph, pub_id)

                    # make the assumption that if it is a PMID, it is a journal
                    if re.match(r'PMID', pub_id):
                        reference.setType(self.globaltt['journal article'])
                        model.makeLeader(pub_id)
                    reference.addRefToGraph()

                    model.addSameIndividual(jid, pub_id)
                else:
                    LOG.warning(
                        "Publication from (%s) not mapped for %s",
                        logical_db, object_key)

                if not self.test_mode and limit is not None and \
                        filereader.line_num > limit:
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

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        line_counter = 0
        geno = Genotype(graph)
        raw = '/'.join((self.rawdir, 'prb_strain_view'))
        LOG.info("getting strains and adding their taxa")
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue
                (strain_key, strain, species) = line

                if self.test_mode is True:
                    if int(strain_key) not in self.test_keys.get('strain'):
                        continue

                strain_id = self.idhash['strain'].get(strain_key)

                if strain_id is not None:
                    self.label_hash[strain_id] = strain

                    # add the species to the graph as a class
                    species = species.strip()
                    sp = self.resolve(species, False)
                    if sp == species:
                        LOG.error("No taxon mapping for " + species)
                        LOG.warning("defaulting to Mus Genus")
                        sp = self.globaltt['Mus']

                    model.addClassToGraph(sp, None)
                    geno.addTaxon(sp, strain_id)
                    model.addIndividualToGraph(strain_id, strain, sp)

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        geno = Genotype(graph)
        line_counter = 0
        raw = '/'.join((self.rawdir, 'mrk_marker_view'))
        LOG.info("getting markers and assigning types")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1

                (marker_key, organism_key, marker_status_key,
                 symbol, name, latin_name, marker_type) = line.split('\t')

                if self.test_mode is True:
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
                        LOG.error(
                            "can't find %s %s in the id hash", marker_key, symbol)
                    # check "Not Applicable" -> "reference_locus"
                    mapped_marker_type = self.resolve(marker_type.strip())

                    # if it's unlocated, or is not a gene,
                    # then don't add it as a class because
                    # it's not added as a gene.
                    # everything except for genes are modeled as individuals

                    if mapped_marker_type in [
                            self.globaltt['gene'],
                            self.globaltt['pseudogene']]:
                        model.addClassToGraph(
                            marker_id, symbol, mapped_marker_type, name)
                        model.addSynonym(
                            marker_id, name, self.globaltt['has_exact_synonym'])
                        self.markers['classes'].append(marker_id)
                    else:
                        model.addIndividualToGraph(
                            marker_id, symbol, mapped_marker_type, name)
                        model.addSynonym(
                            marker_id, name, self.globaltt['has_exact_synonym'])
                        self.markers['indiv'].append(marker_id)

                    self.label_hash[marker_id] = symbol
                    # add the taxon  (default to Mus m.)
                    # latin_name is not always a proper binomial
                    if latin_name == '"Not Applicable':  # localtt conflict
                        latin_name = 'Mus musculus'
                    taxon_id = self.resolve(
                        latin_name, default=self.globaltt['Mus musculus'])
                    geno.addTaxon(taxon_id, marker_id)

                    # make MGI the leader for mouse genes.
                    if taxon_id == self.globaltt['Mus musculus']:
                        model.makeLeader(marker_id)

                    if not self.test_mode and limit is not None \
                            and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        LOG.info("getting markers and equivalent ids from mrk_summary_view")
        line_counter = 0
        raw = '/'.join((self.rawdir, 'mrk_summary_view'))
        with open(raw, 'r') as fh:
            fh.readline()  # read the header row; skip
            for line in fh:
                line = line.rstrip("\n")
                line_counter += 1

                (accid, logicaldb_key, object_key, preferred,
                 mgiid, subtype, short_description) = line.split('\t')

                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                if preferred == '1':

                    if self.idhash['marker'].get(object_key) is None:
                        # can't find the marker in the hash; add it here:
                        self.idhash['marker'][object_key] = mgiid
                        LOG.error(
                            "this marker hasn't been seen before %s %s",
                            mgiid, short_description)

                    if accid == mgiid:
                        # don't need to make equivalences to itself
                        continue

                    mapped_id = None
                    if logicaldb_key == '60':
                        mapped_id = 'ENSEMBL:' + accid
                    elif logicaldb_key == '1':
                        # don't need to add the equivalence to itself.
                        continue
                    elif logicaldb_key == '55':
                        mapped_id = 'NCBIGene:' + accid

                    if mapped_id is not None:
                        if mgiid in self.markers['classes'] \
                                or subtype in ['Gene', 'Pseudogene']:
                            model.addClassToGraph(mapped_id, None)
                            model.addEquivalentClass(mgiid, mapped_id)
                        elif mgiid in self.markers['indiv']:
                            model.addIndividualToGraph(mapped_id, None)
                            model.addSameIndividual(mgiid, mapped_id)

                    # could parse the "subtype" string
                    # to get the kind of thing the marker is

                if not self.test_mode and limit is not None and line_counter > limit:
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
        LOG.info("mapping markers to internal identifiers")
        raw = '/'.join((self.rawdir, 'mrk_acc_view'))
        col = [
            'accid', 'prefix_part', 'logicaldb_key', 'object_key', 'preferred',
            'organism_key']
        with open(raw, 'r') as fh:
            fh.readline()  # read the header row; skip
            for line in fh:
                line = line.rstrip('\n')
                line_counter += 1
                row = line.split('\t')

                accid = row[col.index('accid')]
                prefix_part = row[col.index('prefix_part')]
                logicaldb_key = row[col.index('logicaldb_key')]
                object_key = row[col.index('object_key')]
                preferred = row[col.index('preferred')]
                # organism_key)

                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                # get the hashmap of the identifiers
                if logicaldb_key == '1' and prefix_part == 'MGI:' and preferred == '1':
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        # pass through the file again,
        # and make the equivalence statements to a subset of the idspaces.
        # TODO verify the difference between what the
        # mrk_acc_view vs mrk_summary_view buys us here.
        # if nothing, then we should remove one or the other.
        LOG.info("mapping marker equivalent identifiers in mrk_acc_view")
        line_counter = 0
        with open('/'.join((self.rawdir, 'mrk_acc_view')), 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line = line.rstrip("\n")
                line_counter += 1
                (accid, prefix_part, logicaldb_key, object_key,
                 preferred, organism_key) = line.split('\t')

                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                # right now not caring about other organisms
                if organism_key != 1:
                    continue

                mgiid = self.idhash['marker'].get(object_key)
                if mgiid is None:
                    # presumably we've already added the relevant MGI ids,
                    # so skip those that we can't find
                    LOG.debug("can't find mgiid for %s", object_key)
                    continue
                marker_id = None
                if preferred == '1':  # TODO what does it mean if it's 0?
                    if logicaldb_key == '55':  # entrez/ncbi
                        marker_id = 'NCBIGene:' + accid
                    elif logicaldb_key == '1' and prefix_part != 'MGI:':
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
                        LOG.error("mgiid not in class or indiv hash %s", mgiid)

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        LOG.info("mapping strains to internal identifiers")
        raw = '/'.join((self.rawdir, 'prb_strain_acc_view'))

        tax_id = self.globaltt["Mus musculus"]

        with open(raw, 'r') as fh:
            fh.readline()  # read the header row; skip
            for line in fh:
                line = line.rstrip("\n")
                line_counter += 1
                (accid, prefixpart, logicaldb_key, object_key, preferred) \
                    = line.split('\t')
                # scrub out the backticks from accids
                # TODO notify the source upstream
                accid = re.sub(r'`', '', accid).strip()
                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('strain'):
                        continue

                # get the hashmap of the identifiers
                if logicaldb_key == '1' and prefixpart == 'MGI:' and preferred == '1':
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
        LOG.info("mapping strain equivalent identifiers")
        line_counter = 0
        with open(raw, 'r') as fh:
            fh.readline()  # read the header row; skip
            for line in fh:
                line = line.rstrip("\n")
                line_counter += 1
                (accid, prefixpart, logicaldb_key, object_key, preferred) \
                    = line.split('\t')
                # scrub out the backticks from accids
                # TODO notify the source upstream
                accid = re.sub(r'`', '', accid).strip()
                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('strain'):
                        continue
                mgiid = self.idhash['strain'].get(object_key)
                if mgiid is None:
                    # presumably we've already added the relevant MGI ids,
                    # so skip those that we can't find
                    # LOG.info("can't find mgiid for %s",object_key)
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

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        LOG.info("getting free text descriptions for annotations")
        raw = '/'.join((self.rawdir, 'mgi_note_vocevidence_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (object_key, note) = line

                if self.test_mode is True:
                    if int(object_key) not in self.test_keys.get('notes'):
                        continue
                # object_key == evidence._annotevidence_key
                annotkey = self.idhash['notes'].get(object_key)
                annot_id = self.idhash['annot'].get(annotkey)
                # only add the description for the annotations
                # we have captured through processing

                if annot_id is not None:
                    model.addDescription(annot_id, note.strip())

                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        return

    def _process_mrk_location_cache(self, limit):
        line_counter = 0
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        LOG.info("getting marker locations")
        raw = '/'.join((self.rawdir, 'mrk_location_cache'))
        geno = Genotype(graph)

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

                if self.test_mode is True:
                    if int(marker_key) not in self.test_keys.get('marker'):
                        continue

                # make the chromsomome, and the build-instance
                chrom_id = makeChromID(chromosome, 'NCBITaxon:10090', 'CHR')
                if version is not None and version != '' and version != '(null)':

                    # switch on maptype or mapkey
                    assembly = version
                    build_id = 'NCBIGenome:' + assembly
                    geno.addChromosomeInstance(
                        chromosome, build_id, assembly, chrom_id)
                    chrom_id = makeChromID(chromosome, build_id, 'MONARCH')

                if marker_key in self.idhash['marker']:
                    gene_id = self.idhash['marker'][marker_key]
                    feature = Feature(graph, gene_id, None, None)
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
                            [self.globaltt['FuzzyPosition']])
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
                    LOG.warning('marker key %s not in idhash', str(marker_key))

                if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        LOG.info("getting transgene genes")
        raw = '/'.join((self.rawdir, 'mgi_relationship_transgene_genes'))
        geno = Genotype(graph)
        col = [
            'rel_key', 'allele_key', 'allele_id', 'allele_label', 'category_key',
            'category_name', 'property_key', 'property_name', 'gene_num'
        ]
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            header = next(filereader)
            if header != col:
                LOG.error('expected columns:  %s\n\tBut got:\n%s', col, header)
            for row in filereader:
                # rel_key,
                allele_key = int(row[col.index('allele_key')])
                allele_id = row[col.index('allele_id')]
                # allele_label,
                # category_key,
                # category_name,
                # property_key,
                # property_name,
                gene_num = int(row[col.index('gene_num')])

                if self.test_mode and allele_key not in self.test_keys.get('allele')\
                        and gene_num not in self.test_ids:
                    continue

                gene_id = 'NCBIGene:' + str(gene_num)

                # geno.addParts(gene_id, allele_id, self.globaltt['has_variant_part'])
                seqalt_id = self.idhash['seqalt'].get(allele_key)
                if seqalt_id is None:
                    seqalt_id = allele_id
                geno.addSequenceDerivesFrom(seqalt_id, gene_id)

                if not self.test_mode and limit is not None and \
                        filereader.line_num > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        model = Model(graph)
        LOG.info("Assembling notes on alleles")
        raw = '/'.join((self.rawdir, 'mgi_note_allele_view'))

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
                            int(sequencenum)
                    ):
                        notehash[object_key][notetype].append('')  # ??? I don't get it

                notehash[object_key][notetype][int(sequencenum)-1] = note.strip()

            # finish iteration over notes

        line_counter = 0
        for allele_key in notehash:
            if self.test_mode is True:
                if int(allele_key) not in self.test_keys.get('allele'):
                    continue
            line_counter += 1
            allele_id = self.idhash['allele'].get(allele_key)
            if allele_id is None:
                continue
            for n in notehash[allele_key]:
                LOG.info(
                    "found %d %s notes for %s",
                    len(notehash[allele_key]), n, allele_id)
                notes = ''.join(notehash[allele_key][n])
                notes += ' ['+n+']'
                model.addDescription(allele_id, notes)

            if not self.test_mode and limit is not None and line_counter > limit:
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
        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph
        LOG.info("Getting genotypes for strains")
        raw = '/'.join((self.rawdir, 'prb_strain_genotype_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if line_counter == 1:
                    continue

                (strain_key, genotype_key) = line

                if self.test_mode is True:
                    if int(genotype_key) not in self.test_keys.get('genotype') \
                            and int(strain_key) not in self.test_keys.get('strain'):
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

                graph.addTriple(strain_id, self.globaltt['has_genotype'], genotype_id)
                # TODO
                # verify if this should be contingent on the exactness or not
                # if qualifier == 'Exact':
                #     gu.addTriple(
                #       graph, strain_id,
                #       self.globaltt['has_genotype'],
                #       genotype_id)
                # else:
                #     gu.addXref(graph, strain_id, genotype_id)

                if not self.test_mode and limit is not None and line_counter > limit:
                    break

        return

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
        # these are just blank nodes
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
