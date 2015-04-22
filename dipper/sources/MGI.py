import csv
import os
from datetime import datetime
import logging
import re

from dipper.sources.Source import Source
from dipper.models.Assoc import Assoc
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.utils import CurieUtil
from dipper import config
from dipper import curie_map
from dipper.utils.GraphUtils import GraphUtils

logger = logging.getLogger(__name__)


class MGI(Source):
    """
    This is the Mouse informatics resource, from which we process genotype and phenotype data about laboratory mice.
    Genotypes leverage the GENO genotype model.
    Be sure to have pg user/password connection details in your conf.json file, like:
      dbauth : {
        'mgi' : {'user' : '<username>', 'password' : '<password>'}
      }
    """

    #CONSIDER IF WE NEED:
    # mgi_organism_acc_view: Consider using this for the taxon mapping instead of the hashmap encoded below
    # mgi_reference_allele_view: Don't believe this view is used in either the genotype of phenotype view
    # all_allele_cellline_view: When we want to start dealing with cell lines
    # mgi_note_strain_view: prose descriptons of strains.
    # prb_strain_summary_view: Don't believe this view is used in either the genotype of phenotype view
    # prb_strain_marker_view: eventually i think we want this because it has other relevant markers that are affected

    tables = [
        'mgi_dbinfo',
        'gxd_genotype_view',
        'gxd_genotype_summary_view',
        'gxd_allelepair_view',
        'all_summary_view',
        'all_allele_view',
        'all_allele_mutation_view',
        'mrk_marker_view',
        'voc_annot_view',
        'voc_evidence_view',
        'bib_acc_view',
        'prb_strain_view',
        'mrk_summary_view',
        'mrk_acc_view',
        'prb_strain_acc_view',
        'mgi_note_vocevidence_view'
    ]


    # for testing purposes, this is a list of internal db keys to match and select only portions of the source
    test_keys = {
        'allele' : [1303,56760,816699,51074,14595,816707,246,38139,4334,817387,8567,476,42885,3658,1193,6978,6598,16698],
        'marker' : [38043,305574,444020,34578,9503,38712,17679,445717,38415,12944,377,77197,18436,30157,14252],
        'annot' : [23546035,189443,189450,29645664,29645665,189447,189452,189451,189446,29645667,189449,29645666,
                   189445,189444,189442,29645682,189448,29645663,69611011,93548463,928426,23534428,23535949,43838073,
                   43838073,318424,717023,717025,717027,717028,717029,717026,90942389,717024,90942390,90942392,
                   90942384,90942382,90942386,90942381,90942385,90942391,5647502,69611253,93436975,93436976,93436973,
                   93419639,93436974,93436977,93401080,84201418,6778,93636629,93648755,93614265,93607816,93624755,
                   43803707,43803707,80436328,93484431,93484432,93552054,13487622,13487624,13487623,43804057,43804057,
                   83942519,12035,6620086,93620355,93581890,79642481,93579091,93581841,93584586,93626409,79655585,
                   93618579,93581832,93579870,93576058,93581813,93587213,93604448,93583073,23241933,93583786,93643814,
                   43805682,43805682,92947717,92948518,92947729,92947735,92949200,92947757,92948169,92949301,92948441,
                   93491336,93491334,93491335,93491333,93491337,93551440,24722398,93642680,6173781,93459094,93652704,
                   6173780,93459097,93459095,93613038,93092371,93092368,93092369,93092372,93092373,58485679,93092374,
                   93626918,93643825,6173778,93092375,93647695,62628962,6173775,93459096,93092376,93092377,93092378,
                   93092379,93092380,93621390,43815003,43815003,93092370,93092381,93092382,93510300,93510296,93510297,
                   93510299,93510298,59357866,59357864,59357865,59357867,59357863,60448186,60448185,60448187],
        'genotype' : [81,87,142,206,281,283,286,287,341,350,384,406,407,411,425,457,458,461,476,485,537,546,551,
                      553,11702,12910,13407,13453,14815,26655,28610,37313,38345,59766,60082],
        'pub' : [73197,165659,134151,76922,181903,26681,128938,80054,156949,159965,53672,170462,206876,87798,100777,
                 176693,139205,73199,74017,102010,152095,18062,216614,61933,13385,32366,114625,182408,140802],
        'strain' : [30639,33832,33875,33940,36012,59504,34338,34382,47670,59802,33946,31421,64,40,14,-2,30639,15975,
                    35077,12610,-1,28319,27026,141],
        'notes' : [5114,107221415,1055341,37833486,107158128,107218642,53310,53311,53312,53313,53314,53315,53316,53317,
                   53318,53319,53320,71099,501751,501752,501753,103920323,501754,103920328,501755,103920312,501756,
                   103920325,501757,744108,6049949,103920319,6621213,107251870,6621216,6621218,6621219,7108498,
                   107168838,107248964,14590363,14590364,14590365,25123358,25123360,26688159,32028545,32028546,
                   32028547,32028548,32028549,32028564,47742903,47743253,47744878,47754199,47777269,65105483,
                   66144014,66144015,66144016,66144017,66144018,70046116,78382808,78383050,107154485,107237686,
                   107174867,107218519,107214911,107256603,106949909,106969369,103920318,103920320,103920322,
                   103920324,103920326,103920330,103920331,103920332,103920333,106390006,106390018,106390024,
                   106390046,106390458,106390730,106390807,106391489,106391590,106579450,106579451,106579452,
                   106579453,106579454,106579455,106579456,106579457,106579458,106579459,106579460,106579461,
                   106579462,106579463,106579464,106949910,106969368,106996040,106996041,106996042,106996043,
                   106996044,107022123,107022124,107022125,107022126,107052057,107052058,107058959,107058960,
                   107058961,107058962,107058963,107077922,107077923,107077924,107077925,107077926,107116089,
                   107119066,107119680,107155254,107159385,107160435,107163154,107163183,107163196,107163271,
                   107164877,107165872,107166942,107170557,107194346,107198590,107205179,107206725,107212120,
                   107214364,107215700,107219974,107222064,107222717,107235068,107242709,107244121,107244139,
                   107249091,107250401,107255383]
    }

    def __init__(self):
        Source.__init__(self, 'mgi')
        self.namespaces.update(curie_map.get())

        # update the dataset object with details about this resource
        self.dataset = Dataset('mgi', 'MGI', 'http://www.informatics.jax.org/')

        # check if config exists; if it doesn't, error out and let user know
        if (not (('dbauth' in config.get_config()) and ('mgi' in config.get_config()['dbauth']))):
            logger.error("not configured with PG user/password.")

        # source-specific warnings.  will be cleared when resolved.
        logger.warn("we are ignoring normal phenotypes for now")

        # Set this flag to False for normal running
        self.testMode = True

        # so that we don't have to deal with BNodes, we will create hash lookups for the internal identifiers
        # the hash will hold the type-specific-object-keys to MGI public identifiers.  then, subsequent
        # views of the table will lookup the identifiers in the hash.  this allows us to do the 'joining' on the
        # fly
        self.idhash = {'allele' : {}, 'marker' : {}, 'publication' : {}, 'strain' : {},
                       'genotype' : {}, 'annot' : {}, 'notes' : {}}
        self.markers = {'classes' : [], 'indiv' : []}  # to store if a marker is a class or indiv


        return

    def fetch(self, is_dl_forced):
        """
        For the MGI resource, we connect to the remote database, and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        """

        # create the connection details for MGI
        cxn = config.get_config()['dbauth']['mgi']
        cxn.update({'host' : 'adhoc.informatics.jax.org', 'database' : 'mgd', 'port' : 5432 })

        self.dataset.setFileAccessUrl(('').join(('jdbc:postgresql://', cxn['host'],':', str(cxn['port']), '/',
                                                 cxn['database'])))

        # process the tables
        # self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables,cxn)

        datestamp = ver = None
        # get the resource version information from table mgi_dbinfo, already fetched above
        outfile = ('/').join((self.rawdir,'mgi_dbinfo'))

        if os.path.exists(outfile):
            st = os.stat(outfile)
            with open(outfile, 'r') as f:
                f.readline()  # read the header row; skip
                info = f.readline()
                cols = info.split('\t')
                ver = cols[0]  # col 0 is public_version
                ver = ver.replace('MGI ', '')  # MGI 5.20 --> 5.20
                # MGI has a datestamp for the data within the database; use it instead of the download date
                # datestamp in the table: 2014-12-23 00:14:20
                d = cols[7].strip()  # modification date
                datestamp = datetime.strptime(d, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
                f.close()

        self.dataset.setVersion(datestamp,ver)

        return


    def parse(self, limit=None):
        """
        We process each of the postgres tables in turn.  The order of processing is important here, as we build
        up a hashmap of internal vs external identifers (unique keys by type to MGI id).  These include
        allele, marker (gene), publication, strain, genotype, annotation (association), and descriptive notes.
        :param limit: Only parse this many lines of each table
        :return:
        """
        if (limit is not None):
            logger.info("Only parsing first %d rows of each file", limit)
        logger.info("Parsing files...")

        loops = [True]
        if not self.testOnly:
            loops = [True,False]

        for l in loops:

            self.testMode = l

            # the following will provide us the hash-lookups
            # These must be processed in a specific order
            self._process_prb_strain_acc_view(limit)
            self._process_mrk_acc_view(limit)
            self._process_all_summary_view(limit)
            self._process_bib_acc_view(limit)
            self._process_gxd_genotype_summary_view(limit)

            # The following will use the hash populated above to lookup the ids when filling in the graph
            self._process_prb_strain_view(limit)
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

        logger.info("Finished parsing.")

        #TODO generate report of internal identifiers we created (eg for strains)

        self.load_bindings()
        for g in [self.graph,self.testgraph]:
            Assoc().loadAllProperties(g)

        logger.info("Loaded %d nodes", len(self.graph))
        return


    def _process_gxd_genotype_view(self, limit=None):
        """
        This table indicates the relationship between a genotype and it's background strain.  It leverages the
        Genotype class methods to do this.

        Makes these triples:
        <MGI:genotypeid> GENO:has_reference_part <MGI:strainid>
        <MGI:strainid> a GENO:genomic_background

        If the genotype id isn't in the hashmap, it adds it here (but this shouldn't happen):
        <MGI:genotypeid> a GENO:genotype

        If the strain isn't in the hashmap, it also adds it here with a monarchized identifier using the
        unique key of the strain, formatted like:  :_mgistrainkey12345

        :param limit:
        :return:
        """

        gu = GraphUtils(curie_map.get())
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)
        raw = ('/').join((self.rawdir, 'gxd_genotype_view'))
        logger.info("getting genotypes and their backgrounds")
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1
                (genotype_key, strain_key, isconditional, note, existsas_key, createdby_key, modifiedby_key,
                 creation_date, modification_date, strain, mgiid, dbname, createdbymodifiedby, existsas,
                 empty) = line.split('\t')

                if self.testMode is True:
                    if int(genotype_key) not in self.test_keys.get('genotype'):
                        continue

                if self.idhash['genotype'].get(genotype_key) is None:
                    # just in case we haven't seen it before, catch and add the id mapping here
                    self.idhash['genotype'][genotype_key] = mgiid
                    geno.addGenotype(mgiid, None)
                    # TODO get label

                # if it's in the hash, assume that the individual was created elsewhere
                strain_id = self.idhash['strain'].get(strain_key)
                if (strain_id is None):
                    # some of the strains don't have public identifiers!
                    # so we make one up, and add it to the hash
                    logger.warn("adding background as internal id: %s %s", strain_key, strain)
                    strain_id = ':'+self._makeInternalIdentifier('strain', strain_key)
                    # add it back to the idhash
                    self.idhash['strain'].update({strain_key : strain_id})
                    geno.addGenotype(strain_id, strain)

                geno.addGenomicBackgroundToGenotype(strain_id, mgiid)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_gxd_genotype_summary_view(self,limit=None):
        """
        Add the genotype internal id to mgiid mapping to the idhashmap.  Also, add them as individuals to the graph.
        We re-format the label to put the background strain in brackets after the gvc.

        We must pass through the file once to get the ids and aggregate the vslcs into a hashmap into the genotype

        Triples created:
        <genotype id> a GENO:intrinsic_genotype
        <genotype id> rdfs:label "<gvc> [bkgd]"

        :param limit:
        :return:
        """
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        gu = GraphUtils(curie_map.get())
        line_counter = 0
        geno_hash = {}
        raw = ('/').join((self.rawdir,'gxd_genotype_summary_view'))
        logger.info("building labels for genotypes")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,createdby_key,modifiedby_key,
                 creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('genotype'):
                        continue

                #add the internal genotype to mgi mapping
                self.idhash['genotype'][object_key] = mgiid

                if (preferred == '1'):
                    d = re.sub('\,','/',short_description.strip())
                    if mgiid not in geno_hash:
                        geno_hash[mgiid] = {'vslcs' : [d],'subtype' : subtype}
                    else:
                        vslcs = geno_hash[mgiid].get('vslcs')
                        vslcs.append(d)
                else:
                    pass
                    #TODO what to do with != preferred

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        #now, loop through the hash and add the genotypes as individuals
        gutil = Genotype(g)
        for g in geno_hash:
            geno = geno_hash.get(g)
            gvc = sorted(geno.get('vslcs'))
            label = ('; ').join(gvc) + '[' + geno.get('subtype') +']'
            gutil.addGenotype(g,label)

        #TODO materialize a GVC here?

        return

    #NOTE: might be best to process alleles initially from the all_allele_view, as this does not have any repeats of alleles!
    def _process_all_summary_view(self,limit):
        """
        Here, we get the allele definitions: id, label, description, type
        We add the id to this source's global idhash for lookup later

        <alleleid> a OWL:NamedIndividual
            rdf:label "allele symbol"
            dc:descirption "long allele name"

        :param limit:
        :return:
        """

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        line_counter = 0
        raw = ('/').join((self.rawdir,'all_summary_view'))
        logger.info("getting alleles and their labels and descriptions")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                 createdby_key,modifiedby_key,creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')
                #NOTE:May want to filter alleles based on the preferred field (preferred = 1) or will get duplicates
                ## (24288, to be exact... Reduced to 480 if filtered on preferred = 1)

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('allele'):
                        continue

                #we are setting the allele type to None, so that we can add the type later
                # (since we don't actually know if it's a reference or altered allele)
                altype = None  #temporary; we'll assign the type later

                #If we want to filter on preferred:
                if preferred == '1':
                    #add the allele key to the hash for later lookup
                    self.idhash['allele'][object_key] = mgiid
                    #TODO consider not adding the individuals in this one
                    gu.addIndividualToGraph(g,mgiid,short_description.strip(),altype,description.strip())

                #TODO deal with non-preferreds, are these deprecated?

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return


    def _process_all_allele_view(self,limit):
        """
        Add the allele as a variant locus (or reference locus if wild-type).
        If the marker is specified, we add the link to the marker.
        We assume that the MGI ids are available in the idhash, added in all_summary_view.
        We add the sequence alteration as a BNode here, if there is a marker.  Otherwise, the
        allele itself is a sequence alteration.

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
        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        geno = Genotype(g)
        line_counter = 0
        logger.info("adding alleles, mapping to markers, extracting their sequence alterations")
        raw = ('/').join((self.rawdir,'all_allele_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allele_key,marker_key,strain_key,mode_key,allele_type_key,allele_status_key,transmission_key,
                 collection_key,symbol,name,nomensymbol,iswildtype,isextinct,ismixed,createdby_key,modifiedby_key,
                 approvedby_key,approval_date,creation_date,modification_date,markersymbol,term,statusnum,strain,collection,createdby,modifiedby,approvedby) = line.split('\t')

                if self.testMode is True:
                    if int(allele_key) not in self.test_keys.get('allele'):
                        continue


                allele_id = self.idhash['allele'].get(allele_key)
                if (allele_id is None):
                    logger.error("what to do! can't find allele_id. skipping %s %s", allele_key, symbol)
                    continue

                if (marker_key is not None) and (marker_key != ''):
                    #we make the assumption here that the markers have already been added to the table
                    marker_id = self.idhash['marker'].get(marker_key)
                    if (marker_id is None):
                        logger.error("what to do! can't find marker_id. skipping %s %s",marker_key,symbol)
                        continue

                strain_id = self.idhash['strain'].get(strain_key)
                iseqalt_id = self._makeInternalIdentifier('seqalt',allele_key)
                if (self.testMode):
                    #in test mode, we want to make these identified nodes
                    iseqalt_id = ':'+iseqalt_id
                iseqalt = gu.getNode(iseqalt_id)
                # for non-wild type alleles:
                if iswildtype == '0':
                    locus_type = geno.genoparts['variant_locus']
                    locus_rel = geno.properties['is_sequence_variant_instance_of']
                #for wild type alleles:
                elif iswildtype == '1':
                    locus_type = geno.genoparts['reference_locus']
                    locus_rel = geno.properties['is_reference_instance_of']

                gu.addIndividualToGraph(g,allele_id,symbol,locus_type)
                al = gu.getNode(allele_id)

                #marker_id will be none if the allele is not linked to a marker (as in, it's not mapped to a locus)
                if marker_id is not None:
                    #add link between gene and allele
                    geno.addAlleleOfGene(allele_id,marker_id,locus_rel)

                #sequence alteration in strain
                #FIXME change this to a different relation in_strain, genomically_related_to, sequence_derives_from
                if iswildtype == '0':
                    sa_label = symbol
                    sa_id = iseqalt_id

                    if marker_key is not None and marker_key != '':
                        #sequence alteration has label reformatted(symbol)
                        if re.match(".*<.*>.*", symbol):
                            sa_label = re.sub(".*<", "<", symbol)
                        elif re.match("\+", symbol):
                            #TODO: Check to see if this is the proper handling, as while symbol is just +, marker symbol has entries without any <+>.
                            sa_label = '<+>'
                        geno.addSequenceAlterationToVariantLocus(iseqalt_id,allele_id)
                    else:
                        #make the sequence alteration == allele
                        sa_id = allele_id

                    #else this will end up adding the non-located transgenes as sequence alterations also
                    #removing the < and > from sa
                    sa_label = re.sub('[\<\>]','',sa_label)

                    #gu.addIndividualToGraph(g,sa_id,sa_label,None,name)
                    geno.addSequenceAlteration(sa_id,sa_label,None,name)
                    if strain_id is not None:
                        geno.addDerivesFrom(allele_id,strain_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_gxd_allele_pair_view(self,limit):
        """
        This assumes that the genotype and alleles have already been added to the id hashmap.
        We use the Genotype methods to add all the parts we need.
        Triples added:
        <genotype_id> has_part <vslc>
        <vslc> has_part <allele1>
        <vslc> has_part <allele2>
        <vslc> has_zygosity <zygosity>

        :param limit:
        :return:
        """

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        line_counter = 0
        geno = Genotype(g)
        raw = ('/').join((self.rawdir,'gxd_allelepair_view'))
        logger.info("processing allele pairs (VSLCs) for genotypes")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allelepair_key,genotype_key,allele_key_1,allele_key_2,marker_key,mutantcellline_key_1,mutantcellline_key_2,
                 pairstate_key,compound_key,sequencenum,createdby_key,modifiedby_key,creation_date,modification_date,symbol,
                 chromosome,allele1,allele2,allelestate,compound) = line.split('\t')
                #NOTE: symbol = gene/marker, allele1 + allele2 = VSLC, allele1/allele2 = variant locus, allelestate = zygosity
                #FIXME Need to handle alleles not in the *<*> format, such as many gene traps, induced mutations, and transgenics

                if self.testMode is True:
                    if int(genotype_key) not in self.test_keys.get('genotype'):
                        continue


                genotype_id = self.idhash['genotype'].get(genotype_key)

                if genotype_id is None:
                    logger.error("genotype_id not found for key %s; skipping", genotype_key)
                    continue

                allele1_id = self.idhash['allele'].get(allele_key_1)
                allele2_id = self.idhash['allele'].get(allele_key_2)

                #Need to map the allelestate to a zygosity term
                zygosity_id = self._map_zygosity(allelestate)
                ivslc_id = self._makeInternalIdentifier('vslc',allelepair_key)
                if self.testMode:
                    #make this a real id in test mode
                    ivslc_id = ':'+ivslc_id

                #TODO: VSLC label likely needs processing similar to the processing in the all_allele_view
                #FIXME: handle null alleles
                vslc_label = (allele1+'/'+allele2)

                gu.addIndividualToGraph(g, ivslc_id, vslc_label, geno.genoparts['variant_single_locus_complement'])
                geno.addVSLCtoParent(ivslc_id,genotype_id)
                geno.addPartsToVSLC(ivslc_id,allele1_id,allele2_id,zygosity_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_all_allele_mutation_view(self,limit):
        """
        This fetches the mutation type for the alleles, and maps them to the sequence alteration.
        Note that we create a BNode for the sequence alteration because it isn't publically identified.
        <sequence alteration id> a <SO:mutation_type>

        :param limit:
        :return:
        """

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        line_counter = 0
        raw = ('/').join((self.rawdir,'all_allele_mutation_view'))
        logger.info("getting mutation types for sequence alterations")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allele_key,mutation_key,creation_date,modification_date,mutation) = line.split('\t')
                iseqalt_id = self._makeInternalIdentifier('seqalt',allele_key)

                if self.testMode is True:
                    if int(allele_key) not in self.test_keys.get('allele'):
                        continue
                    iseqalt_id = ':'+iseqalt_id

                #TODO we might need to map the seq alteration to the MGI id for unlocated things; need to use hashmap
                #map the sequence_alteration_type
                seq_alt_type_id = self._map_seq_alt_type(mutation)

                gu.addIndividualToGraph(g,iseqalt_id,None,seq_alt_type_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return



    def _process_voc_annot_view(self,limit):
        """
        This MGI table represents associations between things.
        We currently filter this table on abnormal Genotype-Phenotype associations, but may be expanded in the future.

        We add the internal annotation id to the idhashmap.  It is expected that the genotypes have already
        been added to the idhash

        :param limit:
        :return:
        """

        #TODO also get Strain/Attributes (annottypekey = 1000)
        #TODO what is Phenotype (Derived) vs non-derived?  (annottypekey = 1015)
        #TODO is evidence in this table?  what is the evidence vocab key?

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        line_counter = 0
        logger.info("getting G2P associations")
        raw = ('/').join((self.rawdir,'voc_annot_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:

                (annot_key,annot_type_key,object_key,term_key,qualifier_key,creation_date,modification_date,qualifier,
                 term,sequence_num,accid,logicaldb_key,vocab_key,mgi_type_key,evidence_vocab_key,anot_type) = line.split('\t')

                if self.testMode is True:
                    if int(annot_key) not in self.test_keys.get('annot'):
                        continue

                # Restricting to type 1002, as done in the MousePhenotypes view.
                # Corresponds to 'Mammalian Phenotype/Genotype' and MP terms
                if annot_type_key == '1002':
                    line_counter += 1

                    #TODO add NOT annotations
                    #skip 'normal'
                    if (qualifier=='norm'):
                        logger.info("found normal phenotype: %s",term)
                        continue

                    # We expect the label for the phenotype to be taken care of elsewhere
                    gu.addClassToGraph(g,accid,None)

                    iassoc_id = self._makeInternalIdentifier('annot',annot_key)
                    assoc_id = ':'+self.make_id(iassoc_id)
                    #add the assoc to the hashmap (using the monarch id)
                    self.idhash['annot'][annot_key] = assoc_id
                    genotype_id = self.idhash['genotype'].get(object_key)
                    if (genotype_id is None):
                        logger.error("can't find genotype id for %s",object_key)
                    else:
                        #add the association
                        assoc = G2PAssoc(assoc_id,genotype_id,accid,None,None)
                        assoc.addAssociationNodeToGraph(g)


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return


    def _process_voc_evidence_view(self,limit):
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

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        line_counter = 0
        logger.info("getting evidence and pubs for annotations")
        raw = ('/').join((self.rawdir,'voc_evidence_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (annot_evidence_key,annot_key,evidence_term_key,refs_key,inferred_from,created_by_key,modified_by_key,
                creation_date,modification_date,evidence_code,evidence_seq_num,jnumid,jnum,short_citation,created_by,modified_by)= line.split('\t')

                if self.testMode is True:
                    if int(annot_key) not in self.test_keys.get('annot'):
                        continue

                #add the association id to map to the evidence key (this is to attach the right note to the right assn)
                self.idhash['notes'][annot_evidence_key] = annot_key

                assoc_id = self.idhash['annot'].get(annot_key)

                if (assoc_id is None):
                    #assume that we only want to add the evidence/source for annots that we have in our db
                    continue

                evidence_id = self._map_evidence_id(evidence_code)

                #TODO add it as an instance of what type?
                #add the pub as an individual;
                gu.addIndividualToGraph(g,jnumid,None)

                #add the ECO and citation information to the annot
                Assoc().addEvidence(g,evidence_id,assoc_id)
                Assoc().addSource(g,assoc_id,jnumid)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return


    def _process_bib_acc_view(self,limit):
        """
        This traverses the table twice:
           once to look up the internal key to J number mapping for the id hashmap
           then again to make the equivalences.  All internal keys have both a J and MGI identifier.
        This will make equivalences between the different pub ids
        Triples:
            <pub_id> a owl:NamedIndividual
            <other_pub_id> a owl:NamedIndividual
            <pub_id> owl:sameAs <other_pub_id>
        :param limit:
        :return:
        """

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        #TODO we may consider limiting the publication nodes
        #firstpass, get the J number mapping, and add to the global hash
        line_counter = 0
        logger.info('populating pub id hash')
        raw = ('/').join((self.rawdir,'bib_acc_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if (line_counter == 1):
                    continue #skip header
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                created_by_key,modified_by_key,creation_date,modification_date,logical_db)= line

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('pub'):
                        continue

                #we use the J number here because it is the externally-accessible identifier
                if prefixpart != 'J:':
                    continue
                self.idhash['publication'][object_key] = accid
                gu.addIndividualToGraph(g,accid,None)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        #2nd pass, look up the MGI identifier in the hash
        logger.info("getting pub equivalent ids")
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if (line_counter == 1):
                    continue #skip header
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                created_by_key,modified_by_key,creation_date,modification_date,logical_db)= line

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('pub'):
                        continue

                logical_db = logical_db.strip()

                jid = self.idhash['publication'].get(object_key)

                pub_id = None
                if (logicaldb_key == '29'):  #pubmed
                    pub_id = 'PMID:'+accid
                elif (logicaldb_key == '1' and re.match('MGI:',prefixpart)):
                    #don't get the J numbers, because we dont' need to make the equiv to itself.
                    pub_id = accid
                elif (logical_db == 'Journal Link'):
                    #some DOIs seem to have spaces
                    #FIXME MGI needs to FIX THESE UPSTREAM!!!!
                    #we'll scrub them here for the time being
                    pub_id = 'DOI:'+re.sub('\s+','',accid)
                elif (logicaldb_key == '1' and re.match('J:',prefixpart)):
                    #we can skip the J numbers
                    continue

                if (pub_id is not None):
                    #only add these to the graph if it's mapped to something we understand
                    gu.addIndividualToGraph(g,pub_id,None)
                    gu.addSameIndividual(g,jid,pub_id)
                else:
                    logger.warn("Publication from (%s) not mapped for %s", logical_db, object_key)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_prb_strain_view(self,limit):
        """
        Process a table to get strains (with internal ids), and their labels.
        These strains are created as instances of intrinsic_genotype.
        Triples:
            <strain id> a GENO:intrinsic_genotype
                rdf:label "strain label"
                RO:in_taxon <NCBI taxon id>

        :param limit:
        :return:
        """
        #Only 9 strain types if we want to map them (recombinant congenci, inbred strain, NA, congenic,
        # consomic, coisogenic, recombinant inbred, NS, conplastic)


        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        raw = ('/').join((self.rawdir,'prb_strain_view'))
        logger.info("getting strains and adding their taxa")
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if (line_counter == 1):
                    continue
                (strain_key,species_key,strain_type_key,strain,standard,private,genetic_background,created_by_key,
                modified_by_key,creation_date,modification_date,species,strain_type,created_by,modified_by) = line

                if self.testMode is True:
                    if int(strain_key) not in self.test_keys.get('strain'):
                        continue


                strain_id = self.idhash['strain'].get(strain_key)

                if strain_id is not None:
                    geno.addGenotype(strain_id,strain)

                    #add the species to the graph as a class
                    sp = self._map_strain_species(species)
                    gu.addClassToGraph(g,sp,None)

                    geno.addTaxon(sp,strain_id)

                    #TODO what is mgi's strain type anyway?

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return


    def _process_mrk_marker_view(self,limit):
        """
        This is the definition of markers (as in genes, but other genomic loci types as well).
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

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        geno = Genotype(g)
        line_counter = 0
        raw = ('/').join((self.rawdir,'mrk_marker_view'))
        logger.info("getting markers and assigning types")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (marker_key,organism_key,marker_status_key,marker_type_key,curationstate_key,symbol,name,chromosome,
                cytogenetic_offset,createdby_key,modifiedby_key,creation_date,modification_date,organism,common_name,
                latin_name,status,marker_type,curation_state,created_by,modified_by) = line.split('\t')

                if self.testMode is True:
                    if int(marker_key) not in self.test_keys.get('marker'):
                        continue

                #use only non-withdrawn markers
                if marker_status_key != '2':
                    marker_id = self.idhash['marker'].get(marker_key)

                    #only pull info for mouse genes for now  other species should come from other dbs
                    if (organism_key != '1'):
                        continue

                    if (marker_id is None):
                        logger.error("can't find %s %s in the id hash", marker_key, symbol)

                    mapped_marker_type = self._map_marker_type(marker_type)

                    #if it's unlocated, then don't add it as a class because it's not added as a gene.
                    #everything except for genes are modeled as individuals
                    if (chromosome is not None and chromosome.strip() != 'UN'):
                        gu.addClassToGraph(g,marker_id,symbol,mapped_marker_type,name)
                        gu.addSynonym(g,marker_id,name,Assoc.properties['hasExactSynonym'])
                        self.markers['classes'].append(marker_id)
                    else:
                        gu.addIndividualToGraph(g,marker_id,symbol,mapped_marker_type,name)
                        gu.addSynonym(g,marker_id,name,Assoc.properties['hasExactSynonym'])
                        self.markers['indiv'].append(marker_id)

                    #add the taxon
                    taxon_id = self._map_taxon(latin_name)
                    geno.addTaxon(taxon_id,marker_id)

                    if (not self.testMode) and (limit is not None and line_counter > limit):
                        break

        return


    def _process_mrk_summary_view(self,limit):
        """
        Here we pull the mgiid of the features, and make equivalent (or sameAs) associations to referenced ids.
        Only adding the ENSEMBL genes and NCBI gene ids.  Will wait on other ids later.

        :param limit:
        :return:
        """

        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        logger.info("getting markers and equivalent ids from mrk_summary_view")
        line_counter = 0
        raw = ('/').join((self.rawdir,'mrk_summary_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgi_type_key,private,preferred,
                 created_by_key,modified_by_key,creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                if (preferred=='1'):

                    if self.idhash['marker'].get(object_key) is None:
                        #can't find the marker in the hash; add it here:
                        self.idhash['marker'][object_key] = mgiid
                        logger.error("this marker hasn't been seen before %s %s",mgiid,short_description)


                    if (accid == mgiid):
                        #don't need to make equivalences to itself
                        continue


                    mapped_id = None
                    if logicaldb_key == '60':
                        mapped_id = 'ENSEMBL:'+accid
                    elif logicaldb_key == '1':
                        continue #don't need to add the equivalence to itself.
                    elif logicaldb_key == '55':
                        mapped_id = 'NCBIGene:'+accid
                    #elif logicaldb_key == '41':
                    #    mapped_id = 'TrEMBL:'+accid
                    #elif logicaldb_key == '45':
                    #    mapped_id = 'PDB:'+accid
                    #elif logicaldb_key == '135':
                    #    mapped_id = 'PR:'+accid
                    #ISSUE: MirBase accession ID can map to multiple MGI IDs if the miRNA is also part of a cluster (Mirlet7b is part of cluster Mirc31)
                    #elif logicaldb_key == '83':
                    #    mapped_id = 'miRBase:'+accid
                    #elif logicaldb_key == '27':
                    #    mapped_id = 'RefSeq:'+accid
                    #elif logicaldb_key == '13':
                    #    mapped_id = 'SwissProt:'+accid
                    #elif logicaldb_key == '8':
                    #    mapped_id = 'EC:'+accid
                    #FIXME: The EC IDs are used for multiple genes, resulting in one EC number

                    if (mapped_id is not None):
                        if (mgiid in self.markers['classes'] or subtype == 'Gene'):
                            gu.addClassToGraph(g,mapped_id,None)
                            gu.addEquivalentClass(g,mgiid,mapped_id)
                        elif (mgiid in self.markers['indiv']):
                            gu.addIndividualToGraph(g,mapped_id,None)
                            gu.addSameIndividual(g,mgiid,mapped_id)


                    #TODO add links to other ensembl things: logical_db #134 = protein, #133 = transcript, #60 is gene

                    #could parse the "subtype" string to get the kind of thing the marker is

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_mrk_acc_view(self,limit):
        """
        Use this table to create the idmap between the internal marker id and the public mgiid.
        No triples are produced in this process
        :param limit:
        :return:
        """

        #make a pass through the table first, to create the mapping between the external and internal identifiers
        line_counter = 0
        logger.info("mapping markers to internal identifiers")
        raw = ('/').join((self.rawdir,'mrk_acc_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefix_part,numeric_part,logicaldb_key,object_key,mgi_type_key,private,preferred,
                 created_by_key,modified_by_key,creation_date,modification_date,logicaldb,organism_key) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                #get the hashmap of the identifiers
                if (logicaldb_key == '1') and (prefix_part == 'MGI:') and (preferred == '1'):
                    self.idhash['marker'][object_key] = accid

        return

    def _process_mrk_acc_view_for_equiv(self,limit):
        """
        Add the equivalences, either sameAs or equivalentClass, depending on the nature of the marker.  We
        only process the ENSEMBL genes and NCBI gene ids.
        :param limit:
        :return:
        """
        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        #pass through the file again, and make the equivalence statements to a subset of the idspaces.
        #TODO verify the difference between what the mrk_acc_view vs mrk_summary_view buys us here.
        #if nothing, then we should remove one or the other.
        logger.info("mapping marker equivalent identifiers in mrk_acc_view")
        line_counter = 0
        with open(('/').join((self.rawdir,'mrk_acc_view')), 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefix_part,numeric_part,logicaldb_key,object_key,mgi_type_key,private,preferred,
                 created_by_key,modified_by_key,creation_date,modification_date,logicaldb,organism_key) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('marker'):
                        continue

                #right now not caring about other organisms
                if (organism_key != 1):
                    continue

                mgiid = self.idhash['marker'].get(object_key)
                if (mgiid is None):
                    #presumably we've already added the relevant MGI ids already, so skip those that we can't find
                    logger.debug("can't find mgiid for %s",object_key)
                    continue
                marker_id = None
                if (preferred == '1'):  #what does it mean if it's 0?
                    if logicaldb_key == '55':  #entrez/ncbi
                        marker_id = 'NCBIGene:'+accid
                    elif logicaldb_key == '1' and prefix_part != 'MGI:':  #mgi
                        marker_id = accid
                    elif logicaldb_key == '60':
                        marker_id = 'ENSEMBL:'+accid
                    #TODO get other identifiers
                    #TODO get non-preferred ids==deprecated?

                if (marker_id is not None):
                    if(mgiid in self.markers['classes']):
                        gu.addClassToGraph(g,marker_id,None)
                        gu.addEquivalentClass(g,mgiid,marker_id)
                    elif (mgiid in self.markers['indiv']):
                        gu.addIndividualToGraph(g,marker_id,None)
                        gu.addSameIndividual(g,mgiid,marker_id)
                    else:
                        logger.error("mgiid not in class or indiv hash %s",mgiid)



                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return

    def _process_prb_strain_acc_view(self,limit):
        """
        Use this table to create the idmap between the internal marker id and the public mgiid.
        Also, add the equivalence statements between strains for MGI and JAX
        Triples:
        <strain_id> a GENO:intrinsic_genotype
        <other_strain_id> a GENO:intrinsic_genotype
        <strain_id> owl:sameAs <other_strain_id>

        :param limit:
        :return:
        """

        #make a pass through the table first, to create the mapping between the external and internal identifiers
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph

        geno = Genotype(g)
        logger.info("mapping strains to internal identifiers")
        raw = ('/').join((self.rawdir,'prb_strain_acc_view'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,
                 preferred,createdby_key,modifiedby_key,creation_date,modification_date,logicaldb) = line.split('\t')

                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('strain'):
                        continue

                #get the hashmap of the identifiers
                if (logicaldb_key == '1') and (prefixpart == 'MGI:') and (preferred == '1'):
                    self.idhash['strain'][object_key] = accid
                    gu.addIndividualToGraph(g,accid,None,geno.genoparts['intrinsic_genotype'])

        #pass through the file again, and make the equivalence statements to a subset of the idspaces
        logger.info("mapping strain equivalent identifiers")
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,
                 preferred,createdby_key,modifiedby_key,creation_date,modification_date,logicaldb) = line.split('\t')


                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('strain'):
                        continue


                mgiid = self.idhash['strain'].get(object_key)
                if (mgiid is None):
                    #presumably we've already added the relevant MGI ids already, so skip those that we can't find
                    #logger.info("can't find mgiid for %s",object_key)
                    continue
                strain_id = None
                if (preferred == '1'):  #what does it mean if it's 0?
                    if logicaldb_key == '22':  #JAX
                        #scrub out the backticks from accids
                        #TODO notify the source upstream
                        accid = re.sub('`','',accid)
                        strain_id = 'JAX:'+accid
                    #TODO get non-preferred ids==deprecated?

                if (strain_id is not None):
                    gu.addIndividualToGraph(g,strain_id,None,geno.genoparts['intrinsic_genotype'])
                    gu.addSameIndividual(g,mgiid,strain_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        return


    def _process_mgi_note_vocevidence_view(self,limit):
        """
        Here we fetch the free text descriptions of the phenotype associations.
        Triples:
        <annot_id> dc:description "description text"
        :param limit:
        :return:
        """

        line_counter = 0
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        logger.info("getting free text descriptions for annotations")
        raw = ('/').join((self.rawdir,'mgi_note_vocevidence_view'))
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if (line_counter == 1):
                    continue

                (note_key,object_key,mgitype_key,notetype_key,createdby_key,modifiedby_key,
                 creation_date,modification_date,notetype,note,sequencenum) = line


                if self.testMode is True:
                    if int(object_key) not in self.test_keys.get('notes'):
                        continue

                annotkey = self.idhash['notes'].get(object_key)  # object_key == evidence._annotevidence_key
                annot_id = self.idhash['annot'].get(annotkey)
                #only add the description for the annotations we have captured through processing
                if (annot_id is not None):
                    assoc = Assoc()
                    assoc.addDescription(g,annot_id,note.strip())

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break


        return

    # def _process_genoaggregate(self,raw,limit):
    #     """
    #     Here, we add the aggregated genotype labels to the graph.
    #     This is the table that is created by us in a selection and aggregation of MGI data
    #     in order to get a single label for a given genotype.
    #     :param raw:
    #     :param limit:
    #     :return:
    #     """
    #
    #     gu = GraphUtils(curie_map.get())
    #
    #     line_counter = 0
    #     print("INFO: getting free text descriptions for annotations")
    #     with open(raw, 'r', encoding="utf8") as csvfile:
    #         filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
    #         for line in filereader:
    #             line_counter += 1
    #             if (line_counter == 1):
    #                 continue
    #
    #             (mgiid,genotype_key,gvc,subtype) = line
    #
    #             if self.testMode is True:
    #                 if int(genotype_key) not in self.test_keys.get('genotype'):
    #                     continue
    #
    #            #add the internal genotype to mgi mapping
    #             self.idhash['genotype'][genotype_key] = mgiid
    #
    #             if (preferred == '1'):
    #                 geno_label = gvc.strip() + '[' + subtype + ']'
    #                 gu.addIndividualToGraph(g,mgiid,geno_label,self.terms['intrinsic_genotype'])
    #
    #             #TODO what to do with != preferred
    #             #TODO note the short_description is the GVC  (use this or reason?)
    #
    #             if (limit is not None and line_counter > limit):
    #                 break
    #
    #
    #     return


    #TODO generalize this to a set of utils
    def _getcols(self,cur,table):
        query=(' ').join(("SELECT * FROM",table,"LIMIT 0"))  #for testing
        cur.execute(query)
        colnames = [desc[0] for desc in cur.description]
        logger.info("columns ("+table+"):",colnames)

        return


    def file_len(self,fname):
        with open(fname) as f:
            return sum(1 for line in f)


    #TODO: Finish identifying SO/GENO terms for mappings for those found in MGI
    def _map_seq_alt_type(self, sequence_alteration_type):
        type = 'SO:0001059'  #default to sequence_alteration
        type_map = {
            'Deletion': 'SO:0000159',  # deletion
            'Disruption caused by insertion of vector': 'SO:0000667',  # insertion - correct?
            'Duplication': 'SO:1000035',  # duplication
            'Insertion': 'SO:0000667',  # insertion
            'Insertion of gene trap vector': 'SO:0001218',  # transgenic insertion - correct?  (TODO gene_trap_construct: SO:0001477)
            'Intergenic deletion': 'SO:0000159',  # deletion  #TODO return a list? SO:0001628 intergenic_variant
            'Intragenic deletion': 'SO:0000159',  # deletion  #TODO return a list?  SO:0001564 gene_variant
            'Inversion': 'SO:1000036',  # inversion
            'Not Applicable': 'SO:0001059',
            'Not Specified': 'SO:0001059',
            'Nucleotide repeat expansion': 'SO:1000039',  # tandem duplication  #TODO ask for another term
            'Nucleotide substitutions': 'SO:0002007',  # multiple nucleotide variant
            'Other': 'SO:0001059',
            'Single point mutation': 'SO:1000008',  # point_mutation
            'Translocation': 'SO:0000199',  # translocation
            'Transposon insertion': 'SO:0001837',  # mobile element insertion
            'Undefined': 'SO:0001059',
            'Viral insertion': 'SO:0001838',  # novel sequence insertion (no viral version)
            'wild type': 'SO:0000817'  # wild type
        }
        if (sequence_alteration_type.strip() in type_map):
            type = type_map.get(sequence_alteration_type.strip())
        else:
            # TODO add logging
            logger.error("Sequence Alteration Type (%s) not mapped; defaulting to sequence_alteration", sequence_alteration_type)

        return type

    def _map_zygosity(self, zygosity):
        type = None
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
        if (zygosity.strip() in type_map):
            type = type_map.get(zygosity)
        else:
            logger.error("Zygosity (%s) not mapped", zygosity)

        return type

    def _map_marker_type(self, marker_type):
        type = None
        type_map = {
            'Complex/Cluster/Region': 'SO:0000001',  # region. Something more specific available? #fixme
            'Transgene': 'SO:0000902',  # transgene
            'Gene': 'SO:0000704',  # gene
            'QTL': 'SO:0000771',  # QTL
            'DNA Segment': 'SO:0000110',  # sequence_feature. sequence_motif=SO:0001683? region=SO:0000001
            'Pseudogene': 'SO:0000336',  # pseudogene
            'Cytogenetic Marker': 'SO:0001645',  # genetic_marker?   #fixme
            'Other Genome Feature': 'SO:0000110',  # sequence_feature. Or sequence_motif=SO:0001683?
            'BAC/YAC end': 'SO:0000150',  # BAC_end: SO:0000999, YAC_end: SO:00011498; using parent term
        }
        if (marker_type.strip() in type_map):
            type = type_map.get(marker_type)
        else:
            logger.error("Marker Type (%s) not mapped", marker_type)

        return type

    def _map_allele_type(self,allele_type):
        """
        This makes the assumption that all things are variant_loci (including Not Specified)
        :param allele_type:
        :return:
        """
        if self.testMode:
            g=self.testgraph
        else:
            g=self.graph
        geno = Genotype(g)
        type = geno.genoparts['variant_locus']  #assume it's a variant locus
        type_map = {
            'Not Applicable' : geno.genoparts['reference_locus'],
            'QTL': geno.genoparts['reference_locus'],  # should QTLs be something else?  or SO:QTL?
        }
        if (allele_type.strip() in type_map):
            type = type_map.get(allele_type)

        return type


    def _map_taxon(self, taxon_name):
        type = None
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
            'Mus musculus/domesticus': 'NCBITaxon:10090',  # 10090=Mus musculus, 10092=Mus musculus domesticus
            'Ornithorhynchus anatinus': 'NCBITaxon:9258',
            'Oryctolagus cuniculus': 'NCBITaxon:9986',
            'Ovis aries': 'NCBITaxon:9940',
            'Pan troglodytes': 'NCBITaxon:9598',
            'Pongo pygmaeus': 'NCBITaxon:9600',
            'Rattus norvegicus': 'NCBITaxon:10116',
            'Sus scrofa domestica L.': 'NCBITaxon:9823',  # 9823=Sus scrofa, 9825=Sus scrofa domestica
            'Xenopus (Silurana) tropicalis': 'NCBITaxon:8364',
        }
        if (taxon_name.strip() in type_map):
            type = type_map.get(taxon_name)
        else:
            logger.error("Taxon Name (%s) not mapped", taxon_name)

        return type

    def _map_strain_species(self,species):
        #make the assumption that it is a Mus genus, unless if specified
        tax = '10088'
        id_map = {
            'laboratory mouse and M. m. domesticus (brevirostris)': '116058',
            'laboratory mouse and M. m. bactrianus': '35531',
            'laboratory mouse and M. m. castaneus and M. m. musculus': '477816',  # does not include lab mouse
            'M. m. domesticus and M. m. molossinus and M. m. castaneus': '1266728',
            'M. setulosus': '10102',
            'laboratory mouse and wild-derived': '10088',  #Mus genus
            'laboratory mouse and M. m. musculus (Prague)': '39442',
            'M. m. castaneus and M. m. musculus': '477816',
            'M. m. domesticus (Canada)': '10092',
            'M. m. domesticus and M. m. domesticus poschiavinus': '10092',  #FIXME
            'M. m. musculus and M. spretus or M. m. domesticus': '186842',
            'M. chypre': '862507',  #unclassified mus
            'M. cookii (Southeast Asia)': '10098',
            'M. cervicolor': '10097',
            'M. m. gentilulus': '80274',
            'M. m. domesticus poschiavinus (Tirano, Italy)': '10092',  #FIXME
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
            'M. hortulanus (Austria)': '10103',  #synonymous with Mus spicilegus
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
            'M. hortulanus': '10103',  #synonymous with Mus spicilegus
            'M. dunni': '254704',  #synonymous with Mus terricolor
            'M. m. castaneus (Pathumthani, Thailand)': '10091',
            'M. terricolor': '254704',
            'M. m. domesticus (China)': '10092',
            'M. fragilicauda': '186193',
            'Not Resolved': '10088',  #assume Mus
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
            'Not Specified': '10090',  #OK?
            'Not Specified and M. m. molossinus': '10090',
            'M. m. musculus and M. spretus': '862507',
            'M. minutoides': '10105',
            'M. spretus (Tunisia)': '10096',
            'M. spicilegus': '10103',
            'Peru Coppock': '10088',  #unclassified mus
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
        if (species.strip() in id_map):
            tax = id_map.get(species.strip())

        else:
            logger.warn("Species (%s) not mapped; defaulting to Mus genus.", species)

        return 'NCBITaxon:'+tax



    def _map_evidence_id(self, evidence_code):
        #TODO a default evidence code???  what should it be?
        type = None
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
        if (evidence_code.strip() in type_map):
            type = type_map.get(evidence_code)
        else:
            logger.error("Evidence code (%s) not mapped", evidence_code)

        return type

    # def _map_evidence_label(self, evidence_code):
    #     """
    #     TODO remove this probably
    #     :param evidence_code:
    #     :return:
    #     """
    #     type = None
    #     type_map = {
    #         'EXP': 'experimental evidence',
    #         'IBA': 'biological aspect of ancestor evidence used in manual assertion',
    #         'IC': 'inference from background scientific knowledge',
    #         'IDA': 'direct assay evidence used in manual assertion',
    #         'IEA': 'evidence used in automatic assertion',
    #         'IEP': 'expression pattern evidence',
    #         'IGI': 'genetic interaction evidence used in manual assertion',
    #         'IKR': 'phylogenetic determination of loss of key residues evidence used in manual assertion',
    #         'IMP': 'mutant phenotype evidence used in manual assertion',
    #         'IPI': 'physical interaction evidence used in manual assertion',
    #         'ISA': 'sequence alignment evidence',
    #         'ISM': 'match to sequence model evidence',
    #         'ISO': 'sequence orthology evidence',
    #         'ISS': 'sequence similarity evidence used in manual assertion',
    #         'NAS': 'non-traceable author statement used in manual assertion',
    #         'ND': 'no biological data found',
    #         'RCA': 'computational combinatorial evidence used in manual assertion',
    #         'TAS': 'traceable author statement used in manual assertion'
    #     }
    #     if (evidence_code.strip() in type_map):
    #         type = type_map.get(evidence_code)
    #         # type = 'http://purl.obolibrary.org/obo/' + type_map.get(zygosity)
    #     # print("Mapped: ", allele_type, "to", type)
    #     else:
    #         # TODO add logging
    #         print("ERROR: Taxon Name (", evidence_code, ") not mapped")
    #
    #     return type


    def _makeInternalIdentifier(self,prefix,key):
        '''
        This is a special MGI-to-MONARCH-ism.  MGI tables have unique keys that we use here, but don't want
        to necessarily re-distribute those internal identifiers.  Therefore, we make them into keys in a consistent
        way here.
        :param prefix: the object type to prefix the key with, since the numbers themselves are not unique across tables
        :param key: the number
        :return:
        '''

        return '_mgi'+prefix+'key'+key

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
        #TODO test genotypes
        #from tests.test_genotypes import GenotypeTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MGITestCase)

        return test_suite