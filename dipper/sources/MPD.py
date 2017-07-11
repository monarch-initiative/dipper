import csv
import re
import logging
import io
import gzip
from dipper.models.Provenance import Provenance

from dipper.sources.Source import Source
from dipper.models.Genotype import Genotype
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Model import Model

logger = logging.getLogger(__name__)


class MPD(Source):
    """
    From the [MPD](http://phenome.jax.org/) website:
    This resource is a collaborative standardized collection of measured data
    on laboratory mouse strains and populations. Includes baseline phenotype
    data sets as well as studies of drug, diet, disease and aging effect.
    Also includes protocols, projects and publications, and SNP,
    variation and gene expression studies.

    Here, we pull the data and model the genotypes using GENO and
    the genotype-to-phenotype associations using the OBAN schema.

    MPD provide measurements for particular assays for several strains.
    Each of these measurements is itself mapped to a MP or VT term
    as a phenotype.  Therefore, we can create a strain-to-phenotype association
    based on those strains that lie outside of the "normal" range for the given
    measurements.  We can compute the average of the measurements
    for all strains tested, and then threshold any extreme measurements being
    beyond some threshold beyond the average.

    Our default threshold here, is +/-2 standard deviations beyond the mean.

    Because the measurements are made and recorded at the level of
    a specific sex of each strain, we associate the MP/VT phenotype with
    the sex-qualified genotype/strain.

    """
    MPDDL = 'http://phenomedoc.jax.org/MPD_downloads'
    files = {
        'ontology_mappings': {
            'file': 'ontology_mappings.csv',
            'url': MPDDL + '/ontology_mappings.csv'},
        'straininfo': {
            'file': 'straininfo.csv',
            'url': MPDDL + '/straininfo.csv'},
        'assay_metadata': {
            'file': 'measurements.csv',
            'url': MPDDL + '/measurements.csv'},
        'strainmeans': {
            'file': 'strainmeans.csv.gz',
            'url': MPDDL + '/strainmeans.csv.gz'},
        # 'mpd_datasets_metadata': { #TEC does not seem to be used
        #    'file': 'mpd_datasets_metadata.xml.gz',
        #    'url': MPDDL + '/mpd_datasets_metadata.xml.gz'},
    }

    # the following are strain ids for testing
    # test_ids = [
    #   "MPD:2", "MPD:3", "MPD:5", "MPD:6", "MPD:9", "MPD:11", "MPD:18",
    #   "MPD:20", "MPD:24", "MPD:28", "MPD:30", "MPD:33", "MPD:34", "MPD:36",
    #   "MPD:37", "MPD:39", "MPD:40", "MPD:42", "MPD:47", "MPD:66", "MPD:68",
    #   "MPD:71", "MPD:75", "MPD:78", "MPD:122", "MPD:169", "MPD:438",
    #   "MPD:457","MPD:473", "MPD:481", "MPD:759", "MPD:766", "MPD:770",
    #   "MPD:849",  "MPD:857", "MPD:955", "MPD:964", "MPD:988", "MPD:1005",
    #   "MPD:1017", "MPD:1204", "MPD:1233", "MPD:1235", "MPD:1236", "MPD:1237"]

    test_ids = [
        'MPD:6', 'MPD:849', 'MPD:425', 'MPD:569', "MPD:10", "MPD:1002",
        "MPD:39", "MPD:2319"]

    mgd_agent_id = "MPD:db/q?rtn=people/allinv"
    mgd_agent_label = "Mouse Phenotype Database"
    mgd_agent_type = "foaf:organization"

    def __init__(self, graph_type, are_bnodes_skolemized):
        Source.__init__(self, graph_type, are_bnodes_skolemized, 'mpd')
        # @N, not sure if this step is required
        self.stdevthreshold = 2

        # update the dataset object with details about this resource
        # @N: Note that there is no license as far as I can tell
        self.dataset = Dataset(
            'mpd', 'MPD', 'http://phenome.jax.org', None, None)

        # TODO add a citation for mpd dataset as a whole
        self.dataset.set_citation('PMID:15619963')

        self.assayhash = {}
        self.idlabel_hash = {}
        # to store the mean/zscore of each measure by strain+sex
        self.score_means_by_measure = {}
        # to store the mean value for each measure by strain+sex
        self.strain_scores_by_measure = {}

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        """
        MPD data is delivered in four separate csv files and one xml file,
        which we process iteratively and write out as
        one large graph.

        :param limit:
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %s rows fo each file", str(limit))

        logger.info("Parsing files...")

        self._process_straininfo(limit)
        # the following will provide us the hash-lookups
        # These must be processed in a specific order

        # mapping between assays and ontology terms
        self._process_ontology_mappings_file(limit)
        # this is the metadata about the measurements
        self._process_measurements_file(limit)
        # get all the measurements per strain
        self._process_strainmeans_file(limit)

        # The following will use the hash populated above
        # to lookup the ids when filling in the graph
        self._fill_provenance_graph(limit)

        logger.info("Finished parsing.")
        return

    def _process_ontology_mappings_file(self, limit):

        # line_counter = 0  # TODO unused

        logger.info("Processing ontology mappings...")
        raw = '/'.join((self.rawdir, 'ontology_mappings.csv'))

        with open(raw, 'r') as f:
            reader = csv.reader(f)
            # read the header row; skip
            f.readline()
            for row in reader:
                try:
                    (assay_id, ont_term, descrip) = row
                except ValueError:
                    continue
                assay_id = int(assay_id)
                if re.match(r'(MP|VT)', ont_term):
                    # add the mapping denovo
                    if assay_id not in self.assayhash:
                        self.assayhash[assay_id] = {}
                        self.assayhash[assay_id]['ont_terms'] = set()
                    self.assayhash[assay_id]['ont_terms'].add(ont_term)

        return

    def _process_straininfo(self, limit):
        # line_counter = 0  # TODO unused
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)

        logger.info("Processing measurements ...")
        raw = '/'.join((self.rawdir, self.files['straininfo']['file']))

        tax_id = 'NCBITaxon:10090'

        with open(raw, 'r') as f:
            reader = csv.reader(f, delimiter=',', quotechar='\"')
            f.readline()  # read the header row; skip
            for row in reader:
                (strain_name, vendor, stocknum, panel, mpd_strainid,
                 straintype, n_proj, n_snp_datasets, mpdshortname, url) = row
                # C57BL/6J,J,000664,,7,IN,225,17,,http://jaxmice.jax.org/strain/000664.html
                # create the strain as an instance of the taxon
                if self.testMode and \
                        'MPD:' + str(mpd_strainid) not in self.test_ids:
                    continue
                strain_id = 'MPD-strain:' + str(mpd_strainid)
                model.addIndividualToGraph(strain_id, strain_name, tax_id)
                if mpdshortname.strip() != '':
                    model.addSynonym(strain_id, mpdshortname.strip())
                self.idlabel_hash[strain_id] = strain_name
                # make it equivalent to the vendor+stock
                if stocknum != '':
                    if vendor == 'J':
                        jax_id = 'JAX:'+stocknum
                        model.addSameIndividual(strain_id, jax_id)
                    elif vendor == 'Rbrc':
                        # reiken
                        reiken_id = 'RBRC:'+re.sub(r'RBRC', '', stocknum)
                        model.addSameIndividual(strain_id, reiken_id)
                    else:
                        if url != '':
                            model.addXref(strain_id, url, True)
                        if vendor != '':
                            model.addXref(
                                strain_id, ':'.join((vendor, stocknum)),
                                True)

                # add the panel information
                if panel != '':
                    desc = panel+' [panel]'
                    model.addDescription(strain_id, desc)

                # TODO make the panels as a resource collection

        return

    def _process_measurements_file(self, limit):
        line_counter = 0

        logger.info("Processing measurements ...")
        raw = '/'.join((self.rawdir, 'measurements.csv'))

        with open(raw, 'r') as f:
            reader = csv.reader(f)
            # read the header row; skip
            header = f.readline()
            logger.info("HEADER: %s", header)
            for row in reader:
                # measnum,projsym,varname,descrip,units,cat1,cat2,cat3,
                # intervention,intparm,appmeth,panelsym,datatype,sextested,
                # nstrainstested,ageweeks
                # Again the last row has changed. contains: '(4486 rows)'
                if len(row) != 16:
                    continue
                line_counter += 1
                assay_id = int(row[0])
                assay_label = row[3]
                assay_units = row[4]
                assay_type = row[10] if row[10] is not '' else None

                if assay_id not in self.assayhash:
                    self.assayhash[assay_id] = {}
                description = self.build_measurement_description(row)
                self.assayhash[assay_id]['description'] = description
                self.assayhash[assay_id]['assay_label'] = assay_label
                self.assayhash[assay_id]['assay_type'] = assay_type
                self.assayhash[assay_id]['assay_units'] = assay_units

                # TODO add projectsym property?
                # TODO add intervention?
                # ageweeks might be useful for adding to phenotype assoc

            # end loop on measurement metadata

        return

    def _process_strainmeans_file(self, limit):
        """
        This will store the entire set of strain means in a hash.
        Not the most efficient representation,
        but easy access.
        We will loop through this later to then apply cutoffs
        and add associations
        :param limit:
        :return:

        """
        logger.info("Processing strain means ...")
        line_counter = 0
        raw = '/'.join((self.rawdir, self.files['strainmeans']['file']))
        with gzip.open(raw, 'rb') as f:
            f = io.TextIOWrapper(f)
            reader = csv.reader(f)
            f.readline()  # read the header row; skip
            score_means_by_measure = {}
            strain_scores_by_measure = {}
            for row in reader:
                try:
                    (measnum, varname, strain, strainid, sex, mean, nmice, sd,
                     sem, cv, minval, maxval, logmean, logsd, zscore,
                     logzscore) = row
                except ValueError:
                    continue
                line_counter += 1
                strain_num = int(strainid)
                assay_num = int(measnum)
                # assuming the zscore is across all the items
                # in the same measure+var+strain+sex
                # note: it seems that there is only ever 1 varname per measnum.
                # note: some assays only tested one sex!
                # we split this here by sex
                if assay_num not in score_means_by_measure:
                    score_means_by_measure[assay_num] = {}
                if sex not in score_means_by_measure[assay_num]:
                    score_means_by_measure[assay_num][sex] = list()
                score_means_by_measure[assay_num][sex].append(float(mean))

                if strain_num not in strain_scores_by_measure:
                    strain_scores_by_measure[strain_num] = {}
                if sex not in strain_scores_by_measure[strain_num]:
                    strain_scores_by_measure[strain_num][sex] = {}
                strain_scores_by_measure[strain_num][sex][assay_num] = \
                    {'mean': float(mean), 'zscore': float(zscore)}

            # end loop over strainmeans
        self.score_means_by_measure = score_means_by_measure
        self.strain_scores_by_measure = strain_scores_by_measure

        return

    def _fill_provenance_graph(self, limit):
        logger.info("Building graph ...")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        model = Model(g)
        taxon_id = 'NCBITaxon:10090'  # hardcode to Mus musculus
        model.addClassToGraph(taxon_id, None)

        scores_passing_threshold_count = 0
        scores_passing_threshold_with_ontologies_count = 0
        scores_not_passing_threshold_count = 0

        # loop through all the strains,
        # and make G2P assoc for those with scores beyond threshold
        for strain_num in self.strain_scores_by_measure:
            if self.testMode and 'MPD:'+str(strain_num) not in self.test_ids:
                continue
            strain_id = 'MPD-strain:'+str(strain_num)
            for sex in self.strain_scores_by_measure[strain_num]:
                measures = self.strain_scores_by_measure[strain_num][sex]
                for m in measures:
                    assay_id = 'MPD-assay:'+str(m)
                    # TODO consider using the means
                    # instead of precomputed zscores
                    if 'zscore' in measures[m]:
                        zscore = measures[m]['zscore']
                        if abs(zscore) >= self.stdevthreshold:
                            scores_passing_threshold_count += 1
                            # logger.info(
                            #   "Score passing threshold: %s | %s | %s",
                            #   strain_id, assay_id, zscore)
                            # add the G2P assoc
                            prov = Provenance(self.graph)
                            try:
                                assay_label = self.assayhash[m]['assay_label']
                                assay_description = \
                                    self.assayhash[m]['description']
                                ont_term_ids = self.assayhash[m].get('ont_terms')
                                comment = ' '.join((assay_label,
                                                   '(zscore='+str(zscore)+')'))
                            except KeyError:
                                assay_label = None
                                assay_description = None
                                ont_term_ids = None
                            if assay_label is not None:
                                assay_label += ' ('+str(m)+')'
                            # TODO unused
                            # assay_type = self.assayhash[m]['assay_type']

                            assay_type_id = Provenance.provenance_types['assay']

                            if ont_term_ids is not None:
                                scores_passing_threshold_with_ontologies_count += 1
                                prov.add_assay_to_graph(
                                    assay_id, assay_label, assay_type_id,
                                    assay_description)
                                self._add_g2p_assoc(
                                    g, strain_id, sex, assay_id, ont_term_ids,
                                    comment)
                        else:
                            scores_not_passing_threshold_count += 1

        logger.info("Scores passing threshold: %d",
                    scores_passing_threshold_count)
        logger.info("Scores passing threshold with ontologies: %d",
                    scores_passing_threshold_with_ontologies_count)
        logger.info("Scores not passing threshold: %d",
                    scores_not_passing_threshold_count)

        return

    def _add_g2p_assoc(self, g, strain_id, sex, assay_id, phenotypes, comment):
        """
        Create an association between a sex-specific strain id
        and each of the phenotypes.
        Here, we create a genotype from the strain,
        and a sex-specific genotype.
        Each of those genotypes are created as anonymous nodes.

        The evidence code is hardcoded to be:
            ECO:experimental_phenotypic_evidence.

        :param g:
        :param strain_id:
        :param sex:
        :param assay_id:
        :param phenotypes: a list of phenotypes to association with the strain
        :param comment:
        :return:

        """
        geno = Genotype(g)
        model = Model(g)
        eco_id = "ECO:0000059"  # experimental_phenotypic_evidence
        strain_label = self.idlabel_hash.get(strain_id)
        # strain genotype
        genotype_id = '_'+'-'.join((re.sub(r':', '', strain_id), 'genotype'))
        genotype_label = '[' + strain_label + ']'

        sex_specific_genotype_id = '_'+'-'.join((re.sub(r':', '', strain_id),
                                                 sex, 'genotype'))
        if strain_label is not None:
            sex_specific_genotype_label = strain_label + ' (' + sex + ')'
        else:
            sex_specific_genotype_label = strain_id + '(' + sex + ')'

        genotype_type = Genotype.genoparts['sex_qualified_genotype']
        if sex == 'm':
            genotype_type = Genotype.genoparts['male_genotype']
        elif sex == 'f':
            genotype_type = Genotype.genoparts['female_genotype']

        # add the genotype to strain connection
        geno.addGenotype(
            genotype_id, genotype_label,
            Genotype.genoparts['genomic_background'])
        g.addTriple(
            strain_id, Genotype.object_properties['has_genotype'], genotype_id)

        geno.addGenotype(
            sex_specific_genotype_id, sex_specific_genotype_label,
            genotype_type)

        # add the strain as the background for the genotype
        g.addTriple(
            sex_specific_genotype_id,
            Genotype.object_properties['has_sex_agnostic_genotype_part'],
            genotype_id)

        # #############    BUILD THE G2P ASSOC    #############
        # TODO add more provenance info when that model is completed

        if phenotypes is not None:
            for phenotype_id in phenotypes:
                assoc = G2PAssoc(
                    g, self.name, sex_specific_genotype_id, phenotype_id)
                assoc.add_evidence(assay_id)
                assoc.add_evidence(eco_id)
                assoc.add_association_to_graph()
                assoc_id = assoc.get_association_id()
                model.addComment(assoc_id, comment)

        return

    def getTestSuite(self):
        import unittest
        from tests.test_mpd import MPDTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(MPDTestCase)

        return test_suite

    @staticmethod
    def normalise_units(units):
        # todo:
        return units

    @staticmethod
    def build_measurement_description(row):
        (assay_id, projsym, varname, descrip, units, cat1, cat2, cat3,
         intervention, intparm, appmeth, panelsym, datatype, sextested,
         nstrainstested, ageweeks) = row

        if sextested == 'f':
            sextested = 'female'
        elif sextested == 'm':
            sextested = 'male'
        elif sextested == 'fm':
            sextested = 'male and female'
        else:
            logger.warning("Unknown sex tested key: %s", sextested)
        description = "This is an assay of [" + descrip + "] shown as a [" + \
                      datatype + "] measured in [" + units + "]"

        if intervention is not None and intervention != "":
            description += " in response to [" + intervention + "]"
        if intparm is not None and intervention != "":
            description += \
                ". This represents the [" + intparm + \
                "] arm, using materials and methods that included [" +\
                appmeth + "]"

        description += \
            ".  The overall experiment is entitled [" + projsym + "].  "

        description += \
            "It was conducted in [" + sextested + "] mice at [" + \
            ageweeks + "] of age in" + " [" + nstrainstested + \
            "] different mouse strains. "
        description += "Keywords: " + cat1 + \
                       ((", " + cat2) if cat2.strip() is not "" else "") + \
                       ((", " + cat3) if cat3.strip() is not "" else "") + "."
        return description

    # def _log_missing_ids(self, missing_id, name_of_file_from_which_missing):
    #     if missing_id not in self.missing_assay_hash:
    #         self.missing_assay_hash[missing_id] = set()
    #     self.missing_assay_hash[missing_id].add(name_of_file_from_which_missing)
    #     # todo: remove the offending ids from the hash
    #     return
