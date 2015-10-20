import csv
import gzip
import re
import logging
import io
import math
import zipfile
from zipfile import ZipFile
from dipper.models.Provenance import Provenance

from dipper.sources.Source import Source
from dipper.models.Genotype import Genotype
from dipper.models.Dataset import Dataset
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map
from subprocess import call

logger = logging.getLogger(__name__)


class MPD(Source):
    """
    From the [MPD](http://phenome.jax.org/) website:
    This resource is a collaborative standardized collection of measured data on laboratory mouse strains and
    populations. Includes baseline phenotype data sets as well as studies of drug, diet, disease and aging effect.
    Also includes protocols, projects and publications, and SNP, variation and gene expression studies.

    Here, we pull the data and model the genotypes using GENO and the genotype-to-phenotype associations
    using the OBAN schema.

    We use all identifiers given by MPD.

    """

    files = {

        'ontology_mappings': {'file': 'ontology_mappings.csv',
                              'url': 'http://phenome.jax.org/download/ontology_mappings.csv'},
        'straininfo': {'file': 'straininfo.csv',
                       'url': 'http://phenome.jax.org/download/straininfo.csv'},
        'assay_metadata': {'file': 'measurements.csv',
                           'url': 'http://phenome.jax.org/download/measurements.csv'},
        'strainmeans': {'file': 'strainmeans.zip',
                        'url': 'http://phenome.jax.org/download/strainmeans.zip'},
        'mpd_datasets_metadata': {'file': 'mpd_datasets_metadata.xml',
                                  'url': 'http://phenome.jax.org/download/mpd_datasets_metadata.xml'},

    }

    # the following are strain ids for testing
    test_ids = ["MPD:2", "MPD:3", "MPD:5", "MPD:6", "MPD:9", "MPD:11", "MPD:18", "MPD:20", "MPD:24", "MPD:28", "MPD:30",
                "MPD:33", "MPD:34", "MPD:36", "MPD:37", "MPD:39", "MPD:40", "MPD:42", "MPD:47", "MPD:66", "MPD:68",
                "MPD:71", "MPD:75", "MPD:78", "MPD:122", "MPD:169", "MPD:438", "MPD:457", "MPD:473", "MPD:481",
                "MPD:759", "MPD:766", "MPD:770", "MPD:849", "MPD:857", "MPD:955", "MPD:964", "MPD:988", "MPD:1005",
                "MPD:1017", "MPD:1204", "MPD:1233", "MPD:1235", "MPD:1236", "MPD:1237"]

    mgd_agent_id = "MPD:db/q?rtn=people/allinv"
    mgd_agent_label = "Mouse Phenotype Database"
    mgd_agent_type = "foaf:organization"

    def __init__(self):
        Source.__init__(self, 'mpd')
        # @N, not sure if this step is required
        self.namespaces.update(curie_map.get())
        self.stdevthreshold = 2

        # update the dataset object with details about this resource
        # @N: Note that there is no license as far as I can tell
        self.dataset = Dataset('mpd', 'MPD', 'http://phenome.jax.org', None, None)

        # TODO add a citation for mpd dataset as a whole
        self.dataset.set_citation('PMID:15619963')

        # so that we don't have to deal with BNodes, we will create hash lookups for the internal identifiers
        # the hash will hold the type-specific-object-keys to MPD public identifiers.  then, subsequent
        # views of the table will lookup the identifiers in the hash.  this allows us to do the 'joining' on the
        # fly
        self.assayhash = {}
        # self.selectedAssayHash = {}
        # produces objects that look like:
        # {
        #     'metadata': {'ont_terms': [], 'description': None, 'assay_label': None, 'assay_type': None,
        #                  'measurement_unit': None},
        #     # have to repeat this otherwise both f and m will correspond to same object in memory
        #     # storing these in a less hierarchical way as it facilitates manual computation of zscores
        #     'f': {'strain_ids': [], 'strain_labels': [], 'assay_means': [], 'assay_zscores': []},
        #     'm': {'strain_ids': [], 'strain_labels': [], 'assay_means': [], 'assay_zscores': []}
        # }

        self.assays_missing_phenotypes = list()

        # create a hash to hold problem identifiers (eg. those missing from some tables but not others)
        # the keys will be the assay ids and the value would be a list of files from which known to be missing
        self.missing_assay_hash = {}

        return

    def fetch(self, is_dl_forced=False):

        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        """
        MPD data is delivered in four separate csv files and one xml file, which we process iteratively and write out as
        one large graph.

        :param limit:
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %s rows fo each file", str(limit))

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_straininfo(limit)

        # the following will provide us the hash-lookups
        # These must be processed in a specific order
        self._process_ontology_mappings_file(limit)
        self._process_measurements_file(limit)
        self._process_strainmeans_file(limit)
        # self._calculate_stats()

        # The following will use the hash populated above to lookup the ids when filling in the graph
        self._fill_provenance_graph(limit)
        # self._fill_association_graph(limit)

        logger.info("Finished parsing.")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return

    def _process_ontology_mappings_file(self, limit):

        line_counter = 0

        logger.info("Processing ontology mappings...")
        raw = '/'.join((self.rawdir, 'ontology_mappings.csv'))

        with open(raw, 'r') as f:
            reader = csv.reader(f)
            f.readline()  # read the header row; skip
            for row in reader:
                line_counter += 1
                (assay_id, ont_term, descrip) = row
                assay_id = int(assay_id)

                if re.match('(MP|VT)', ont_term):
                    # add the mapping denovo
                    if assay_id not in self.assayhash:
                        self.assayhash[assay_id] = {
                            'metadata': {'ont_terms': [], 'description': None, 'assay_label': None, 'assay_type': None,
                                         'measurement_unit': None},
                            # have to repeat this otherwise both f and m will correspond to same object in memory
                            # storing these in a less hierarchical way as it facilitates manual computation of zscores
                            'f': {'strain_ids': [], 'strain_labels': [], 'assay_means': [], 'assay_zscores': []},
                            'm': {'strain_ids': [], 'strain_labels': [], 'assay_means': [], 'assay_zscores': []}
                        }
                    self.assayhash[assay_id]['metadata']['ont_terms'] += [ont_term]

                else:
                    self.assays_missing_phenotypes.append(assay_id)
                    # the assays themselves will be added to the graph later as part of provenance
        return

    def _process_straininfo(self, limit):
        line_counter = 0
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        logger.info("Processing measurements ...")
        raw = '/'.join((self.rawdir, self.files['straininfo']['file']))

        tax_id = 'NCBITaxon:10090'

        gu = GraphUtils(curie_map.get())

        with open(raw, 'r') as f:
            reader = csv.reader(f, delimiter=',', quotechar='\"')
            f.readline()  # read the header row; skip
            for row in reader:
                (strain_name, vendor, stocknum, panel, mpd_strainid, straintype,
                 n_proj, n_snp_datasets, mpdshortname, url) = row
                # C57BL/6J,J,000664,,7,IN,225,17,,http://jaxmice.jax.org/strain/000664.html
                # create the strain as an instance of the taxon
                strain_id = 'MPD-strain:'+str(mpd_strainid)
                gu.addIndividualToGraph(g, strain_id, strain_name, tax_id)
                gu.addSynonym(g, strain_id, mpdshortname)

                # make it equivalent to the vendor+stock
                if vendor == 'J':
                    jax_id = 'JAX:'+stocknum
                    gu.addSameIndividual(g, strain_id, jax_id)
                elif vendor == 'Rbrc':
                    # reiken
                    reiken_id = 'RBRC:'+re.sub('RBRC', '', stocknum)
                    gu.addSameIndividual(g, strain_id, reiken_id)
                else:
                    if url != '':
                        gu.addXref(g, strain_id, url, True)
                    if vendor != '' and stocknum != '':
                        gu.addXref(g, strain_id, ':'.join((vendor, stocknum)), True)

                # add the panel information
                if panel != '':
                    desc = panel+' [panel]'
                    gu.addDescription(g, strain_id, desc)

                # TODO make the panels as a resource collection

        return

    def _process_measurements_file(self, limit):
        line_counter = 0

        logger.info("Processing measurements ...")
        raw = '/'.join((self.rawdir, 'measurements.csv'))

        with open(raw, 'r') as f:
            reader = csv.reader(f)
            f.readline()  # read the header row; skip
            for row in reader:
                line_counter += 1
                # measnum,projsym,varname,descrip,units,cat1,cat2,cat3,intervention,intparm,appmeth,panelsym,datatype,sextested,nstrainstested,ageweeks
                assay_id = int(row[0])
                assay_label = row[3]
                assay_units = row[4]
                assay_type = row[10] if row[10] is not '' else None
                if assay_id in self.assayhash and len(self.assayhash[assay_id]['metadata']['ont_terms']) != 0:
                    description = self._build_measurement_description(row)
                    self.assayhash[assay_id]['metadata']['description'] = description
                    self.assayhash[assay_id]['metadata']['assay_label'] = assay_label
                    self.assayhash[assay_id]['metadata']['assay_type'] = assay_type
                    self.assayhash[assay_id]['metadata']['assay_units'] = assay_units
                else:
                    # else, if we haven't already discarded this assay_id due to the lack of an MP/VT ontology mapping...
                    if assay_id not in self.assays_missing_phenotypes:
                        # we should log the fact that it didn't have a corresponding mapping at all
                        if assay_id not in self.missing_assay_hash:
                            self._log_missing_ids(assay_id, "ontology_mappings.csv")

        return

    @staticmethod
    def _build_measurement_description(row):
        (assay_id, projsym, varname, descrip, units, cat1, cat2, cat3, intervention, intparm, appmeth, panelsym,
         datatype, sextested, nstrainstested, ageweeks) = row

        sextested = 'female' if (sextested is 'fm') else 'male'
        if sextested == 'f':
            sextested = 'female'
        elif sextested == 'm':
            sextested = 'male'
        elif sextested == 'fm':
            sextested = 'male and female'
        else:
            logger.warn("Unknown sex tested key: %s", sextested)
        description = "This is an assay of [" + descrip + "] shown as a [" + datatype + "] measured in [" + units + "]"

        if intervention is not None and intervention != "":
            description += " in response to [" + intervention + "]"
        if intparm is not None and intervention != "":
            description += ". This represents the [" + intparm + "] arm, using materials and methods that " \
                                                                 "included [" + appmeth + "]"

        description += ".  The overall experiment is entitled [" + projsym + "].  "

        description += "It was conducted in [" + sextested + "] mice at [" + ageweeks + "] of age in" \
                       + " [" + nstrainstested + "] different mouse strains. "
        description += "Keywords: " + cat1 + \
                       ((", " + cat2) if cat2.strip() is not "" else "") + \
                       ((", " + cat3) if cat3.strip() is not "" else "") + "."
        return description

    @staticmethod
    def normalise_units(units):
        # todo:
        return units

    def _process_strainmeans_file(self, limit):
        logger.info("Processing strain means ...")
        line_counter = 0
        f = '/'.join((self.rawdir, self.files['strainmeans']['file']))
        myzip = ZipFile(f, 'r')
        # assume that the first entry is the item
        # fname = myzip.namelist()[0]

        with myzip.open(self.files['strainmeans']['file'].replace('.zip', '.csv'), 'r') as f:
            f = io.TextIOWrapper(f)
            reader = csv.reader(f)
            f.readline()  # read the header row; skip
            for row in reader:
                line_counter += 1
                # measnum,varname,strain,strainid,sex,mean,nmice,sd,sem,cv,minval,maxval,logmean,logsd,zscore,logzscore
                vals = (assay_id, strain_label, strainid, sex, mean, zscore) = int(row[0]), row[2], int(row[3]), \
                                                                               row[4], float(row[5]), float(row[14])
                for val in vals:
                    self.check_vals(val)

                if self.testMode is True:
                    if int(strainid) not in self.test_ids:
                        continue

                # if the assay_id is missing from the hash, it is because it corresponds to data that was leaked
                # and retracted ontology_mappings or measurements
                if assay_id not in self.assayhash:
                    self._log_missing_ids(assay_id, "ontology_mappings.csv")
                    self._log_missing_ids(assay_id, "measurements.csv")
                    # make a map of maps to hold strains for each (assay_id+sex)...

                else:
                    if assay_id not in self.assays_missing_phenotypes:
                        self.assayhash[assay_id][sex]['strain_ids'].append(strainid)
                        self.assayhash[assay_id][sex]['strain_labels'].append(strain_label)
                        self.assayhash[assay_id][sex]['assay_means'].append(mean)
                        self.assayhash[assay_id][sex]['assay_zscores'].append(zscore)
                        logger.debug(
                            str(assay_id) + ':' + sex + ' Strain IDs: '
                            + str(self.assayhash[assay_id][sex]['strain_ids'])
                            + ' Assay means: ' + str(self.assayhash[assay_id][sex]['assay_means']))
        return

    @staticmethod
    def check_vals(value):

        assert (value is not None and value is not '')

        return value

    # # loop over the resulting hash to calculate the z scores for each assay_id + sex
    # def _calculate_stats(self):
    #     logger.info("Calculating stats ...")
    #
    #     sexes = ['f', 'm']
    #     for sex in sexes:
    #         for assay_id in self.assayhash:
    #             strains = self.assayhash[assay_id][sex]['strain_id']
    #             # if strains is not an empty list
    #             if (strains):
    #                 logger.info("Present: " + str(assay_id) + ":" + sex)
    #
    #                 # if (assay_id == 23201):
    #                 #     logger.info()
    #
    #                 meanslist = self.assayhash[assay_id][sex]['assay_means']
    #                 zscores = self.zify(meanslist)
    #                 self.assayhash[sex]['assay_zscores'] = dict(zip(strains, zscores))
    #             # if strains is empty list it is because the corresponding assay_id isn't in strainmeans.csv
    #             else:
    #                 self._log_missing_ids(assay_id, "strainmeans.csv")
    #                 logger.info("Missing: " + str(assay_id) + ":" + sex)
    #
    #     return

    # # @N:
    # # I keep getting these errors of positional arguments (must be my misunderstanding of 'self')
    # # All stats code below is based on the statistics module in Python 3.4.
    # def zify(self, meanslist):
    #     pop_mean = self.mean(meanslist)
    #     pop_stdev = self.pstdev(meanslist)
    #
    #     zscores = list()
    #
    #     for mean in meanslist:
    #         zscore = (mean - pop_mean) / pop_stdev
    #         zscores.append(zscore)
    #     return zscores
    #
    # def mean(self, data):
    #     """Return the sample arithmetic mean of data."""
    #     n = len(data)
    #     if n < 1:
    #         raise ValueError('mean requires at least one data point')
    #     return sum(data) / float(n)  # in Python 2 use sum(data)/float(n)
    #
    # def _ss(self, data):
    #     """Return sum of square deviations of sequence data."""
    #     c = self.mean(data)
    #     ss = sum((x - c) ** 2 for x in data)
    #     return ss
    #
    # def pstdev(self, data):
    #     """Calculates the population standard deviation."""
    #     n = len(data)
    #     if n < 2:
    #         raise ValueError('variance requires at least two data points')
    #     ss = self._ss(data)
    #     pvar = ss / n  # the population variance
    #     return pvar ** 0.5

    def _log_missing_ids(self, missing_id, name_of_file_from_which_missing):
        if missing_id not in self.missing_assay_hash:
            self.missing_assay_hash[missing_id] = set()
        self.missing_assay_hash[missing_id].add(name_of_file_from_which_missing)
        # todo: remove the offending ids from the hash
        return

    def _fill_provenance_graph(self, limit):
        logger.info("Building graph ...")
        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)
        prov = Provenance()

        sexes = ['m', 'f']

        # loop through the hashmap
        for assay_id in self.assayhash:
            if assay_id not in self.missing_assay_hash:
                assay_label = self.assayhash[assay_id]['metadata']['assay_label']
                assay_type = self.assayhash[assay_id]['metadata']['assay_type']

                measurement_unit = self.assayhash[assay_id]['metadata']['measurement_unit']
                for sex in sexes:
                    index = 0

                    # First check that the assay_id is in the measurement hash
                    # and that it has a corresponding ontology mapping we are interested in
                    if len(self.assayhash[assay_id][sex]['strain_ids']) != 0:

                        for strain_id in self.assayhash[assay_id][sex]['strain_ids']:

                            ##############    ADD THE STRAIN AS GENOTYPE    #############
                            effective_genotype_id = self.make_id('-'.join((str(strain_id), sex)))
                            # todo: make this a method in geno that other classes can access
                            effective_genotype_label = str(strain_id) + ' (' + sex + ')'
                            geno.addGenotype(effective_genotype_id, effective_genotype_label,
                                             geno.genoparts['sex_qualified_genotype'])
                            geno.addParts(effective_genotype_id, effective_genotype_id,
                                          geno.object_properties['has_alternate_part'])

                            ##############    ADD THE TAXON AS CLASS    #############
                            taxon_id = 'NCBITaxon:10090'  # map to Mus musculus
                            gu.addClassToGraph(g, taxon_id, None)

                            ##############    BUILD THE G2P ASSOC    #############
                            # Phenotypes associations are made to limits strainid+gender+overall study

                            zscore = self.assayhash[assay_id][sex]['assay_zscores'][index]
                            if (zscore <= -self.stdevthreshold or zscore >= self.stdevthreshold):
                                logger.debug(str(assay_id) + sex + "\t" + str(zscore) + " significant!")

                                ont_term_ids = self.assayhash[assay_id]['metadata']['ont_terms']

                                for phenotype_id in ont_term_ids:
                                    eco_id = "ECO:0000059"  # experimental_phenotypic_evidence This was used in ZFIN

                                    # the association comes as a result of a g2p from a procedure in a pipeline at a center
                                    # and parameter tested

                                    assoc = G2PAssoc(self.name, effective_genotype_id, phenotype_id)
                                    assay_desc = self.assayhash[assay_id]['metadata']['description']

                                    measurement_value = self.assayhash[assay_id][sex]['assay_means'][index]
                                    significance = prov.get_zscore(zscore)

                                    logger.debug(
                                        "ASSAY: " + str(assay_id) +
                                        "\tLABEL:" + assay_label +
                                        "\tTYPE:" + assay_type +
                                        "\tDESC:" + assay_desc
                                    )

                                    effective_assay_type = \
                                        (self.make_id(str(assay_type.strip().replace(' ', '-')))) \
                                            if (assay_type is not None) \
                                            else None

                                prov.add_assay_to_graph(g, 'MPD-assay:' + str(assay_id), assay_label,
                                                        effective_assay_type, assay_desc)
                                prov.add_measurement_data(assay_id, measurement_unit, significance)
                                logger.debug(
                                    "AGENT ID: " + str(self.mgd_agent_id) + "\tAGENT LABEL: " + self.mgd_agent_label
                                    + "\tAGENT TYPE: " + self.mgd_agent_type
                                )
                                prov.add_agent_to_graph(g, self.mgd_agent_id, self.mgd_agent_label, self.mgd_agent_type)
                                prov.add_provenance_to_graph(g)

                                assoc.add_evidence(eco_id)
                                assoc.set_score(float(zscore))
                        else:
                            logger.debug(str(assay_id) + sex + "\t" + str(strain_id) + "\t" + str(
                                zscore) + " Not significant; skipping phenotype association.")
                    index += 1

                    # else:
                    #     logger.debug("skipping assay " + str(assay_id) + sex + " as it is missing any associated data.")

        return


def getTestSuite(self):
    import unittest
    from tests.test_mpd import MPDTestCase
    # TODO test genotypes

    test_suite = unittest.TestLoader().loadTestsFromTestCase(MPDTestCase)

    return test_suite
