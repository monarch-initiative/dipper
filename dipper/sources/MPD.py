import csv
import gzip
import re
import logging
import os

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

        'measurements': {'file': 'measurements.csv',
                         'url': 'http://phenome.jax.org/download/measurements.csv'},

        'strainmeans': {'file': 'strainmeans_cleanzipped.gz',
                        'url': 'http://phenome.jax.org/download/strainmeans_cleanzipped.gz'},

        'mpd_datasets_metadata': {'file': 'mpd_datasets_metadata.xml',
                                  'url': 'http://phenome.jax.org/download/mpd_datasets_metadata.xml'},

    }

    # # TODO move these into the conf.json
    # # the following are strain ids for testing
    test_ids = ["MPD:2", "MPD:3", "MPD:5", "MPD:6", "MPD:9", "MPD:11", "MPD:18", "MPD:20", "MPD:24", "MPD:28", "MPD:30",
                "MPD:33", "MPD:34", "MPD:36", "MPD:37", "MPD:39", "MPD:40", "MPD:42", "MPD:47", "MPD:66", "MPD:68",
                "MPD:71", "MPD:75", "MPD:78", "MPD:122", "MPD:169", "MPD:438", "MPD:457", "MPD:473", "MPD:481",
                "MPD:759", "MPD:766", "MPD:770", "MPD:849", "MPD:857", "MPD:955", "MPD:964", "MPD:988", "MPD:1005",
                "MPD:1017", "MPD:1204", "MPD:1233", "MPD:1235", "MPD:1236", "MPD:1237"]

    def __init__(self):
        Source.__init__(self, 'mpd')
        # @N, not sure if this step is required
        self.namespaces.update(curie_map.get())

        # update the dataset object with details about this resource
        # @N: Note that there is no license as far as I can tell
        self.dataset = Dataset('mpd', 'MPD', 'http://phenome.jax.org', None, None)

        # TODO add a citation for mpd dataset as a whole
        # :mpd cito:citesAsAuthority PMID:15619963

        # so that we don't have to deal with BNodes, we will create hash lookups for the internal identifiers
        # the hash will hold the type-specific-object-keys to MPD public identifiers.  then, subsequent
        # views of the table will lookup the identifiers in the hash.  this allows us to do the 'joining' on the
        # fly
        self.idhash = {'strainid': {}, 'sex': {}, 'measurements': {}}
        self.label_hash = {}  # use this to store internally generated labels for various features

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

        # the following will provide us the hash-lookups
        # These must be processed in a specific order
        self._process_ontology_mappings(limit)
        self._process_measurements(limit)
        self._process_strainmeans(limit)

        # The following will use the hash populated above to lookup the ids when filling in the graph


        logger.info("Finished parsing.")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return

    def _process_ontology_mappings(self, limit):

        line_counter = 0

        logger.info("Processing ontology mappings...")
        raw = '/'.join((self.rawdir, 'ontology_mappings.csv'))

        with open(raw, 'r') as f:
            reader = csv.reader(f)
            f.readline()  # read the header row; skip
            for row in reader:
                line_counter += 1
                (measnum, ont_term, descrip) = row

                if re.match('MP', ont_term) or re.match('VT', ont_term):
                    # add the mapping denovo
                    if measnum not in self.idhash['measurements']:
                        # create space for the measurement description to go later
                        self.idhash['measurements'][measnum] = [[ont_term],'placeholder_description']
                    # add the new mapping to the existing list of mappings
                    else:
                        self.idhash['measurements'][measnum][0] += [ont_term]

                    # the assays will be added to the graph later as part of provenance
        return

    def _process_measurements(self, limit):
        line_counter = 0

        logger.info("Processing measurements ...")
        raw = '/'.join((self.rawdir, 'measurements.csv'))

        with open(raw, 'r') as f:
            reader = csv.reader(f)
            f.readline()  # read the header row; skip
            for row in reader:
                line_counter += 1
                (measnum, projsym, varname, descrip, units, cat1, cat2, cat3, intervention, intparm, appmeth, panelsym,
                 datatype, sextested, nstrainstested, ageweeks) = row

                if measnum in self.idhash['measurements']:
                    sextested = 'female' if (sextested is 'fm') else 'male'
                    units = (units)
                    description = "This is an assay of " + descrip + " as a " + datatype + " measured in " + units
                    if intervention is not None and intervention != "":
                        description += " in response to " + intervention

                    description += ". This represents the " + intparm + " arm of the overall experiment entitled " \
                                   + projsym + " using materials and methods that included " + appmeth + "."
                    description += "It was conducted in " + sextested + " mice at " + ageweeks + " of age in" \
                                   + nstrainstested + " different mouse strains. "
                    description += "Keywords: " + cat1 + \
                                   ((", " + cat2) if cat2 is not "" else "") + \
                                   ((", " + cat3) if cat3 is not "" else "") + "."

                    self.idhash['measurements'][measnum][1] = description
        return

    def normalise_units(units):
        # todo:

        # &deg;C
        # &micro;g/dL
        # &micro;g/g
        # &micro;g/mL
        # &micro;L
        # &micro;L/&micro;L
        # &micro;L/cm
        # &micro;L/g
        # &micro;m
        # &micro;m<sup>2</sup>
        # &micro;m<sup>2</sup>/n
        # &micro;mol/g
        # &micro;mol/L
        # &micro;mol/mg
        # &micro;mol/min/hPa
        # &micro;mol/min/hPa/mL
        # &micro;mol/mL/h
        # &micro;s
        # &micro;V
        # %
        # 1/cm
        # amplitude
        # AUC
        # cal
        # cm
        # cm/mL
        # cm/mL/s
        # cm/s
        # cm/s<sup>2</sup>
        # cm/sec
        # cm&sdot;s
        # cm<sup>2</sup>
        # cm<sup>3</sup>
        # d
        # daPa
        # dB
        # deg
        # deg/cm
        # deg/s
        # designation
        # fL
        # frequency rate
        # g
        # g/cm<sup>2</sup>
        # g/cm<sup>3</sup>
        # g/dL
        # g/kg
        # g/wk
        # h
        # Hz
        # IU/L
        # kcal
        # kcal/d
        # kcal/g/d
        # kcal/h
        # kg/m<sup>2</sup>
        # kg&sdot;m
        # km
        # km/d
        # L/g
        # log(&micro;g/dL)
        # m
        # m/d
        # m/min
        # mA
        # mEq/L
        # mg
        # mg/cm<sup>2</sup>
        # mg/dL
        # mg/g
        # mg/kg
        # mg/mL
        # min
        # min/d
        # mL
        # mL/kg
        # mL/kg/h
        # mL/kg/s
        # mL/min
        # mm
        # mm/mL
        # mm/s
        # mm<sup>2</sup>
        # mm<sup>2</sup>/g
        # mm<sup>3</sup>
        # mmHg
        # mmol/kg
        # mmol/L
        # mmol/L&sdot;min
        # mmol&sdot;h
        # mN
        # mol/m<sup>3</sup>
        # ms
        # ms*mV
        # mV
        # n
        # n/&micro;L
        # n/&micro;m<sup>2</sup>
        # n/cm
        # n/d
        # n/g
        # n/L
        # n/min
        # n/mL
        # N/mm
        # N/mm<sup>2</sup>
        # n/mm<sup>3</sup>
        # n/s
        # n/wks
        # N&sdot;mm
        # ng/g
        # ng/mL
        # nm
        # nmol/g
        # nmol/mg
        # nmol/mL
        # pg
        # pg/mL
        # pH
        # pmol/g
        # pmol/L
        # pmol/mg
        # pmol/mL
        # ratio
        # s
        # score
        # slope
        # U/g
        # U/L
        # wks

        return units


    def _process_strainmeans(self, limit):
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

        # make a pass through the table first, to create the mapping between the external and internal identifiers
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph

        geno = Genotype(g)
        logger.info("mapping strains to internal identifiers")
        raw = '/'.join((self.rawdir, 'strainmeans.csv'))
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (measnum, varname, strain, strainid, sex, mean, nmice, sd, sem, cv, minval, maxval, logmean, logsd,
                 zscore, logzscore, measinfo, mp_mapping, ma_mapping) = line.split('\t')

                if self.testMode is True:
                    if int(strainid) not in self.test_keys.get('strain'):
                        continue

                # we are only interested in values of significance; 2 is a weak cutoff but <.1% were significant at 3
                # standard deviations above or below mean of means. Todo: calculate these properly based on the raw data
                if zscore >= 2 or zscore < -2:
                    continue

                # # get the hashmap of the identifiers
                # self.idhash['strain'][strainid] = accid
                # gu.addIndividualToGraph(g, accid, None, geno.genoparts['intrinsic_genotype'])

                # pass through the file again, and make the equivalence statements to a subset of the idspaces
                logger.info("mapping strain equivalent identifiers")
                line_counter = 0
                with open(raw, 'r') as f:
                    f.readline()  # read the header row; skip
                    for line in f:
                        line_counter += 1
                        (accession_key, accid, prefixpart, numericpart, logicaldb_key, object_key, mgitype_key, private,
                         preferred, createdby_key, modifiedby_key,
                         creation_date, modification_date, logicaldb) = line.split('\t')

                        if self.testMode is True:
                            if int(object_key) not in self.test_keys.get('strain'):
                                continue

                        mgiid = self.idhash['strain'].get(object_key)
                        if mgiid is None:
                            # presumably we've already added the relevant MGI ids already, so skip those that we can't find
                            # logger.info("can't find mgiid for %s",object_key)
                            continue
                        strain_id = None
                        if preferred == '1':  # what does it mean if it's 0?
                            if logicaldb_key == '22':  # JAX
                                # scrub out the backticks from accids
                                # TODO notify the source upstream
                                accid = re.sub('`', '', accid)
                                strain_id = 'JAX:' + accid
                                # TODO get non-preferred ids==deprecated?

                        if strain_id is not None:
                            gu.addIndividualToGraph(g, strain_id, None, geno.genoparts['intrinsic_genotype'])
                            gu.addSameIndividual(g, mgiid, strain_id)

                        if not self.testMode and limit is not None and line_counter > limit:
                            break

        return


def getTestSuite(self):
    import unittest
    from tests.test_mpd import MPDTestCase
    # TODO test genotypes

    test_suite = unittest.TestLoader().loadTestsFromTestCase(MPDTestCase)

    return test_suite
