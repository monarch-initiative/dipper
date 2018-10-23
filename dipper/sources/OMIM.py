import logging
import re
import json
import urllib
from urllib.error import HTTPError

from dipper.sources.Source import Source, USER_AGENT
from dipper.models.Model import Model
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper.models.Reference import Reference
from dipper import config
from dipper.utils.romanplus import romanNumeralPattern, fromRoman, toRoman

LOG = logging.getLogger(__name__)


# omimftp key EXPIRES
# get a new one here: https://omim.org/help/api

OMIMURL = 'https://data.omim.org/downloads/'

OMIMFTP = OMIMURL + config.get_config()['keys']['omim']

OMIMAPI = 'https://api.omim.org/api/entry?format=json&apiKey=' + \
    config.get_config()['keys']['omim'] + '&'


class OMIM(Source):
    """
    The only anonymously obtainable data from the ftp site is mim2gene.
    However, more detailed information is available via their API.
    So, we pull the omim identifiers from their ftp site,
    then query their API in batchs of 20.
    Their prescribed rate limits have been mecurial
     one per two seconds or  four per second,
     in  2017 November all mention of api rate limits have vanished
     (save 20 IDs per call if any include is used)

    Note this ingest requires an api Key which is not stored in the repo,
    but in a separate conf.yaml file.

    Processing this source serves two purposes:
    1.  the creation of the OMIM classes for merging into the disease ontology
    2.  add annotations such as disease-gene associations

    When creating the disease classes, we pull from their REST-api
    id/label/definition information.
    Additionally we pull the Orphanet and UMLS mappings
    (to make equivalent ids).
    We also pull the phenotypic series annotations as grouping classes.


    """

    files = {
        'all': {
            'file': 'mim2gene.txt',
            'url': 'https://omim.org/static/omim/data/mim2gene.txt',
            'clean': OMIMURL
        },
        'morbidmap': {
            'file': 'morbidmap.txt',
            'url': OMIMFTP + '/morbidmap.txt',
            'clean': OMIMURL
        },
        'phenotypicSeries': {
            'file': 'phenotypic_series_title_all.txt',
            'url': 'https://omim.org/phenotypicSeriesTitle/all?format=tsv',
            'headers': {'User-Agent': USER_AGENT},
            'clean': OMIMURL
        },
        'mimTitles': {
            'file': 'mimTitles.txt',
            'url':  OMIMFTP + '/mimTitles.txt',
            'headers': {'User-Agent': USER_AGENT},
            'clean': OMIMURL,
            'columns': (  # expected
                'Prefix',
                'Mim Number',
                'Preferred Title; symbol',
                'Alternative Title(s); symbol(s)',
                'Included Title(s); symbols',
            ),
        },
    }
    resources = {'test_ids': '../../resources/test_ids.yaml'}

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'omim',
            ingest_title='Online Mendelian Inheritance in Man',
            ingest_url='http://www.omim.org',
            # ingest_desc=None,
            license_url=None,
            data_rights='http://omim.org/help/agreement',
            # file_handle=None
        )

        self.omim_ncbigene_idmap = {}

        # check if config exists; if it doesn't, error out and let user know
        if 'keys' not in config.get_config() and \
                'omim' not in config.get_config()['keys']:
            LOG.error("not configured with API key.")

        all_test_ids = self.open_and_parse_yaml(self.resources['test_ids'])
        # integer portion of omim identifier
        self.test_ids = [x[5:] for x in all_test_ids['disease'] if x[:5] == 'OMIM:']

        self.omim_type = {}

        return

    def fetch(self, is_dl_forced=True):
        """
        Get the preconfigured static files.
        This DOES NOT fetch the individual records via REST...that is handled
        in the parsing function.  (To be refactored.)
        over riding Source.fetch()  calling Source.get_files()
        :param is_dl_forced:
        :return:

        """
        self.get_files(is_dl_forced)

        return

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %d rows", limit)

        LOG.info("Parsing files...")

        if self.testOnly:
            self.testMode = True

        self._process_all(limit)
        self._process_morbidmap(limit)
        self._process_phenotypicseries(limit)

        LOG.info("Done parsing.")

        return

    def _get_omim_ids(self):
        '''
            side effect:
                populate omim_type map from a omim number to an ontology term
                the ontology terms's labels as
                 -  'gene'
                    when they declare it as a gene

                -   'Phenotype'
                    Phenotype, molecular basis known

                -   'heritable_phenotypic_marker'
                    Phenotype or locus, molecular basis unknown

                -   'obsolete'
                    when Removed or moved to another entry

                -   'has_affected_feature'
                    "when declared as  "Gene and phenotype, combined"
                    hope being it could be detected and used as either

            :return a unique list of omim numbers
        '''

        omim_nums = set()        # all types
        line_counter = 0
        omimfile = '/'.join((self.rawdir, self.files['all']['file']))
        LOG.info("Obtaining OMIM record identifiers from: %s", omimfile)
        # TODO check to see if the file is there
        with open(omimfile, "r") as fh:
            # f.readline()            # copyright
            # line = f.readline()     # Generated: YYYY-MM-DD
            # f.readline()            # discription
            # f.readline()            # disclaimer
            # f.readline()            # column headers
            for line in fh:
                line_counter += 1
                if line[0] == '#':  # skip omments
                    continue

                (omim_num, mimtype, ncbigene, hgnc, ensembl) = line.split('\t')
                omim_nums.update({omim_num})
                if mimtype == 'gene':
                    self.omim_type[omim_num] = self.globaltt['gene']

                # Phenotype, molecular basis known
                elif mimtype == 'phenotype':
                    self.omim_type[omim_num] = self.globaltt['Phenotype']

                # Phenotype or locus, molecular basis unknown
                elif mimtype == 'predominantly phenotypes':
                    self.omim_type[omim_num] = self.globaltt[
                        'heritable_phenotypic_marker']  # ?

                # Removed or moved to another entry
                elif mimtype == 'moved/removed':
                    self.omim_type[omim_num] = self.globaltt['obsolete']

                # "Gene and phenotype, combined"  works as both/either.
                elif mimtype == 'gene/phenotype':
                    self.omim_type[omim_num] = self.globaltt['has_affected_feature']

        LOG.info("Done. found %d omim ids", len(omim_nums))

        return list(omim_nums)

    def process_entries(
            self, omimids, transform, included_fields=None, graph=None, limit=None,
            globaltt=None
    ):
        """
        Given a list of omim ids,
        this will use the omim API to fetch the entries, according to the
        ```included_fields``` passed as a parameter.
        If a transformation function is supplied,
        this will iterate over each entry,
        and either add the results to the supplied ```graph```
        or will return a set of processed entries that the calling function
        can further iterate.

        If no ```included_fields``` are provided, this will simply fetch
        the basic entry from omim,
        which includes an entry's:  prefix, mimNumber, status, and titles.

        :param omimids: the set of omim entry ids to fetch using their API
        :param transform: Function to transform each omim entry when looping
        :param included_fields: A set of what fields are required to retrieve
         from the API
        :param graph: the graph to add the transformed data into
        :return:
        """

        omimparams = {}

        # add the included_fields as parameters
        if included_fields is not None and len(included_fields) > 0:
            omimparams['include'] = ','.join(included_fields)

        processed_entries = list()

        # scrub any omim prefixes from the omimids before processing
        # cleanomimids = set()
        # for omimid in omimids:
        #    scrubbed = str(omimid).split(':')[-1]
        #    if re.match(r'^\d+$', str(scrubbed)):
        #        cleanomimids.update(scrubbed)
        # omimids = list(cleanomimids)

        cleanomimids = [o.split(':')[-1] for o in omimids]
        diff = set(omimids) - set(cleanomimids)
        if len(diff) > 0:
            LOG.warning('OMIM has %i dirty bits see"\n %s', len(diff), str(diff))
            omimids = cleanomimids
        else:
            cleanomimids = list()

        it = 0  # for counting

        # note that you can only do request batches of 20
        # see info about "Limits" at http://omim.org/help/api
        # TODO 2017 May seems a majority of many groups of 20
        # are producing python None for RDF triple Objects

        groupsize = 20
        if not self.testMode and limit is not None:
            # just in case the limit is larger than the number of records,
            maxit = limit
            if limit > len(omimids):
                maxit = len(omimids)
        else:
            maxit = len(omimids)

        while it < maxit:
            end = min((maxit, it + groupsize))
            # iterate through the omim ids list,
            # and fetch from the OMIM api in batches of 20

            if self.testMode:
                intersect = list(
                    set([str(i) for i in self.test_ids]) & set(omimids[it:end]))
                # some of the test ids are in the omimids
                if len(intersect) > 0:
                    LOG.info("found test ids: %s", intersect)
                    omimparams.update({'mimNumber': ','.join(intersect)})
                else:
                    it += groupsize
                    continue
            else:
                omimparams.update({'mimNumber': ','.join(omimids[it:end])})

            url = OMIMAPI + urllib.parse.urlencode(omimparams)
            LOG.info('fetching: %s', url)

            try:
                req = urllib.request.urlopen(url)
            except HTTPError as e:  # URLError?
                error_msg = e.read()
                if re.search(r'The API key: .* is invalid', str(error_msg)):
                    msg = "API Key not valid"
                    raise HTTPError(url, e.code, msg, e.hdrs, e.fp)
                else:
                    LOG.warning("url %s returned 404, skipping", url)
                    break

            resp = req.read().decode()
            it += groupsize

            myjson = json.loads(resp)
            # snag a copy, hopefully we will get the ftp dl soon
            with open('./raw/omim/_' + str(it) + '.json', 'w') as fp:
                json.dump(myjson, fp)

            entries = myjson['omim']['entryList']

            for e in entries:
                # apply the data transformation, and save it to the graph
                processed_entry = transform(e, graph, globaltt)
                if processed_entry is not None:
                    processed_entries.append(processed_entry)

                # ### end iterating over batch of entries
        return processed_entries

    def _process_all(self, limit):
        """
        This takes the list of omim identifiers from the omim.txt.Z file,
        and iteratively queries the omim api for the json-formatted data.
        This will create OMIM classes, with the label,
        definition, and some synonyms.
        If an entry is "removed",
            it is added as a deprecated class.
        If an entry is "moved",
            it is deprecated and consider annotations are added.

        Additionally, we extract:
        *phenotypicSeries ids as superclasses
        *equivalent ids for Orphanet and UMLS

        If set to testMode,
            it will write only those items in the test_ids to the testgraph.

        :param limit:
        :return:
        """

        omimids = self._get_omim_ids()
        LOG.info('Have %i omim numbers to fetch records from their API', len(omimids))
        LOG.info('Have %i omim types ', len(self.omim_type))

        if self.testMode:
            graph = self.testgraph
        else:
            graph = self.graph
        geno = Genotype(graph)
        model = Model(graph)
        tax_label = 'Homo sapiens'
        tax_id = self.globaltt[tax_label]

        # add genome and taxon
        geno.addGenome(tax_id, tax_label)   # tax label can get added elsewhere
        model.addClassToGraph(tax_id, None)   # label added elsewhere

        includes = set()
        includes.add('all')

        self.process_entries(
            omimids, self._transform_entry, includes, graph, limit, self.globaltt)

        return

    def _transform_entry(self, e, graph, globaltt):
        self.graph = graph
        model = Model(graph)
        geno = Genotype(graph)
        self.globaltt = globaltt
        tax_num = '9606'
        tax_label = 'Human'
        build_num = "GRCh38"
        build_id = "NCBIGenome:"+build_num

        # get the numbers, labels, and descriptions
        omim_num = str(e['entry']['mimNumber'])
        titles = e['entry']['titles']
        label = titles['preferredTitle']

        other_labels = []
        if 'alternativeTitles' in titles:
            other_labels += self._get_alt_labels(titles['alternativeTitles'])
        if 'includedTitles' in titles:
            other_labels += self._get_alt_labels(titles['includedTitles'])

        # add synonyms of alternate labels
        # preferredTitle": "PFEIFFER SYNDROME",
        # "alternativeTitles":
        #   "ACROCEPHALOSYNDACTYLY, TYPE V; ACS5;;\nACS V;;\nNOACK SYNDROME",
        # "includedTitles":
        #   "CRANIOFACIAL-SKELETAL-DERMATOLOGIC DYSPLASIA, INCLUDED"

        # remove the abbreviation (comes after the ;) from the preferredTitle,
        # and add it as a synonym
        abbrev = None
        if len(re.split(r';', label)) > 1:
            abbrev = (re.split(r';', label)[1].strip())
        newlabel = self._cleanup_label(label)

        omim_curie = 'OMIM:' + omim_num

        if e['entry']['status'] == 'removed':
            model.addDeprecatedClass(omim_curie)
        else:
            if omim_num in self.omim_type:
                omimtype = self.omim_type[omim_num]
            else:
                LOG.error('No type found for %s', omim_num)
                omimtype = None
            nodelabel = newlabel
            # this uses our cleaned-up label
            if omimtype == self.globaltt['heritable_phenotypic_marker']:
                if abbrev is not None:
                    nodelabel = abbrev
                # in this special case,
                # make it a disease by not declaring it as a gene/marker
                # ??? with none?

                model.addClassToGraph(
                    omim_curie, nodelabel, self.globaltt['disease or disorder'],
                    newlabel)
            elif omimtype == self.globaltt['gene']:
                if abbrev is not None:
                    nodelabel = abbrev
                model.addClassToGraph(omim_curie, nodelabel, omimtype, newlabel)
            else:
                model.addClassToGraph(omim_curie, newlabel, omimtype)

            # add the original screaming-caps OMIM label as a synonym
            model.addSynonym(omim_curie, label)

            # add the alternate labels and includes as synonyms
            for label in other_labels:
                model.addSynonym(omim_curie, label, model.globaltt['hasRelatedSynonym'])

            # KS: commenting out, we will get disease descriptions
            # from MONDO, and gene descriptions from the mygene API

            if abbrev is not None:
                model.addSynonym(
                    omim_curie, abbrev, model.globaltt['hasRelatedSynonym'])

            # if this is a genetic locus (but not sequenced)
            #   then add the chrom loc info
            # but add it to the ncbi gene identifier,
            # not to the omim id (we reserve the omim id to be the phenotype)
            feature_id = None
            feature_label = None
            if 'geneMapExists' in e['entry'] and e['entry']['geneMapExists']:
                genemap = e['entry']['geneMap']
                is_gene = False

                if omimtype == self.globaltt['heritable_phenotypic_marker']:
                    # get the ncbigene ids
                    ncbifeature = self._get_mapped_gene_ids(e['entry'], graph)
                    if len(ncbifeature) == 1:
                        feature_id = 'NCBIGene:' + str(ncbifeature[0])
                        # add this feature as a cause for the omim disease
                        # TODO SHOULD I EVEN DO THIS HERE?
                        assoc = G2PAssoc(graph, self.name, feature_id, omim_curie)
                        assoc.add_association_to_graph()

                    elif len(ncbifeature) > 1:
                        LOG.info(
                            "Its ambiguous when %s maps to >1 gene id: %s",
                            omim_curie, str(ncbifeature))
                    else:  # no ncbi feature, make an anonymous one
                        feature_id = self._make_anonymous_feature(omim_num)
                        feature_label = abbrev

                elif omimtype in {
                        self.globaltt['gene'], self.globaltt['has_affected_feature']}:
                    feature_id = omim_curie
                    is_gene = True
                else:
                    # 158900 falls into this category
                    feature_id = self._make_anonymous_feature(omim_num)
                    if abbrev is not None:
                        feature_label = abbrev
                    omimtype = self.globaltt['heritable_phenotypic_marker']

                if feature_id is not None:
                    if 'comments' in genemap:
                        # add a comment to this feature
                        comment = genemap['comments']
                        if comment.strip() != '':
                            model.addDescription(feature_id, comment)
                    if 'cytoLocation' in genemap:
                        cytoloc = genemap['cytoLocation']
                        # parse the cytoloc.
                        # add this omim thing as
                        # a subsequence of the cytofeature
                        # 18p11.3-p11.2
                        # FIXME
                        # add the other end of the range,
                        # but not sure how to do that
                        # not sure if saying subsequence of feature
                        # is the right relationship

                        f = Feature(graph, feature_id, feature_label, omimtype)
                        if 'chromosomeSymbol' in genemap:
                            chrom_num = str(genemap['chromosomeSymbol'])
                            chrom = makeChromID(chrom_num, tax_num, 'CHR')
                            geno.addChromosomeClass(
                                chrom_num, self.globaltt['Homo sapiens'], tax_label)

                            # add the positional information, if available
                            fstart = fend = -1
                            if 'chromosomeLocationStart' in genemap:
                                fstart = genemap['chromosomeLocationStart']
                            if 'chromosomeLocationEnd' in genemap:
                                fend = genemap['chromosomeLocationEnd']
                            if fstart >= 0:
                                # make the build-specific chromosome
                                chrom_in_build = makeChromID(
                                    chrom_num, build_num, 'MONARCH')
                                # then, add the chromosome instance
                                # (from the given build)
                                geno.addChromosomeInstance(
                                    chrom_num, build_id, build_num, chrom)
                                if omimtype == self.globaltt[
                                        'heritable_phenotypic_marker']:
                                    postypes = [self.globaltt['FuzzyPosition']]
                                else:
                                    postypes = None
                                # NOTE that no strand information
                                # is available in the API
                                f.addFeatureStartLocation(
                                    fstart, chrom_in_build, None, postypes)
                                if fend >= 0:
                                    f.addFeatureEndLocation(
                                        fend, chrom_in_build, None, postypes)
                                if fstart > fend:
                                    LOG.info(
                                        "start>end (%d>%d) for %s",
                                        fstart, fend, omim_curie)
                            # add the cytogenic location too
                            # for now, just take the first one
                            cytoloc = cytoloc.split('-')[0]
                            loc = makeChromID(cytoloc, tax_num, 'CHR')
                            model.addClassToGraph(loc, None)
                            f.addSubsequenceOfFeature(loc)
                            f.addFeatureToGraph(True, None, is_gene)

                # end adding causative genes/features

            # check if moved, if so,
            # make it deprecated and
            # replaced consider class to the other thing(s)
            # some entries have been moved to multiple other entries and
            # use the joining raw word "and"
            # 612479 is movedto:  "603075 and 603029"  OR
            # others use a comma-delimited list, like:
            # 610402 is movedto: "609122,300870"
            if e['entry']['status'] == 'moved':
                if re.search(r'and', str(e['entry']['movedTo'])):
                    # split the movedTo entry on 'and'
                    newids = re.split(r'and', str(e['entry']['movedTo']))
                elif len(str(e['entry']['movedTo']).split(',')) > 0:
                    # split on the comma
                    newids = str(e['entry']['movedTo']).split(',')
                else:
                    # make a list of one
                    newids = [str(e['entry']['movedTo'])]
                # cleanup whitespace and add OMIM prefix to numeric portion
                fixedids = []
                for i in newids:
                    fixedids.append('OMIM:'+i.strip())

                model.addDeprecatedClass(omim_curie, fixedids)

            self._get_phenotypicseries_parents(e['entry'], graph)
            self._get_mappedids(e['entry'], graph)
            self._get_mapped_gene_ids(e['entry'], graph)

            self._get_pubs(e['entry'], graph)

            self._get_process_allelic_variants(e['entry'], graph)  # temp gag

        return

    def _process_morbidmap(self, limit):
        """
        This will process the morbidmap file to get the links between
        omim genes and diseases. Here, we create anonymous nodes for some
        variant loci that are variants of the gene that causes the disease.
        Triples created:
        <some_anonymous_variant_locus>
            is_allele_of
                    <omim_gene_id>
        <some_anonymous_variant_locus> causes condition <omim_disease_id>
        <assoc> hasSubject <some_anonymous_variant_locus>
        <assoc> hasObject <omim_disease_id>
        <assoc> hasPredicate <causes condition>
        <assoc> DC:evidence <eco_id>
        :param limit:
        :return:
        """
        if self.testMode:
            graph = self.testgraph
        else:
            graph = self.graph
        line_counter = 0
        assoc_count = 0
        with open('/'.join((self.rawdir, self.files['morbidmap']['file']))) as fh:
            # Copyright
            # Generated: 2016-04-11
            # See end of file for additional documentation on specific fields
            # Phenotype	Gene Symbols	MIM Number	Cyto Location
            # since there are comments at the end of the file as well,
            # filter both header & footer as lines beginning with octothorp

            for line in fh:
                line = line.strip()
                if line.startswith(r'#') or line.startswith(r'\t'):
                    continue  # header/footer/ empty phenotype
                row = line.split('\t')
                if len(row) != 4:
                    LOG.warning("Expected 4 columns got %i columns.", len(row))
                    LOG.warning(row)
                    continue

                (disorder, gene_symbols, gene_num, loc) = row
                line_counter += 1

                # LOG.info("morbidmap disorder:  %s", disorder)  # too verbose

                # disorder = disorder label , number (mapping key)
                # 3-M syndrome 1, 273750 (3)|CUL7, 3M1|609577|6p21.1

                # but note that for those diseases where they are genomic loci
                # (not genes though), the omim id is only listed as the gene
                # Alopecia areata 1 (2)|AA1|104000|18p11.3-p11.2
                # when there's a gene and disease

                disorder_match = re.match(
                    r'(.*), (\d{6})\s*(?:\((\d+)\))?', disorder)
                nogene_match = re.match(r'(.*)\s+\((\d+)\)', disorder)

                if disorder_match is not None:
                    disorder_parts = disorder_match.groups()
                    (disorder_label, disorder_num, phene_key) = disorder_parts

                    if self.testMode and (
                            int(disorder_num) not in self.test_ids or
                            int(gene_num) not in self.test_ids):
                        continue
                    assoc_count += 1
                    gene_symbols = gene_symbols.split(', ')
                    gene_id = 'OMIM:' + str(gene_num)
                    self._make_pheno_assoc(
                        graph, gene_id, gene_symbols[0], disorder_num, disorder_label,
                        phene_key)
                elif nogene_match is not None:
                    # this is a case where the disorder
                    # a blended gene/phenotype
                    # we lookup the NCBIGene feature and make the association
                    (disorder_label, phene_key) = nogene_match.groups()
                    disorder_num = gene_num
                    # make what's in the gene column the disease
                    disorder_id = 'OMIM:' + str(disorder_num)
                    if self.testMode and int(disorder_num) not in self.test_ids:
                        continue
                    if disorder_id in self.omim_ncbigene_idmap:
                        # get the gene ids
                        gene_ids = self.omim_ncbigene_idmap[disorder_id]
                        if gene_ids is None:
                            continue
                        for gene_num in gene_ids:
                            # TODO add gene filter for testMode and NCBIGenes
                            gene_id = 'NCBIGene:' + str(gene_num).strip()
                            assoc_count += 1
                            self._make_pheno_assoc(
                                graph, gene_id, gene_symbols[0],
                                disorder_num, disorder_label, phene_key)
                    else:
                        # we can create an anonymous feature
                        # to house this thing for example, 158900
                        feature_id = self._make_anonymous_feature(gene_num)
                        assoc_count += 1
                        self._make_pheno_assoc(
                            graph, feature_id, gene_symbols[0], disorder_num,
                            disorder_label, phene_key)

                        LOG.info(
                            "We don't have an NCBIGene feature id to link %s with %s",
                            disorder_id, disorder_label)

                    if self.testMode and int(gene_num) not in self.test_ids:
                        continue

                else:
                    LOG.warning(
                        "There are misformatted rows %d:%s", line_counter, str(line))
                if not self.testMode and limit is not None and line_counter > limit:
                    break
            LOG.info("Added %d G2P associations", assoc_count)

        return

    @staticmethod
    def _make_anonymous_feature(omim_num):
        ''' more blank nodes '''
        return '_:feature' + omim_num

    def _make_pheno_assoc(
            self, graph, gene_id, gene_symbol, disorder_num, disorder_label, phene_key
    ):

        """
        From the docs:
        Brackets, "[ ]", indicate "nondiseases," mainly genetic variations
        that lead to apparently abnormal laboratory test values
        (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).

        Braces, "{ }", indicate mutations that contribute to susceptibility
        to multifactorial disorders (e.g., diabetes, asthma) or to
        susceptibility to infection (e.g., malaria).

        A question mark, "?", before the phenotype name indicates that the
        relationship between the phenotype and gene is provisional.
        More details about this relationship are provided in the comment
        field of the map and in the gene and phenotype OMIM entries.

        Phene key:
        The number in parentheses after the name of each disorder indicates
        the following:
          (1) the disorder was positioned by mapping of the wildtype gene;
          (2) the disease phenotype itself was mapped;
          (3) the molecular basis of the disorder is known;
          (4) the disorder is a chromosome deletion or duplication syndrome.

        reference: https://omim.org/help/faq#1_6

        :param graph: graph object of type dipper.graph.Graph
        :param gene_id: str, gene id as curie
        :param gene_symbol: str, symbol
        :param disorder_num: str, disorder id
        :param disorder_label: str, disorder label
        :param phene_key: int or str, 1-4, see docstring
        :return:
        """

        disorder_id = ':'.join(('OMIM', disorder_num))
        rel_label = 'causes condition'
        rel_id = self.globaltt[rel_label]

        if disorder_label.startswith('['):
            rel_id = self.globaltt['is marker for']
            # rel_label = 'is a marker for'
        elif disorder_label.startswith('{'):
            rel_id = self.globaltt['contributes to']
            # rel_label = 'contributes to'
        elif disorder_label.startswith('?'):
            # this is a questionable mapping!  skip?
            rel_id = self.globaltt['contributes to']

        assoc = G2PAssoc(graph, self.name, gene_id, disorder_id, rel_id)

        if phene_key is not None:
            evidence = self.resolve(phene_key, False)
            if evidence != phene_key:
                assoc.add_evidence(evidence)  # evidence is Found

        assoc.add_association_to_graph()

        return

    @staticmethod
    def _get_description(entry):
        """
        Get the description of the omim entity
        from the textSection called 'description'.
        Note that some of these descriptions have linebreaks.
        If printed in turtle syntax, they will appear to be triple-quoted.
        :param entry:
        :return:

        """
        description = None
        if entry is not None and 'textSectionList' in entry:
            textsectionlist = entry['textSectionList']
            for ts in textsectionlist:
                if ts['textSection']['textSectionName'] == 'description':
                    description = ts['textSection']['textSectionContent']
                    # there are internal references to OMIM identifiers in
                    # the description, I am formatting them in our style.
                    description = re.sub(r'{(\d+)}', r'OMIM:\1', description)

                    # TODO
                    # reformat the citations in the description with PMIDs
                    break

        return description

    def _get_process_allelic_variants(self, entry, graph):
        model = Model(graph)
        reference = Reference(graph)
        geno = Genotype(graph)
        if entry is not None:
            # to hold the entry-specific publication mentions
            # for the allelic variants
            publist = {}
            entry_num = entry['mimNumber']

            # process the ref list just to get the pmids
            ref_to_pmid = self._get_pubs(entry, graph)

            if 'allelicVariantList' in entry:
                allelicVariantList = entry['allelicVariantList']
                for al in allelicVariantList:
                    al_num = al['allelicVariant']['number']
                    al_id = 'OMIM:'+str(entry_num)+'.'+str(al_num).zfill(4)
                    al_label = None
                    al_description = None
                    if al['allelicVariant']['status'] == 'live':
                        publist[al_id] = set()
                        if 'mutations' in al['allelicVariant']:
                            al_label = al['allelicVariant']['mutations']
                        if 'text' in al['allelicVariant']:
                            al_description = al['allelicVariant']['text']
                            m = re.findall(r'\{(\d+)\:', al_description)
                            publist[al_id] = set(m)
                        geno.addAllele(
                            al_id, al_label, self.globaltt['variant_locus'],
                            al_description)
                        geno.addAlleleOfGene(
                            al_id, 'OMIM:' + str(entry_num),
                            self.globaltt['is_allele_of'])
                        for r in publist[al_id]:
                            pmid = ref_to_pmid[int(r)]
                            graph.addTriple(pmid, self.globaltt['is_about'], al_id)
                        # look up the pubmed id in the list of references
                        if 'dbSnps' in al['allelicVariant']:
                            dbsnp_ids = re.split(r',', al['allelicVariant']['dbSnps'])
                            for dnum in dbsnp_ids:
                                did = 'dbSNP:'+dnum.strip()
                                model.addIndividualToGraph(did, None)
                                model.addSameIndividual(al_id, did)

                        # Note that RCVs are variant to disease associations
                        # in ClinVar, rather than variant entries
                        # so we make these xrefs instead of equivalents
                        if 'clinvarAccessions' in al['allelicVariant']:
                            # clinvarAccessions triple semicolon delimited
                            # each >1 like RCV000020059;;;
                            rcv_ids = re.split(
                                r';;;', al['allelicVariant']['clinvarAccessions'])
                            rcv_ids = [
                                (re.match(r'(RCV\d+);*', r)).group(1) for r in rcv_ids]
                            for rnum in rcv_ids:
                                rid = 'ClinVar:' + rnum
                                model.addXref(al_id, rid)
                        reference.addPage(
                            al_id, "http://omim.org/entry/" +
                            str(entry_num)+"#" + str(al_num).zfill(4))
                    elif re.search(
                            r'moved', al['allelicVariant']['status']):
                        # for both 'moved' and 'removed'
                        moved_ids = None
                        if 'movedTo' in al['allelicVariant']:
                            moved_id = 'OMIM:'+al['allelicVariant']['movedTo']
                            moved_ids = [moved_id]
                        model.addDeprecatedIndividual(al_id, moved_ids)
                    else:
                        LOG.error(
                            'Uncaught alleleic variant status %s',
                            al['allelicVariant']['status'])
                # end loop allelicVariantList
        return

    @staticmethod
    def _cleanup_label(label):
        """
        Reformat the ALL CAPS OMIM labels to something more pleasant to read.
        This will:
        1.  remove the abbreviation suffixes
        2.  convert the roman numerals to integer numbers
        3.  make the text title case,
            except for suplied conjunctions/prepositions/articles
        :param label:
        :return:
        """
        conjunctions = ['and', 'but', 'yet', 'for', 'nor', 'so']
        little_preps = [
            'at', 'by', 'in', 'of', 'on', 'to', 'up', 'as', 'it', 'or']
        articles = ['a', 'an', 'the']

        # remove the abbreviation
        l = re.split(r';', label)[0]

        fixedwords = []
        i = 0
        for w in l.split():
            i += 1
            # convert the roman numerals to numbers,
            # but assume that the first word is not
            # a roman numeral (this permits things like "X inactivation"
            if i > 1 and re.match(romanNumeralPattern, w):
                n = fromRoman(w)
                # make the assumption that the number of syndromes are <100
                # this allows me to retain "SYNDROME C"
                # and not convert it to "SYNDROME 100"
                if 0 < n < 100:
                    # get the non-roman suffix, if present.
                    # for example, IIIB or IVA
                    suffix = w.replace(toRoman(n), '', 1)
                    fixed = ''.join((str(n), suffix))
                    w = fixed

            # capitalize first letter
            w = w.title()

            # replace interior conjunctions, prepositions,
            # and articles with lowercase
            if w.lower() in (conjunctions+little_preps+articles) and i != 1:
                w = w.lower()

            fixedwords.append(w)

        l = ' '.join(fixedwords)
        # print (label,'-->',l)
        return l

    def _process_phenotypicseries(self, limit):
        """
        Creates classes from the OMIM phenotypic series list.
        These are grouping classes to hook the more granular OMIM diseases.
        :param limit:
        :return:

        """
        if self.testMode:
            graph = self.testgraph
        else:
            graph = self.graph
        LOG.info("getting phenotypic series titles")
        model = Model(graph)
        line_counter = 0
        with open('/'.join((self.rawdir,
                            self.files['phenotypicSeries']['file']))) as fh:

            for i in range(5):  # yep, header blurb is not commented
                fh.readline()
                line_counter += 1

            for line in fh:
                line = line.strip()
                line_counter += 1
                if re.match(r'^\w*$', line) or line[0] == '#':
                    # skip blank lines and comments,
                    continue
                row = line.split('\t')
                if len(row) < 2:
                    LOG.warning(
                        'Unexpected input on line: %i  got: %s', line_counter, line)
                    continue
                ps_label = row[0]
                ps_num = row[1]

                omim_curie = 'OMIM:' + ps_num
                model.addClassToGraph(omim_curie, ps_label)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        return

    @staticmethod
    def _get_phenotypicseries_parents(entry, graph):
        """
        Extract the phenotypic series parent relationship out of the entry
        :param entry:
        :return:
        """
        model = Model(graph)
        omim_num = str(entry['mimNumber'])
        omim_curie = 'OMIM:' + omim_num
        # the phenotypic series mappings
        serieslist = []
        if 'phenotypicSeriesExists' in entry:
            if entry['phenotypicSeriesExists'] is True:
                if 'phenotypeMapList' in entry:
                    phenolist = entry['phenotypeMapList']
                    for p in phenolist:
                        serieslist.append(
                            p['phenotypeMap']['phenotypicSeriesNumber'])
                if 'geneMap' in entry and 'phenotypeMapList' in entry['geneMap']:
                    phenolist = entry['geneMap']['phenotypeMapList']
                    for p in phenolist:
                        if 'phenotypicSeriesNumber' in p['phenotypeMap']:
                            serieslist.append(
                                p['phenotypeMap']['phenotypicSeriesNumber'])
        # add this entry as a subclass of the series entry
        for ser in serieslist:
            series_id = 'OMIM:' + ser
            model.addClassToGraph(series_id, None)
            model.addSubClass(omim_curie, series_id)

        return

    # TODO PYLINT Method could be a function
    def _get_mappedids(self, entry, graph):
        """
        Extract the Orphanet and UMLS ids as equivalences from the entry
        :param entry:
        :return:
        """
        model = Model(graph)
        omim_num = str(entry['mimNumber'])
        omim_curie = 'OMIM:' + omim_num
        orpha_mappings = []
        if 'externalLinks' in entry:
            links = entry['externalLinks']
            if 'orphanetDiseases' in links:
                # triple semi-colon delimited list of
                # double semi-colon delimited orphanet ID/disease pairs
                # 2970;;566;;Prune belly syndrome
                items = links['orphanetDiseases'].split(';;;')
                for item in items:
                    (orpha_num, internal_num, orpha_label) = item.split(';;')
                    orpha_curie = 'ORPHA:' + orpha_num.strip()
                    orpha_mappings.append(orpha_curie)
                    model.addClassToGraph(orpha_curie, orpha_label.strip())
                    model.addXref(omim_curie, orpha_curie)

            if 'umlsIDs' in links:
                umls_mappings = links['umlsIDs'].split(',')
                for umls in umls_mappings:
                    umls_curie = 'UMLS:' + umls
                    model.addClassToGraph(umls_curie, None)
                    model.addXref(omim_curie, umls_curie)
        return

    def _get_mapped_gene_ids(self, entry, graph):

        gene_ids = []
        model = Model(graph)
        omim_num = str(entry['mimNumber'])
        omim_curie = 'OMIM:' + omim_num
        if 'externalLinks' in entry:
            links = entry['externalLinks']
            omimtype = self.omim_type[omim_num]
            if 'geneIDs' in links:
                entrez_mappings = links['geneIDs']
                gene_ids = entrez_mappings.split(',')
                self.omim_ncbigene_idmap[omim_curie] = gene_ids
                if omimtype == self.globaltt['gene']:
                    for ncbi in gene_ids:
                        model.addEquivalentClass(omim_curie, 'NCBIGene:' + str(ncbi))

        return gene_ids

    def _get_alt_labels(self, titles):
        """
        From a string of delimited titles, make an array.
        This assumes that the titles are double-semicolon (';;') delimited.
        This will additionally pass each through the _cleanup_label method to
        convert the screaming ALL CAPS to something more pleasant to read.
        :param titles:
        :return: an array of cleaned-up labels
        """

        labels = []
        # "alternativeTitles": "
        #   ACROCEPHALOSYNDACTYLY, TYPE V; ACS5;;\nACS V;;\nNOACK SYNDROME",
        # "includedTitles":
        #   "CRANIOFACIAL-SKELETAL-DERMATOLOGIC DYSPLASIA, INCLUDED"

        for title in titles.split(';;'):
            # remove ', included', if present
            label = re.sub(r',\s*INCLUDED', '', title.strip(), re.IGNORECASE)
            label = self._cleanup_label(label)
            labels.append(label)

        return labels

    def _get_pubs(self, entry, graph):
        """
        Extract mentioned publications from the reference list
        :param entry:
        :return:
        """

        ref_to_pmid = {}
        entry_num = entry['mimNumber']
        if 'referenceList' in entry:
            reflist = entry['referenceList']
            for rlst in reflist:
                if 'pubmedID' in rlst['reference']:
                    pub_id = 'PMID:' + str(rlst['reference']['pubmedID'])
                    ref = Reference(
                        graph, pub_id, self.globaltt['journal article'])
                else:
                    # make blank node for internal reference
                    pub_id = '_:OMIM' + str(entry_num) + 'ref' + str(
                        rlst['reference']['referenceNumber'])

                    ref = Reference(graph, pub_id)
                    title = author_list = source = citation = None
                    if 'title' in rlst['reference']:
                        title = rlst['reference']['title']
                        ref.setTitle(title)
                    if 'authors' in rlst['reference']:
                        author_list = rlst['reference']['authors']
                        ref.setAuthorList(author_list)
                        citation = re.split(r'\.\,', author_list)[0] + ' et al'
                    if 'source' in rlst['reference']:
                        source = rlst['reference']['source']
                    citation = '; '.join(
                        [tok for tok in [citation, title, source] if tok is not None])
                    ref.setShortCitation(citation)
                ref.addRefToGraph()
                ref_to_pmid[rlst['reference']['referenceNumber']] = pub_id

                # add is_about for the pub
                omim_id = 'OMIM:' + str(entry_num)
                graph.addTriple(omim_id, self.globaltt['mentions'], pub_id)

        return ref_to_pmid

    @staticmethod
    def _get_omimtype(entry, globaltt):
        """
        (note: there is anlaternative using mimTitle in omia)


        Here, we look at the omim 'prefix' to help to type the entry.
        For now, we only classify omim entries as genes;
        the rest we leave alone.
        :param entry:
        :return:
        """

        # An asterisk (*) before an entry number indicates a gene.
        # A number symbol (#) before an entry number indicates
        # that it is a descriptive entry, usually of a phenotype,
        # and does not represent a unique locus.
        # The reason for the use of the number symbol
        # is given in the first paragraph of the entry.
        # Discussion of any gene(s) related to the phenotype resides in
        # another entry(ies) as described in the first paragraph.
        #
        # A plus sign (+) before an entry number indicates that the
        # entry contains the description of a gene of
        # known sequence and a phenotype.
        #
        # A percent sign (%) before an entry number indicates that the
        # entry describes a confirmed mendelian phenotype or phenotypic locus
        # for which the underlying molecular basis is not known.
        #
        # No symbol before an entry number generally indicates a
        # description of a phenotype for which the mendelian basis,
        # although suspected, has not been clearly established
        # or that the separateness of this phenotype
        # from that in another entry is unclear.
        #
        # A caret (^) before an entry number means the
        # entry no longer exists because it was removed from the database
        # or moved to another entry as indicated.
        prefix = None
        type_id = None
        if 'prefix' in entry:
            prefix = entry['prefix']

        if prefix == '*':
            # gene, may not have a known sequence or a phenotype
            # note that some genes are also phenotypes,
            # even in this class, like 102480
            # examples: 102560,102480,100678,102750
            type_id = globaltt['gene']
        elif prefix == '#':
            # phenotype/disease -- indicate that here?
            # examples: 104200,105400,114480,115300,121900
            # type_id = globaltt['Phenotype']  # 'UPHENO_0001001' # species agnostic
            # type_id = globaltt['human phenotypic abnormality']
            pass

        elif prefix == '+':
            # gene of known sequence and has a phenotype
            # examples: 107670,110600,126453
            type_id = globaltt['gene']  # doublecheck this
        elif prefix == '%':
            # this is a disease (with a known locus).
            # examples include:  102150,104000,107200,100070
            type_id = globaltt['heritable_phenotypic_marker']
        elif prefix == '':
            # this is probably just a phenotype
            pass

        return type_id

    # def getTestSuite(self):
    #   ''' this should find a home under /test , if it is needed'''
    #        import unittest
    #        # TODO PYLINT  Unable to import 'tests.test_omim'
    #   from tests.test_omim import OMIMTestCase
    #
    #   test_suite = unittest.TestLoader().loadTestsFromTestCase(OMIMTestCase)
    #   return test_suite


def get_omim_id_from_entry(entry):
    if entry is not None and 'mimNumber' in entry:
        omimid = 'OMIM:' + str(entry['mimNumber'])
    else:
        omimid = None
    return omimid
