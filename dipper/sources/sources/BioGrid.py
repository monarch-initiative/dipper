import os
import logging
import re

from datetime import datetime
from stat import ST_CTIME
from zipfile import ZipFile

from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.assoc.InteractionAssoc import InteractionAssoc

__author__ = 'nicole'

LOG = logging.getLogger(__name__)
BGDL = 'http://thebiogrid.org/downloads/archives/Latest%20Release'


class BioGrid(Source):
    """
    Biogrid interaction data

    """
    # TODO write up class summary for docstring

    files = {
        'interactions': {
            'file': 'interactions.mitab.zip',
            'url': BGDL + '/BIOGRID-ALL-LATEST.mitab.zip'},
        'identifiers':  {
            'file': 'identifiers.tab.zip',
            'url': BGDL + '/BIOGRID-IDENTIFIERS-LATEST.tab.zip'}
    }

    # biogrid-specific identifiers for use in subsetting identifier mapping
    biogrid_ids = [
        106638, 107308, 107506, 107674, 107675, 108277, 108506, 108767, 108814,
        108899, 110308, 110364, 110678, 111642, 112300, 112365, 112771, 112898,
        199832, 203220, 247276, 120150, 120160, 124085]

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'biogrid',
            ingest_title='Biological General Repository for Interaction Datasets',
            ingest_url='http://thebiogrid.org',
            license_url='https://wiki.thebiogrid.org/doku.php/terms_and_conditions'
            # data_rights=None,
            # file_handle=None
        )

        self.tax_ids = tax_ids
        # Defaults
        # our favorite animals
        # taxids = [9606,10090,10116,7227,7955,6239,8355]
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]

        if 'gene' not in self.all_test_ids:
            LOG.warning("not configured with gene test ids.")
        else:
            self.test_ids = self.all_test_ids['gene']

        # data-source specific warnings
        # (will be removed when issues are cleared)
        LOG.warning(
            "several MI experimental codes do not exactly map to ECO; "
            "using approximations.")
        return

    def fetch(self, is_dl_forced=False):
        """

        :param is_dl_forced:
        :return:  None
        """

        self.get_files(is_dl_forced)

        # the version number is encoded in the filename in the zip.
        # for example, the interactions file may unzip to
        # BIOGRID-ALL-3.2.119.mitab.txt, where the version number is 3.2.119
        f = '/'.join((self.rawdir, self.files['interactions']['file']))
        st = os.stat(f)
        filedate = datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")
        with ZipFile(f, 'r') as myzip:
            flist = myzip.namelist()
            # assume that the first entry is the item
            fname = flist[0]
            # get the version from the filename
            version = re.match(
                r'BIOGRID-ALL-(\d+\.\d+\.\d+)\.mitab.txt', fname)
        myzip.close()

        self.dataset.setVersion(filedate, str(version.groups()[0]))

        return

    def parse(self, limit=None):
        """

        :param limit:
        :return:

        """
        if self.testOnly:
            self.testMode = True

        self._get_interactions(limit)
        self._get_identifiers(limit)

        LOG.info("Loaded %d test graph nodes", len(self.testgraph))
        LOG.info("Loaded %d full graph nodes", len(self.graph))

        return

    def _get_interactions(self, limit):
        LOG.info("getting interactions")
        line_counter = 0
        f = '/'.join((self.rawdir, self.files['interactions']['file']))
        myzip = ZipFile(f, 'r')
        # assume that the first entry is the item
        fname = myzip.namelist()[0]
        matchcounter = 0

        with myzip.open(fname, 'r') as csvfile:
            for line in csvfile:
                # skip comment lines
                if re.match(r'^#', line.decode()):
                    LOG.debug("Skipping header line")
                    continue
                line_counter += 1
                line = line.decode().strip()
                # print(line)
                (interactor_a, interactor_b, alt_ids_a, alt_ids_b, aliases_a,
                 aliases_b, detection_method, pub_author, pub_id, taxid_a,
                 taxid_b, interaction_type, source_db, interaction_id,
                 confidence_val) = line.split('\t')
                taxid_a = taxid_a.rstrip()
                taxid_b = taxid_b.rstrip()

                # get the actual gene ids,
                # typically formated like: gene/locuslink:351|BIOGRID:106848
                gene_a_num = re.search(r'locuslink\:(\d+)\|?', interactor_a).groups()[0]
                gene_b_num = re.search(r'locuslink\:(\d+)\|?', interactor_b).groups()[0]

                if self.testMode:
                    graph = self.testgraph
                    # skip any genes that don't match our test set
                    if (int(gene_a_num) not in self.test_ids) or\
                            (int(gene_b_num) not in self.test_ids):
                        continue
                else:
                    graph = self.graph
                    # when not in test mode, filter by taxon
                    if int(taxid_a.split(':')[-1]) not in self.tax_ids or \
                            int(taxid_b.split(':')[-1]) not in self.tax_ids:
                        continue
                    else:
                        matchcounter += 1

                gene_a = 'NCBIGene:'+gene_a_num
                gene_b = 'NCBIGene:'+gene_b_num

                # get the interaction type
                # psi-mi:"MI:0407"(direct interaction)
                int_type = re.search(r'MI:\d+', interaction_type).group()
                rel = self.resolve(int_type, False)
                if rel == int_type:
                    rel = self.globaltt['interacts with']

                # scrub pubmed-->PMID prefix
                pub_id = re.sub(r'pubmed', 'PMID', pub_id)
                # remove bogus whitespace
                pub_id = pub_id.strip()

                # get the method, and convert to evidence code
                det_code = re.search(r'MI:\d+', detection_method).group()
                evidence = self.resolve(det_code, False)
                if evidence == det_code:
                    evidence = self.globaltt["experimental evidence"]

                # note that the interaction_id is some kind of internal biogrid
                # identifier that does not map to a public URI.
                # we will construct a monarch identifier from this

                assoc = InteractionAssoc(graph, self.name, gene_a, gene_b, rel)
                assoc.add_evidence(evidence)
                assoc.add_source(pub_id)
                assoc.add_association_to_graph()

                if not self.testMode and (
                        limit is not None and line_counter > limit):
                    break

        myzip.close()

        return

    def _get_identifiers(self, limit):
        """
        This will process the id mapping file provided by Biogrid.
        The file has a very large header, which we scan past,
        then pull the identifiers, and make equivalence axioms

        :param limit:
        :return:

        """

        LOG.info("getting identifier mapping")
        line_counter = 0
        f = '/'.join((self.rawdir, self.files['identifiers']['file']))
        myzip = ZipFile(f, 'r')
        # assume that the first entry is the item
        fname = myzip.namelist()[0]
        foundheader = False

        # TODO align this species filter with the one above
        # speciesfilters = 'Homo sapiens,Mus musculus,Drosophila melanogaster,
        # Danio rerio, Caenorhabditis elegans,Xenopus laevis'.split(',')

        speciesfilters = 'Homo sapiens,Mus musculus'.split(',')
        with myzip.open(fname, 'r') as csvfile:
            for line in csvfile:
                # skip header lines
                if not foundheader:
                    if re.match(r'BIOGRID_ID', line.decode()):
                        foundheader = True
                    continue

                line = line.decode().strip()
                # BIOGRID_ID
                # IDENTIFIER_VALUE
                # IDENTIFIER_TYPE
                # ORGANISM_OFFICIAL_NAME
                # 1	814566	ENTREZ_GENE	Arabidopsis thaliana
                (biogrid_num, id_num, id_type,
                 organism_label) = line.split('\t')

                if self.testMode:
                    graph = self.testgraph
                    # skip any genes that don't match our test set
                    if int(biogrid_num) not in self.biogrid_ids:
                        continue
                else:
                    graph = self.graph

                model = Model(graph)

                # for each one of these,
                # create the node and add equivalent classes
                biogrid_id = 'BIOGRID:'+biogrid_num
                prefix = self.localtt[id_type]

                # TODO make these filters available as commandline options
                # geneidtypefilters='NCBIGene,OMIM,MGI,FlyBase,ZFIN,MGI,HGNC,
                #                   WormBase,XenBase,ENSEMBL,miRBase'.split(',')
                geneidtypefilters = 'NCBIGene,MGI,ENSEMBL,ZFIN,HGNC'.split(',')
                # proteinidtypefilters='HPRD,Swiss-Prot,NCBIProtein'
                if (speciesfilters is not None) \
                        and (organism_label.strip() in speciesfilters):
                    line_counter += 1
                    if (geneidtypefilters is not None) \
                            and (prefix in geneidtypefilters):
                        mapped_id = ':'.join((prefix, id_num))
                        model.addEquivalentClass(biogrid_id, mapped_id)
                    # this symbol will only get attached to the biogrid class
                    elif id_type == 'OFFICIAL_SYMBOL':
                        model.addClassToGraph(biogrid_id, id_num)
                    # elif (id_type == 'SYNONYM'):
                    #   FIXME - i am not sure these are synonyms, altids?
                    #   gu.addSynonym(g,biogrid_id,id_num)

                if not self.testMode and limit is not None and line_counter > limit:
                    break

        myzip.close()

        return

    def getTestSuite(self):
        import unittest
        from tests.test_biogrid import BioGridTestCase
        # TODO add InteractionAssoc tests
        # TODO add test about if all prefixes are mapped?

        test_suite = unittest.TestLoader().loadTestsFromTestCase(BioGridTestCase)

        return test_suite
