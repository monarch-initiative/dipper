import os
import logging
import re

from datetime import datetime
from stat import ST_CTIME
from zipfile import ZipFile

from dipper import config
from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.assoc.InteractionAssoc import InteractionAssoc
from dipper.models.Dataset import Dataset

__author__ = 'nicole'

logger = logging.getLogger(__name__)
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
        super().__init__(graph_type, are_bnodes_skolemized, 'biogrid')

        self.tax_ids = tax_ids

        self.dataset = Dataset(
            'biogrid', 'The BioGrid', 'http://thebiogrid.org/', None,
            'http://wiki.thebiogrid.org/doku.php/terms_and_conditions')

        # Defaults
        # our favorite animals
        # taxids = [9606,10090,10116,7227,7955,6239,8355]
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090, 7955]

        if 'test_ids' not in config.get_config() or \
                'gene' not in config.get_config()['test_ids']:
            logger.warning("not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['gene']

        # data-source specific warnings
        # (will be removed when issues are cleared)
        logger.warning(
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
            version = \
                re.match(r'BIOGRID-ALL-(\d+\.\d+\.\d+)\.mitab.txt', fname)
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

        logger.info("Loaded %d test graph nodes", len(self.testgraph))
        logger.info("Loaded %d full graph nodes", len(self.graph))

        return

    def _get_interactions(self, limit):
        logger.info("getting interactions")
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
                    logger.debug("Skipping header line")
                    continue
                line_counter += 1
                line = line.decode().strip()
                # print(line)
                (interactor_a, interactor_b, alt_ids_a, alt_ids_b, aliases_a,
                 aliases_b, detection_method, pub_author, pub_id, taxid_a,
                 taxid_b, interaction_type, source_db, interaction_id,
                 confidence_val) = line.split('\t')

                # get the actual gene ids,
                # typically formated like: gene/locuslink:351|BIOGRID:106848
                gene_a_num = re.search(
                    r'locuslink\:(\d+)\|?', interactor_a).groups()[0]
                gene_b_num = re.search(
                    r'locuslink\:(\d+)\|?', interactor_b).groups()[0]

                if self.testMode:
                    g = self.testgraph
                    # skip any genes that don't match our test set
                    if (int(gene_a_num) not in self.test_ids) or\
                            (int(gene_b_num) not in self.test_ids):
                        continue
                else:
                    g = self.graph
                    # when not in test mode, filter by taxon
                    if int(re.sub(r'taxid:', '', taxid_a.rstrip())) not in\
                            self.tax_ids or\
                            int(re.sub(
                                r'taxid:', '', taxid_b.rstrip())) not in\
                            self.tax_ids:
                        continue
                    else:
                        matchcounter += 1

                gene_a = 'NCBIGene:'+gene_a_num
                gene_b = 'NCBIGene:'+gene_b_num

                # get the interaction type
                # psi-mi:"MI:0407"(direct interaction)
                int_type = re.search(r'MI:\d+', interaction_type).group()
                rel = self._map_MI_to_RO(int_type)

                # scrub pubmed-->PMID prefix
                pub_id = re.sub(r'pubmed', 'PMID', pub_id)
                # remove bogus whitespace
                pub_id = pub_id.strip()

                # get the method, and convert to evidence code
                det_code = re.search(r'MI:\d+', detection_method).group()
                evidence = self._map_MI_to_ECO(det_code)

                # note that the interaction_id is some kind of internal biogrid
                # identifier that does not map to a public URI.
                # we will construct a monarch identifier from this

                assoc = InteractionAssoc(g, self.name, gene_a, gene_b, rel)
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

        logger.info("getting identifier mapping")
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
                    g = self.testgraph
                    # skip any genes that don't match our test set
                    if int(biogrid_num) not in self.biogrid_ids:
                        continue
                else:
                    g = self.graph

                model = Model(g)

                # for each one of these,
                # create the node and add equivalent classes
                biogrid_id = 'BIOGRID:'+biogrid_num
                prefix = self._map_idtype_to_prefix(id_type)

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

                if not self.testMode and limit is not None \
                        and line_counter > limit:
                    break

        myzip.close()

        return

    @staticmethod
    def _map_MI_to_RO(mi_id):
        rel = InteractionAssoc.interaction_object_properties
        mi_ro_map = {
            # colocalization
            'MI:0403': rel['colocalizes_with'],
            # direct interaction
            'MI:0407': rel['interacts_with'],
            # synthetic genetic interaction defined by inequality
            'MI:0794': rel['genetically_interacts_with'],
            # suppressive genetic interaction defined by inequality
            'MI:0796': rel['genetically_interacts_with'],
            # additive genetic interaction defined by inequality
            'MI:0799': rel['genetically_interacts_with'],
            # association
            'MI:0914': rel['interacts_with'],
            # physical association
            'MI:0915': rel['interacts_with']
        }

        ro_id = rel['interacts_with']  # default
        if mi_id in mi_ro_map:
            ro_id = mi_ro_map.get(mi_id)

        return ro_id

    @staticmethod
    def _map_MI_to_ECO(mi_id):
        eco_id = 'ECO:0000006'  # default to experimental evidence
        mi_to_eco_map = {
            'MI:0018': 'ECO:0000068',  # yeast two-hybrid
            'MI:0004': 'ECO:0000079',  # affinity chromatography
            'MI:0047': 'ECO:0000076',  # far western blotting
            'MI:0055': 'ECO:0000021',  # should be FRET, but using physical_interaction FIXME
            'MI:0090': 'ECO:0000012',  # desired: protein complementation, using: functional complementation
            'MI:0096': 'ECO:0000085',  # desired: pull down, using: immunoprecipitation
            'MI:0114': 'ECO:0000324',  # desired: x-ray crystallography, using: imaging assay
            'MI:0254': 'ECO:0000011',  # desired: genetic interference, using: genetic interaction evidence
            'MI:0401': 'ECO:0000172',  # desired: biochemical, using: biochemical trait evidence
            'MI:0415': 'ECO:0000005',  # desired: enzymatic study, using: enzyme assay evidence
            'MI:0428': 'ECO:0000324',  # imaging
            'MI:0686': 'ECO:0000006',  # desired: unspecified, using: experimental evidence
            'MI:1313': 'ECO:0000006'   # None?
        }
        if mi_id in mi_to_eco_map:
            eco_id = mi_to_eco_map.get(mi_id)
        else:
            logger.warning(
                "unmapped code %s. Defaulting to experimental_evidence", mi_id)

        return eco_id

    @staticmethod
    def _map_idtype_to_prefix(idtype):
        """
        Here we need to reformat the BioGrid source prefixes
        to standard ones used in our curie-map.
        :param idtype:
        :return:

        """
        prefix = idtype
        idtype_to_prefix_map = {
            'XENBASE': 'XenBase',
            'TREMBL': 'TrEMBL',
            'MGI': 'MGI',
            'REFSEQ_DNA_ACCESSION': 'RefSeqNA',
            'MAIZEGDB': 'MaizeGDB',
            'BEEBASE': 'BeeBase',
            'ENSEMBL': 'ENSEMBL',
            'TAIR': 'TAIR',
            'GENBANK_DNA_GI': 'NCBIgi',
            'CGNC': 'CGNC',
            'RGD': 'RGD',
            'GENBANK_GENOMIC_DNA_GI': 'NCBIgi',
            'SWISSPROT': 'Swiss-Prot',
            'MIM': 'OMIM',
            'FLYBASE': 'FlyBase',
            'VEGA': 'VEGA',
            'ANIMALQTLDB': 'AQTLDB',
            'ENTREZ_GENE_ETG': 'ETG',
            'HPRD': 'HPRD',
            'APHIDBASE': 'APHIDBASE',
            'GENBANK_PROTEIN_ACCESSION': 'NCBIProtein',
            'ENTREZ_GENE': 'NCBIGene',
            'SGD': 'SGD',
            'GENBANK_GENOMIC_DNA_ACCESSION': 'NCBIGenome',
            'BGD': 'BGD',
            'WORMBASE': 'WormBase',
            'ZFIN': 'ZFIN',
            'DICTYBASE': 'dictyBase',
            'ECOGENE': 'ECOGENE',
            'BIOGRID': 'BIOGRID',
            'GENBANK_DNA_ACCESSION': 'NCBILocus',
            'VECTORBASE': 'VectorBase',
            'MIRBASE': 'miRBase',
            'IMGT/GENE-DB': 'IGMT',
            'HGNC': 'HGNC',
            'SYSTEMATIC_NAME': None,
            'OFFICIAL_SYMBOL': None,
            'REFSEQ_GENOMIC_DNA_ACCESSION': 'NCBILocus',
            'GENBANK_PROTEIN_GI': 'NCBIgi',
            'REFSEQ_PROTEIN_ACCESSION': 'RefSeqProt',
            'SYNONYM': None,
            'GRID_LEGACY': None,
            # the following showed up in 3.3.124
            'UNIPROT-ACCESSION': 'UniprotKB',
            'SWISS-PROT': 'Swiss-Prot',
            'OFFICIAL SYMBOL': None,
            'ENSEMBL RNA': None,
            'GRID LEGACY': None,
            'ENSEMBL PROTEIN': None,
            'REFSEQ-RNA-GI': None,
            'REFSEQ-RNA-ACCESSION': None,
            'REFSEQ-PROTEIN-GI': None,
            'REFSEQ-PROTEIN-ACCESSION-VERSIONED': None,
            'REFSEQ-PROTEIN-ACCESSION': None,
            'REFSEQ-LEGACY': None,
            'SYSTEMATIC NAME': None,
            'ORDERED LOCUS': None,
            'UNIPROT-ISOFORM': 'UniprotKB',
            'ENSEMBL GENE': 'ENSEMBL',
            'CGD': None,  # Not sure what this is?
            'WORMBASE-OLD': 'WormBase'

        }
        if idtype in idtype_to_prefix_map:
            prefix = idtype_to_prefix_map.get(idtype)
        else:
            logger.warning("unmapped prefix %s", prefix)

        return prefix

    def getTestSuite(self):
        import unittest
        from tests.test_biogrid import BioGridTestCase
        # TODO add InteractionAssoc tests
        # TODO add test about if all prefixes are mapped?

        test_suite = \
            unittest.TestLoader().loadTestsFromTestCase(BioGridTestCase)

        return test_suite
