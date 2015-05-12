import csv
import re
import gzip
import os
import logging
import sys

from dipper import curie_map
from dipper import config
from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.Chem2DiseaseAssoc import Chem2DiseaseAssoc
from dipper.models.Genotype import Genotype
from dipper.models.G2PAssoc import G2PAssoc
from dipper.utils.GraphUtils import GraphUtils


logger = logging.getLogger(__name__)


class CTD(Source):
    """
    The Comparative Toxicogenomics Database (CTD) includes curated data describing cross-species chemical–gene/protein
    interactions and chemical– and gene–disease associations to illuminate molecular mechanisms underlying variable
    susceptibility and environmentally influenced diseases.

    Here, we fetch, parse, and convert data from CTD into triples, leveraging only the associations based on
    direct evidence (not using the inferred associations).  We currently only process the following associations:
        *chemical-disease
        *gene-pathway
        *gene-disease

    As part of a separate process through human-disease-ontology, we pull in the disease identifiers as mapped to
    MESH.

    """

    files = {
        'chemical_disease_interactions': {
            'file': 'CTD_chemicals_diseases.tsv.gz',
            'url': 'http://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz'
        },
        'gene_pathway': {
            'file': 'CTD_genes_pathways.tsv.gz',
            'url': 'http://ctdbase.org/reports/CTD_genes_pathways.tsv.gz'
        },
        'gene_disease': {
            'file': 'CTD_genes_diseases.tsv.gz',
            'url': 'http://ctdbase.org/reports/CTD_genes_diseases.tsv.gz'
        }
    }
    static_files = {
        'publications': {'file': 'CTD_curated_references.tsv.gz'}
    }

    def __init__(self):
        Source.__init__(self, 'ctd')
        self.dataset = Dataset('ctd', 'CTD', 'http://ctdbase.org', None, 'http://ctdbase.org/about/legal.jsp')

        if 'test_ids' not in config.get_config() or 'gene' not in config.get_config()['test_ids']:
            logger.warn("not configured with gene test ids.")
            self.test_geneids = []
        else:
            self.test_geneids = config.get_config()['test_ids']['gene']

        if 'test_ids' not in config.get_config() or 'disease' not in config.get_config()['test_ids']:
            logger.warn("not configured with disease test ids.")
            self.test_diseaseids = []
        else:
            self.test_diseaseids = config.get_config()['test_ids']['disease']

        return

    def fetch(self, is_dl_forced):
        """
        Override Source.fetch()
        Fetches resources from CTD using the CTD.files dictionary
        Args:
            :param is_dl_forced (bool): Force download
        Returns:
            :return None
        """
        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        """
        Override Source.parse()
        Parses version and interaction information from CTD
        Args:
            :param limit (int, optional) limit the number of rows processed
        Returns:
            :return None
        """
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        logger.info("Parsing files...")
        pub_map = dict()
        file_path = '/'.join((self.rawdir,
                              self.static_files['publications']['file']))
        if os.path.exists(file_path) is True:
            pub_map = self._parse_publication_file(
                self.static_files['publications']['file']
            )

        if self.testOnly:
            self.testMode = True

        self._parse_ctd_file(limit, self.files['chemical_disease_interactions']['file'], pub_map)
        self._parse_ctd_file(limit, self.files['gene_pathway']['file'])
        self._parse_ctd_file(limit, self.files['gene_disease']['file'])
        self.load_bindings()
        logger.info("Done parsing files.")

        return

    def _parse_ctd_file(self, limit, file, pub_map=None):
        """
        Parses files in CTD.files dictionary
        Args:
            :param limit (int): limit the number of rows processed
            :param file (str): file name (must be defined in CTD.file)
            :param pub_map(dict, optional):
                    publication mapping dictionary
                    generated with _split_pub_ids_by_evidence()
        Returns:
            :return None
        """
        row_count = 0
        version_pattern = re.compile('^# Report created: (.+)$')
        is_versioned = False
        file_path = '/'.join((self.rawdir, file))
        with gzip.open(file_path, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                # Scan the header lines until we get the version
                # There is no official version sp we are using
                # the upload timestamp instead
                if is_versioned is False:
                    match = re.match(version_pattern, ' '.join(row))
                    if match:
                        version = re.sub(r'\s|:', '-', match.group(1))
                        self.dataset.setVersion(version)
                        is_versioned = True
                elif re.match('^#', ' '.join(row)):
                    next
                else:
                    if file == self.files['chemical_disease_interactions']['file']:
                        self._process_interactions(row, pub_map)
                    elif file == self.files['gene_pathway']['file']:
                        self._process_pathway(row)
                    elif file == self.files['gene_disease']['file']:
                        #self._process_disease2gene(row)
                        continue
                    row_count += 1
                    if (not self.testMode) and (limit is not None and row_count >= limit):
                        break
        return

    def _process_pathway(self, row):
        """
        Process row of CTD data from CTD_genes_pathways.tsv.gz
        and generate triples
        Args:
            :param row (list): row of CTD data
        Returns:
            :return None
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        self._check_list_len(row, 4)
        gu = GraphUtils(curie_map.get())
        (gene_symbol, gene_id, pathway_name, pathway_id) = row

        if self.testMode and (int(gene_id) not in self.test_geneids):
            return

        entrez_id = 'NCBIGene:'+gene_id
        # Adding all pathways as classes, may refactor in the future
        # to only add Reactome pathways from CTD
        gu.addClassToGraph(g, pathway_id, pathway_name)
        gu.addSubclass(g, self._get_class_id('signal transduction'), pathway_id)
        gu.addClassToGraph(g, entrez_id, gene_symbol)
        gu.addInvolvedIn(g, entrez_id, pathway_id)

        # if re.match(re.compile('^REACT'),pathway_id):
        #    gu.addClassToGraph(g, pathway_id, pathway_name)

        return

    def _process_interactions(self, row, pub_map):
        """
        Process row of CTD data from CTD_genes_pathways.tsv.gz
        and generate triples.
        Only create associations based on direct evidence (not using the inferred-via-gene)
        Args:
            :param row (list): row of CTD data
            :param pub_map(dict, optional): publication mapping dictionary
        Returns:
            :return None
        """
        self._check_list_len(row, 10)
        (chem_name, chem_id, cas_rn, disease_name, disease_id, direct_evidence,
         inferred_gene_symbol, inference_score, omim_ids, pubmed_ids) = row
        evidence_pattern = re.compile('^therapeutic|marker\/mechanism$')
        dual_evidence = re.compile('^marker\/mechanism\|therapeutic$')

        # filter on those diseases that are mapped to omim ids in the test set
        intersect = list(set(['OMIM:'+str(i) for i in omim_ids.split('|')]+[disease_id]) & set(self.test_diseaseids))
        if self.testMode and len(intersect) < 1:
            return

        if direct_evidence == '':
            next
        elif re.match(evidence_pattern, direct_evidence):
            reference_list = self._process_pubmed_ids(pubmed_ids)
            self._make_association(chem_id, disease_id, direct_evidence, reference_list)
        elif ((re.match(dual_evidence, direct_evidence))
              and (len(pub_map.keys()) > 0)):
            reference_list = self._process_pubmed_ids(pubmed_ids)
            pub_evidence_map = self._split_pub_ids_by_evidence(reference_list, chem_id, disease_id, pub_map)
            for val in direct_evidence.split('|'):
                self._make_association(chem_id, disease_id, val, pub_evidence_map[val])
        else:
            # there's dual evidence, but haven't mapped the pubs
            logger.info("Dual evidence for %s and %s but no pub map", chem_name, disease_id)

        return

    def _process_disease2gene(self, row):
        """
        Process row of CTD data from CTD_genes_pathways.tsv.gz
        and generate triples
        Args:
            :param row (list): row of CTD data
        Returns:
            :return None
        """
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        self._check_list_len(row, 9)
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())
        (gene_symbol, gene_id, disease_name, disease_id, direct_evidence,
         inference_chemical_name, inference_score, omim_ids, pubmed_ids) = row

        if self.testMode and (int(gene_id) not in self.test_geneids):
            return

        #
        if direct_evidence != 'marker/mechanism':
            next

        rel = gu.object_properties['correlates_with']
        gene_id = 'NCBIGene:'+gene_id

        # If there is only one OMIM ID for the Disease ID or in the omim_ids list,
        # use the OMIM ID preferentially over any MeSH ID.
        if re.match('OMIM:.*', disease_id) and re.match('.*\|.*', omim_ids) and omim_ids != '' and omim_ids is not None:
            # Currently no entries with an OMIM ID as the disease ID and additional OMIM IDs in the omim_ids field.
            # This is when the disease ID is an OMIM ID and there is more than one OMIM entry in omim_ids.
            preferred_disease_id = disease_id
            #print('Disease ID is OMIM ID: '+disease_id+' & More than one OMIM ID in omim_ids: '+omim_ids)
        elif re.match('OMIM:.*', disease_id) and not re.match('.*\|.*', omim_ids) and omim_ids != '' and omim_ids is not None:
            # This is when the disease ID is an OMIM ID and there is only one OMIM entry in omim_ids.
            preferred_disease_id = disease_id
            if disease_id != ('OMIM:'+omim_ids):
                #print('IDs do not match! '+disease_id+' != '+omim_ids)
                alternate_disease_ids = ['OMIM:'+omim_ids]
            #print('Disease ID is OMIM ID: '+disease_id+' & Single OMIM ID in omim_ids: '+omim_ids)
        elif not re.match('OMIM:.*', disease_id) and omim_ids != '' and omim_ids is not None and not re.match('.*\|.*', omim_ids):
            # This is when the disease ID is not an OMIM ID and there is only one OMIM entry in omim_ids.
            #print('Disease ID is not OMIM ID: '+disease_id+' & Single OMIM ID in omim_ids: '+omim_ids)
            preferred_disease_id = 'OMIM:'+omim_ids
        elif not re.match('OMIM:.*', disease_id) and omim_ids != '' and omim_ids is not None and re.match('.*\|.*', omim_ids):
            # This is when the disease ID is an OMIM ID and there is more than one OMIM entry in omim_ids.
            #print('Disease ID is not OMIM ID: '+disease_id+' & More than one OMIM ID in omim_ids: '+omim_ids)
            #FIXME:Correct handling, or avoid given the existence of multiple OMIM IDs?
            preferred_disease_id = disease_id
        else:
            # This is when the omim_ids field is empty. May be either a MeSH ID or an OMIM ID.
            preferred_disease_id = disease_id



        # Make an association ID.
        assoc_id = self.make_id((preferred_disease_id+gene_id))

        # we actually want the association between the gene and the disease to be via an alternate locus
        # not the "wildtype" gene itself.
        # so we make an anonymous alternate locus, and put that in the association.
        alt_locus = '_'+gene_id+'-'+preferred_disease_id+'VL'
        alt_label = 'some variant of '+gene_symbol+' that causes '+disease_name
        gu.addIndividualToGraph(g, alt_locus, alt_label, geno.genoparts['variant_locus'])
        geno.addAlleleOfGene(alt_locus, gene_id)
        # Add the disease to gene relationship.
        assoc = G2PAssoc(assoc_id, alt_locus, disease_id, None, None)
        assoc.loadAllProperties(g)
        assoc.addAssociationToGraph(g)

        return

    def _make_association(self, chem_id, disease_id, direct_evidence, pubmed_ids):
        """
        Make a reified chemical-phenotype association
        using the Chem2DiseaseAssoc class
        Args:
            :param chem_id
            :param disease_id
            :param direct_evidence
            :param pubmed_ids
        Returns:
            :return None
        """

        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        assoc_id = self.make_ctd_chem_disease_assoc_id(chem_id,disease_id,direct_evidence)
        evidence_code = self._get_evidence_code('TAS')
        chem_mesh_id = 'MESH:'+chem_id
        relationship = self._get_relationship_id(direct_evidence)
        assoc = Chem2DiseaseAssoc(assoc_id, chem_mesh_id, disease_id,
                                  pubmed_ids, relationship, evidence_code)
        assoc.loadAllProperties(g)
        assoc.addAssociationNodeToGraph(g)

        return g

    def make_ctd_chem_disease_assoc_id(self, chem_id, disease_id, direct_evidence):
        assoc_id = self.make_id('ctd' + chem_id + disease_id + direct_evidence)

        return assoc_id

    def _process_pubmed_ids(self, pubmed_ids):
        """
        Take a list of pubmed IDs and add PMID prefix
        Args:
            :param pubmed_ids -  string representing publication
                                 ids seperated by a | symbol
        Returns:
            :return list: Pubmed curies
        """
        id_list = pubmed_ids.split('|')
        for (i, val) in enumerate(id_list):
            id_list[i] = 'PMID:'+val
        return id_list

    def _get_evidence_code(self, evidence):
        """
        Get curie for evidence class label
        Args:
            :param evidence (str): evidence label
        Label:
            :return str: curie for evidence label from ECO
        """
        ECO_MAP = {
            'TAS': 'ECO:0000033'
        }
        return ECO_MAP[evidence]

    def _get_relationship_id(self, rel):
        """
        Get curie from relationship property label
        Args:
            :param rel (str): relationship label
        Returns:
            :return str: curie for relationship label
        """
        REL_MAP = {
            'therapeutic': 'MONARCH:treats',
            'marker/mechanism': 'MONARCH:causes'
        }
        return REL_MAP[rel]

    def _get_class_id(self, cls):
        """
        Fet curie from CLASS_MAP dictionary
        Args:
            :param cls (str): class label
        Returns:
            :return str: curie for class label
        """
        CLASS_MAP = {
            'pathway': 'PW:0000001',
            'signal transduction': 'GO:0007165'
        }
        return CLASS_MAP[cls]

    def _split_pub_ids_by_evidence(self, pub_ids, chem_id, disease_id, pub_map):
        """
        Split a list of ambiguous sources to chemical-phenotype relationship
        Args:
            :param pub_ids (list): list of publication IDs
            :param disease_id (str): disease curie
            :param chem_id (str): chemical curie
            :param pub_map (dict): publication dictionary
        Returns:
            :return dictionary in the following structure:
                     {
                      'therapeutic': [1234,2345],
                      'marker/mechanism': [4567,5678]
                     }
        """
        publication = {'therapeutic': [], 'marker/mechanism': []}
        for val in pub_ids:
            key = disease_id+'-'+chem_id+'-'+val
            if pub_map.get(key):
                publication[pub_map[key]].append(val)
            else:
                logger.error("Could not disambiguate publication "
                             "for disease id: %s"
                             "\nchemical id: %s"
                             "\npublication id: %s", disease_id, chem_id, str(val))
                #sys.exit(1)

        return publication

    def _parse_publication_file(self, file):
        """
        Parse publication file found in CTD.static_files
        Args:
            :param file (str): file name
        Returns:
            :return dict: key containing the chemID, phenotypeID, and pubID
                           mapped to relationship
        """
        row_count = 0
        pub_map = dict()
        file_path = '/'.join((self.rawdir, file))
        with gzip.open(file_path, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                self._check_list_len(row, 10)
                # catch comment lines
                if re.match('^#', ' '.join(row)):
                    next
                else:
                    (pub_id, disease_label, disease_id, disease_cat, evidence,
                     chem_label, chem_id, cas_rn, gene_symbol, gene_acc) = row
                    key = disease_id+'-'+chem_id+'-PMID:'+pub_id
                    if chem_id is '' or disease_id is '':
                        next
                    elif pub_map.get(key) is not None:
                        logger.error("Ambiguous publication mapping for"
                                     " key: %s\n "
                                     "This is not yet handled by Dipper", key)
                        sys.exit(1)
                    else:
                        pub_map[key] = evidence

        return pub_map

    def getTestSuite(self):
        import unittest
        from tests.test_ctd import CTDTestCase
        from tests.test_interactions import InteractionsTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(CTDTestCase)
        test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(InteractionsTestCase))

        return test_suite
