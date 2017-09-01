#!/usr/bin/env python3

import unittest
import logging
from tests.test_source import SourceTestCase
from dipper.sources.GWASCatalog import GWASCatalog
from dipper.graph.RDFGraph import RDFGraph
from rdflib import URIRef

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class GWASCatalogTestCase(SourceTestCase):

    def setUp(self):
        self.source = GWASCatalog('rdf_graph', True)
        self.source.test_ids = self._get_conf()['test_ids']
        self.source.settestonly(True)
        self._setDirToSource()
        return

    def tearDown(self):
        self.source = None
        return


class TestGwasSNPModel(unittest.TestCase):
    """
    Test the modelling of a  SNP to trait association
    from sample GWAS catalog data
    """

    def setUp(self):
        self.source = GWASCatalog('rdf_graph', True)
        self.test_data = {
            'snp_label': 'rs1491921-C',
            'chrom_num': '5',
            'chrom_pos': '21259029',
            'context': 'intergenic_variant',
            'allele_freq': '0.013',
            'trait': 'Diisocyanate-induced asthma',
            'trait_uri': 'http://www.ebi.ac.uk/efo/EFO_0006995, http://www.ebi.ac.uk/efo/EFO_0000270',
            'pvalue': '0.0000007',
            'merged': '0',
            'snp_id_current': '1491921',
            'mapped_gene': 'LOC102723561 - GUSBP1',
            'snp_gene_nums': '',
            'upstream_gene_num': '107986179',
            'downstream_gene_num': '107986180',
            'init_sample_desc': '74 European ancestry cases, 824 European ancestry controls',
            'replicated_sample_desc': 'NA',
            'platform': 'Illumina [1556551]',
            'pubmed': '25918132'
        }

    def tearDown(self):
        self.source = None
        self.efo_ontology = None

    def test_snp_type_resolution(self):
        """
        Given the label: rs1491921-C
        return dbSNP:rs1491921, snp
        """
        variant_curie, variant_type = \
            self.source._get_curie_and_type_from_id(self.test_data['snp_label'])

        self.assertEqual(variant_curie, "dbSNP:rs1491921")
        self.assertEqual(variant_type, 'snp')

    def test_snp_model(self):
        """
        Test output model of _add_snp_to_graph()
        """
        variant_curie, variant_type = \
            self.source._get_curie_and_type_from_id(self.test_data['snp_label'])

        self.source._add_snp_to_graph(
            variant_curie, self.test_data['snp_label'], self.test_data['chrom_num'],
            self.test_data['chrom_pos'], self.test_data['context'],
            self.test_data['allele_freq'])

        sparql_query = """
            SELECT ?snp
            WHERE {
                ?snp a OBO:SO_0000694,
                    OBO:SO_0001628 ;
                    rdfs:label "rs1491921-C" ;
                    faldo:location
                    <https://monarchinitiative.org/.well-known/genid/GRCh38chr5-21259029-21259029-Region> ;
                    OBO:RO_0002162 OBO:NCBITaxon_9606 ;
                    dc:description "0.013 [risk allele frequency]" .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr5-21259029-21259029-Region> a faldo:Region ;
                    faldo:begin <https://monarchinitiative.org/.well-known/genid/GRCh38chr5-21259029> ;
                    faldo:end <https://monarchinitiative.org/.well-known/genid/GRCh38chr5-21259029> .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr5-21259029> a faldo:Position ;
                    faldo:position 21259029 ;
                    faldo:reference OBO:CHR_GRCh38chr5 .
            }
        """
        # To debug
        #print(self.source.graph.serialize(format="turtle").decode("utf-8"))
        #self.assertTrue(False)

        sparql_output = self.source.graph.query(sparql_query)
        # Test that query passes and returns one row
        results = list(sparql_output)
        expected = [(URIRef(self.source.graph._getNode("dbSNP:rs1491921")),)]
        self.assertEqual(results, expected)

    def test_snp_gene_relation(self):
        """
        test the _add_snp_gene_relation function
        :return:
        """

        variant_curie, variant_type = \
            self.source._get_curie_and_type_from_id(self.test_data['snp_label'])

        self.source._add_snp_gene_relation(
            variant_curie, self.test_data['snp_gene_nums'],
            self.test_data['upstream_gene_num'],
            self.test_data['downstream_gene_num'])

        sparql_query = """
            SELECT ?snp
            WHERE {
                ?snp OBO:RO_0002528 <http://www.ncbi.nlm.nih.gov/gene/107986180> ;
                     OBO:RO_0002529 <http://www.ncbi.nlm.nih.gov/gene/107986179> .
            }
        """
        sparql_output = self.source.graph.query(sparql_query)
        # Test that query passes and returns one row
        results = list(sparql_output)
        expected = [(URIRef(self.source.graph._getNode("dbSNP:rs1491921")),)]
        self.assertEqual(results, expected)

    def test_deprecated_snp(self):
        """
        test the _add_deprecated_snp
        :return:
        """
        #fake data
        snp_id_current = '12345'
        merged = '1'

        variant_curie, variant_type = \
            self.source._get_curie_and_type_from_id(self.test_data['snp_label'])

        self.source._add_deprecated_snp(
            variant_curie, snp_id_current, merged,
            self.test_data['chrom_num'], self.test_data['chrom_pos'])

        sparql_query = """
            SELECT ?snp
            WHERE {
                ?snp a owl:NamedIndividual ;
                    OBO:IAO_0100001 dbSNP:rs12345 ;
                    owl:deprecated true .

                dbSNP:rs12345 MONARCH:cliqueLeader true .
            }
        """
        sparql_output = self.source.graph.query(sparql_query)
        # Test that query passes and returns one row
        results = list(sparql_output)
        expected = [(URIRef(self.source.graph._getNode("dbSNP:rs1491921")),)]
        self.assertEqual(results, expected)

    def test_snp_trait_association(self):
        """
        test the _add_variant_trait_association
        :return:
        """
        efo_ontology = RDFGraph()
        logger.info("Loading EFO ontology in separate rdf graph")
        efo_ontology.parse(self.source.files['efo']['url'], format='xml')
        efo_ontology.bind_all_namespaces()
        logger.info("Finished loading EFO ontology")

        variant_curie, variant_type = \
            self.source._get_curie_and_type_from_id(self.test_data['snp_label'])

        description = self.source._make_description(
            self.test_data['trait'], self.test_data['init_sample_desc'],
            self.test_data['replicated_sample_desc'],
            self.test_data['platform'], self.test_data['pvalue'])

        self.source._add_variant_trait_association(
            variant_curie, self.test_data['trait_uri'], efo_ontology,
            self.test_data['pubmed'], description)

        sparql_query = """
            SELECT ?snp
            WHERE {{
                <https://monarchinitiative.org/MONARCH_b630ef6046547f10> a OBAN:association ;
                    dc:description "{}" ;
                    OBO:RO_0002558 OBO:ECO_0000213 ;
                    dc:source PMID:25918132 ;
                    OBAN:association_has_object EFO:0000270 ;
                    OBAN:association_has_predicate OBO:RO_0002326 ;
                    OBAN:association_has_subject ?snp .

                <https://monarchinitiative.org/MONARCH_b70a05d8eb1c3d4b> a OBAN:association ;
                    OBO:RO_0002558 OBO:ECO_0000213 ;
                    dc:source PMID:25918132 ;
                    OBAN:association_has_object EFO:0006995 ;
                    OBAN:association_has_predicate OBO:RO_0002326 ;
                    OBAN:association_has_subject ?snp .

                EFO:0000270 a owl:Class ;
                    rdfs:label "asthma"^^xsd:string ;
                    rdfs:subClassOf OBO:DOID_4 .

                ?snp OBO:RO_0002326 EFO:0000270,
                        EFO:0006995 .

                PMID:25918132 a OBO:IAO_0000013 .
            }}
        """.format(description)
        sparql_output = self.source.graph.query(sparql_query)
        # Test that query passes and returns one row
        results = list(sparql_output)
        expected = [(URIRef(self.source.graph._getNode("dbSNP:rs1491921")),)]
        self.assertEqual(results, expected)


class TestGwasHaplotypeModel(unittest.TestCase):
    """
    Test the modelling of a  SNP to trait association
    from sample GWAS catalog data
    """

    def setUp(self):
        self.source = GWASCatalog('rdf_graph', True)
        self.test_data = {
            'snp_label': 'rs1329573-?; rs7020413-?; rs3824344-?; rs3758171-?',
            'chrom_num': '9;9;9;9',
            'chrom_pos': '36998996;37002118;37000690;36997420',
            'context': 'intron_variant; intron_variant; intron_variant; intron_variant',
            'allele_freq': 'NR',
            'trait': 'Intelligence',
            'trait_uri': 'http://www.ebi.ac.uk/efo/EFO_0004337',
            'pvalue': '0.00000004',
            'merged': '0',
            'snp_id_current': '',
            'mapped_gene': 'PAX5; PAX5; PAX5; PAX5',
            'snp_gene_nums': '',
            'upstream_gene_num': '107986179',
            'downstream_gene_num': '107986180',
            'init_sample_desc': '656 European ancestry individuals from ADHD families',
            'replicated_sample_desc': 'NA',
            'platform': 'Illumina [795637]',
            'pubmed': '22449649'
        }

    def tearDown(self):
        self.source = None

    def test_snp_model(self):
        """
        Test output model of _process_haplotype()
        self._process_haplotype(
                            variant_curie, strongest_snp_risk_allele,
                            chrom_num, chrom_pos, context,
                            risk_allele_frequency, mapped_gene, so_ontology)
        """
        variant_curie, variant_type = \
            self.source._get_curie_and_type_from_id(self.test_data['snp_label'])

        so_ontology = RDFGraph()
        logger.info("Loading SO ontology in separate rdf graph")
        so_ontology.parse(self.source.files['so']['url'], format='xml')
        so_ontology.bind_all_namespaces()
        logger.info("Finished loading SO ontology")

        self.source._process_haplotype(
            variant_curie, self.test_data['snp_label'], self.test_data['chrom_num'],
            self.test_data['chrom_pos'], self.test_data['context'],
            self.test_data['allele_freq'], self.test_data['mapped_gene'], so_ontology)

        sparql_query = """
            SELECT ?snp
            WHERE {
                :haplotype_bcb627b1f64039b0 a OBO:GENO_0000871 ;
                    rdfs:label "rs1329573-?; rs7020413-?; rs3824344-?; rs3758171-?" ;
                    OBO:GENO_0000382 ?snp,
                        dbSNP:rs3758171,
                        dbSNP:rs3824344,
                        dbSNP:rs7020413 ;
                OBO:GENO_0000418 <http://www.ncbi.nlm.nih.gov/gene/5079> ;
                OBO:RO_0002162 OBO:NCBITaxon_9606 .

                ?snp a OBO:SO_0000694,
                        OBO:SO_0001627 ;
                    rdfs:label "rs1329573-?" ;
                    faldo:location <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36998996-36998996-Region> ;
                    OBO:GENO_0000418 <http://www.ncbi.nlm.nih.gov/gene/5079> ;
                    OBO:RO_0002162 OBO:NCBITaxon_9606 .

                dbSNP:rs3758171 a OBO:SO_0000694,
                        OBO:SO_0001627 ;
                    rdfs:label "rs3758171-?" ;
                    faldo:location <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36997420-36997420-Region> ;
                    OBO:GENO_0000418 <http://www.ncbi.nlm.nih.gov/gene/5079> ;
                    OBO:RO_0002162 OBO:NCBITaxon_9606 .

                dbSNP:rs3824344 a OBO:SO_0000694,
                        OBO:SO_0001627 ;
                    rdfs:label "rs3824344-?" ;
                    faldo:location <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37000690-37000690-Region> ;
                    OBO:GENO_0000418 <http://www.ncbi.nlm.nih.gov/gene/5079> ;
                    OBO:RO_0002162 OBO:NCBITaxon_9606 .

                dbSNP:rs7020413 a OBO:SO_0000694,
                        OBO:SO_0001627 ;
                    rdfs:label "rs7020413-?" ;
                    faldo:location <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37002118-37002118-Region> ;
                    OBO:GENO_0000418 <http://www.ncbi.nlm.nih.gov/gene/5079> ;
                    OBO:RO_0002162 OBO:NCBITaxon_9606 .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36997420-36997420-Region> a faldo:Region ;
                    faldo:begin <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36997420> ;
                    faldo:end <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36997420> .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36998996-36998996-Region> a faldo:Region ;
                    faldo:begin <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36998996> ;
                    faldo:end <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36998996> .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37000690-37000690-Region> a faldo:Region ;
                    faldo:begin <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37000690> ;
                    faldo:end <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37000690> .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37002118-37002118-Region> a faldo:Region ;
                    faldo:begin <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37002118> ;
                    faldo:end <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37002118> .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36997420> a faldo:Position ;
                    faldo:position 36997420 ;
                    faldo:reference OBO:CHR_GRCh38chr9 .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-36998996> a faldo:Position ;
                    faldo:position 36998996 ;
                    faldo:reference OBO:CHR_GRCh38chr9 .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37000690> a faldo:Position ;
                    faldo:position 37000690 ;
                    faldo:reference OBO:CHR_GRCh38chr9 .

                <https://monarchinitiative.org/.well-known/genid/GRCh38chr9-37002118> a faldo:Position ;
                    faldo:position 37002118 ;
                    faldo:reference OBO:CHR_GRCh38chr9 .
            }
        """
        sparql_output = self.source.graph.query(sparql_query)
        # Test that query passes and returns one row
        results = list(sparql_output)
        expected = [(URIRef(self.source.graph._getNode("dbSNP:rs1329573")),)]
        self.assertEqual(results, expected)


if __name__ == '__main__':
    unittest.main()
