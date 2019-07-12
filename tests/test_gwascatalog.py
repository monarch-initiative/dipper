#!/usr/bin/env python3

import unittest
import logging
from dipper.sources.GWASCatalog import GWASCatalog
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.TestUtils import TestUtils

logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)


class TestGwasSNPModel(unittest.TestCase):
    """
    Test the modelling of a  SNP to trait association
    from sample GWAS catalog data
    """

    def setUp(self):
        self.test_util = TestUtils()
        self.source = GWASCatalog('rdf_graph', True)
        self.source.graph = RDFGraph(True)  # Reset graph
        self.source.graph.bind_all_namespaces()
        self.test_data = {
            'snp_label': 'rs1491921-C',
            'chrom_num': '5',
            'chrom_pos': '21259029',
            'context': 'intergenic_variant',
            'allele_freq': '0.013',
            'trait': 'response to diisocyanate, asthma',
            'trait_uri': 'http://www.ebi.ac.uk/efo/EFO_0006995, http://www.ebi.ac.uk/efo/EFO_0003949',
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
        self.assertTrue(len(list(self.source.graph)) == 0)
        variant_curie, variant_type = self.source._get_curie_and_type_from_id(
            self.test_data['snp_label'])

        self.assertEqual(variant_curie, "dbSNP:rs1491921")
        self.assertEqual(variant_type, 'snp')

    def test_snp_model(self):
        """
        Test output model of _add_snp_to_graph()
        """
        self.assertTrue(len(list(self.source.graph)) == 0)
        variant_curie, variant_type = self.source._get_curie_and_type_from_id(
            self.test_data['snp_label'])

        self.source._add_snp_to_graph(
            variant_curie, self.test_data['snp_label'], self.test_data['chrom_num'],
            self.test_data['chrom_pos'], self.test_data['context'],
            self.test_data['allele_freq'])

        triples = """
    dbSNP:rs1491921 a OBO:SO_0000694, OBO:SO_0001628 ;
        rdfs:label "rs1491921-C" ;
        faldo:location  <https://monarchinitiative.org/.well-known/genid/bee62aebc378782c23ad> ;
        OBO:RO_0002162 OBO:NCBITaxon_9606 ;
        dc:description "0.013 [risk allele frequency]" .

    <https://monarchinitiative.org/.well-known/genid/bee62aebc378782c23ad> a faldo:Region ;
        rdfs:label 'GRCh38chr5-21259029-21259029-Region' ;
        faldo:begin <https://monarchinitiative.org/.well-known/genid/bed3341497b6bee94927> ;
        faldo:end <https://monarchinitiative.org/.well-known/genid/bed3341497b6bee94927> .

    <https://monarchinitiative.org/.well-known/genid/bed3341497b6bee94927> a faldo:Position ;
        rdfs:label 'GRCh38chr5-21259029' ;
        faldo:position 21259029 ;
        faldo:reference OBO:CHR_GRCh38chr5 .

    OBO:NCBITaxon_9606 biolink:category biolink:OrganismTaxon . 
    dbSNP:rs1491921 biolink:category biolink:SequenceVariant .
    <https://monarchinitiative.org/.well-known/genid/bed3341497b6bee94927> biolink:category biolink:GenomicSequenceLocalization .
    <https://monarchinitiative.org/.well-known/genid/bee62aebc378782c23ad> biolink:category biolink:GenomicSequenceLocalization .
"""
        # dbg
        LOG.debug(
            "Reference graph: %s",
            self.source.graph.serialize(format="turtle").decode("utf-8"))

        self.assertTrue(self.test_util.test_graph_equality(triples, self.source.graph))

    def test_snp_gene_relation(self):
        """
        test the _add_snp_gene_relation function
        :return:
        """
        self.assertTrue(len(list(self.source.graph)) == 0)
        variant_curie, variant_type = self.source._get_curie_and_type_from_id(
            self.test_data['snp_label'])

        self.source._add_snp_gene_relation(
            variant_curie, self.test_data['snp_gene_nums'],
            self.test_data['upstream_gene_num'],
            self.test_data['downstream_gene_num'])

        triples = """
        dbSNP:rs1491921 OBO:RO_0002528 ENSEMBL:107986180 ;
            OBO:RO_0002529 ENSEMBL:107986179 .
            
        ENSEMBL:107986180 biolink:category biolink:Gene .
        ENSEMBL:107986179 biolink:category biolink:Gene .
        dbSNP:rs1491921 biolink:category biolink:SequenceVariant .
        """
        self.assertTrue(self.test_util.test_graph_equality(triples, self.source.graph))

    def test_deprecated_snp(self):
        """
        test the _add_deprecated_snp
        :return:
        """
        self.assertTrue(len(list(self.source.graph)) == 0)
        # fake data
        snp_id_current = '12345'
        merged = '1'

        variant_curie, variant_type = self.source._get_curie_and_type_from_id(
            self.test_data['snp_label'])

        self.source._add_deprecated_snp(
            variant_curie, snp_id_current, merged,
            self.test_data['chrom_num'], self.test_data['chrom_pos'])

        triples = """
        dbSNP:rs1491921 a owl:NamedIndividual ;
            OBO:IAO_0100001 dbSNP:rs12345 ;
            owl:deprecated true .

        dbSNP:rs12345 MONARCH:cliqueLeader true .
        
        dbSNP:rs12345 biolink:category biolink:SequenceVariant .
        dbSNP:rs1491921 biolink:category biolink:SequenceVariant .
        """
        self.assertTrue(self.test_util.test_graph_equality(triples, self.source.graph))

    def test_snp_trait_association(self):
        """
        test the _add_variant_trait_association
        :return:
        """
        self.assertTrue(len(list(self.source.graph)) == 0)
        mondo_data = {"graphs" : [ {"nodes" : [ {"meta": {"xrefs": [{"val": "EFO:0003949"}]}}]}]}

        variant_curie, variant_type = self.source._get_curie_and_type_from_id(
            self.test_data['snp_label'])

        description = self.source._make_description(
            self.test_data['trait'], self.test_data['init_sample_desc'],
            self.test_data['replicated_sample_desc'],
            self.test_data['platform'], self.test_data['pvalue'])

        self.source._add_variant_trait_association(
            variant_curie, self.test_data['trait_uri'],
            self.test_data['trait'], mondo_data,
            self.test_data['pubmed'], description)

        triples = """

    MONARCH:bffc7a930c08cc8fe931 a OBAN:association ;
        dc:description "{0}" ;
        OBO:RO_0002558 OBO:ECO_0000213 ;
        dc:source PMID:25918132 ;
        OBAN:association_has_object EFO:0003949 ;
        OBAN:association_has_predicate RO:0003304 ;
        OBAN:association_has_subject dbSNP:rs1491921 .

    MONARCH:bff9b97458d67ed7f517 a OBAN:association ;
        dc:description "{0}" ;
        OBO:RO_0002558 OBO:ECO_0000213 ;
        dc:source PMID:25918132 ;
        OBAN:association_has_object EFO:0006995 ;
        OBAN:association_has_predicate RO:0003304 ;
        OBAN:association_has_subject dbSNP:rs1491921 .

    EFO:0006995 a owl:Class ;
        rdfs:label "response to diisocyanate" ;
        rdfs:subClassOf UPHENO:0001001 .

    dbSNP:rs1491921 RO:0003304 EFO:0003949,
            EFO:0006995 .

    PMID:25918132 a OBO:IAO_0000013 ;
        biolink:category biolink:Publication .

    ECO:0000213 biolink:category biolink:EvidenceType .
    EFO:0003949 biolink:category biolink:PhenotypicFeature .
    EFO:0006995 biolink:category biolink:PhenotypicFeature .
    dbSNP:rs1491921 biolink:category biolink:Genotype .
    PMID:25918132 biolink:category biolink:EvidenceType .
        
        """.format(description)

        # dbg
        LOG.debug(
           "Reference graph: %s",
           self.source.graph.serialize(format="turtle").decode("utf-8"))
        self.assertTrue(self.test_util.test_graph_equality(triples, self.source.graph))


class TestGwasHaplotypeModel(unittest.TestCase):
    """
    Test the modelling of a  SNP to trait association
    from sample GWAS catalog data
    """

    def setUp(self):
        self.test_util = TestUtils()
        self.source = GWASCatalog('rdf_graph', True)
        self.source.graph = RDFGraph(True)
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
        self.assertTrue(len(list(self.source.graph)) == 0)
        variant_curie, variant_type = self.source._get_curie_and_type_from_id(
            self.test_data['snp_label'])

        so_ontology = RDFGraph()
        LOG.info("Loading SO ontology in separate rdf graph")
        so_ontology.parse(self.source.files['so']['url'], format='xml')
        so_ontology.bind_all_namespaces()
        LOG.info("Finished loading SO ontology")

        self.source._process_haplotype(
            variant_curie, self.test_data['snp_label'], self.test_data['chrom_num'],
            self.test_data['chrom_pos'], self.test_data['context'],
            self.test_data['allele_freq'], self.test_data['mapped_gene'], so_ontology)

        triples = """
NCBITaxon:9606 biolink:category biolink:OrganismTaxon .

:haplotype_bb627b1f64039b0f751a a SO:0001024 ;
    rdfs:label "rs1329573-?; rs7020413-?; rs3824344-?; rs3758171-?" ;
    GENO:0000382 dbSNP:rs1329573,
        dbSNP:rs3758171,
        dbSNP:rs3824344,
        dbSNP:rs7020413 ;
    GENO:0000418 HGNC:8619 ;
    RO:0002162 NCBITaxon:9606 .

dbSNP:rs1329573 a SO:0000694,
        SO:0001627 ;
    rdfs:label "rs1329573-?" ;
    faldo:location <https://monarchinitiative.org/.well-known/genid/b3fad5df82cdfb283329> ;
    GENO:0000418 HGNC:8619 ;
    RO:0002162 NCBITaxon:9606 .

dbSNP:rs3758171 a SO:0000694,
        SO:0001627 ;
    rdfs:label "rs3758171-?" ;
    faldo:location <https://monarchinitiative.org/.well-known/genid/b25a2da36647bdd71be3> ;
    GENO:0000418 HGNC:8619 ;
    RO:0002162 NCBITaxon:9606 .

dbSNP:rs3824344 a SO:0000694,
        OBO:SO_0001627 ;
    rdfs:label "rs3824344-?" ;
    faldo:location <https://monarchinitiative.org/.well-known/genid/b096a3e94e32fe23374a> ;
    GENO:0000418 HGNC:8619 ;
    RO:0002162 NCBITaxon:9606 .

dbSNP:rs7020413 a SO:0000694,
        SO:0001627 ;
    rdfs:label "rs7020413-?" ;
    faldo:location <https://monarchinitiative.org/.well-known/genid/bbb252d9b6cd02e9880a> ;
    GENO:0000418 HGNC:8619 ;
    RO:0002162 NCBITaxon:9606 .

<https://monarchinitiative.org/.well-known/genid/b25a2da36647bdd71be3> a faldo:Region ;
    rdfs:label "GRCh38chr9-36997420-36997420-Region";
    faldo:begin <https://monarchinitiative.org/.well-known/genid/b21985847fe0774084eb> ;
    faldo:end <https://monarchinitiative.org/.well-known/genid/b21985847fe0774084eb> .

<https://monarchinitiative.org/.well-known/genid/b3fad5df82cdfb283329> a faldo:Region ;
    rdfs:label "GRCh38chr9-36998996-36998996-Region";
    faldo:begin <https://monarchinitiative.org/.well-known/genid/b55051762f8d5a3dbeb5> ;
    faldo:end <https://monarchinitiative.org/.well-known/genid/b55051762f8d5a3dbeb5> .

<https://monarchinitiative.org/.well-known/genid/b096a3e94e32fe23374a> a faldo:Region ;
    rdfs:label "GRCh38chr9-37000690-37000690-Region";
    faldo:begin <https://monarchinitiative.org/.well-known/genid/b5d61dbc7958a979d046> ;
    faldo:end <https://monarchinitiative.org/.well-known/genid/b5d61dbc7958a979d046> .

<https://monarchinitiative.org/.well-known/genid/bbb252d9b6cd02e9880a> a faldo:Region ;
    rdfs:label "GRCh38chr9-37002118-37002118-Region";
    faldo:begin <https://monarchinitiative.org/.well-known/genid/bb870c3d7606a3e0fc3c> ;
    faldo:end <https://monarchinitiative.org/.well-known/genid/bb870c3d7606a3e0fc3c> .

<https://monarchinitiative.org/.well-known/genid/b21985847fe0774084eb> a faldo:Position ;
    rdfs:label "GRCh38chr9-36997420";
    faldo:position 36997420 ;
    faldo:reference OBO:CHR_GRCh38chr9 .

<https://monarchinitiative.org/.well-known/genid/b55051762f8d5a3dbeb5> a faldo:Position ;
    rdfs:label "GRCh38chr9-36998996";
    faldo:position 36998996 ;
    faldo:reference CHR:GRCh38chr9 .

<https://monarchinitiative.org/.well-known/genid/b5d61dbc7958a979d046> a faldo:Position ;
    rdfs:label "GRCh38chr9-37000690";
    faldo:position 37000690 ;
    faldo:reference CHR:GRCh38chr9 .

<https://monarchinitiative.org/.well-known/genid/bb870c3d7606a3e0fc3c> a faldo:Position ;
    rdfs:label "GRCh38chr9-37002118";
    faldo:position 37002118 ;
    faldo:reference CHR:GRCh38chr9 .

HGNC:8619 biolink:category biolink:Gene .
dbSNP:rs1329573 biolink:category biolink:SequenceVariant .
dbSNP:rs3758171 biolink:category biolink:SequenceVariant .
dbSNP:rs3824344 biolink:category biolink:SequenceVariant .
dbSNP:rs7020413 biolink:category biolink:SequenceVariant .
:haplotype_bb627b1f64039b0f751a biolink:category biolink:SequenceVariant .
<https://monarchinitiative.org/.well-known/genid/b096a3e94e32fe23374a>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/b21985847fe0774084eb>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/b25a2da36647bdd71be3>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/b3fad5df82cdfb283329>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/b55051762f8d5a3dbeb5>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/b5d61dbc7958a979d046>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/bb870c3d7606a3e0fc3c>  biolink:category biolink:GenomicSequenceLocalization .
<https://monarchinitiative.org/.well-known/genid/bbb252d9b6cd02e9880a>  biolink:category biolink:GenomicSequenceLocalization .
        """

        # dbg
        LOG.debug(
            "Reference graph: %s",
            self.source.graph.serialize(format="turtle").decode("utf-8"))

        self.assertTrue(self.test_util.test_graph_equality(triples, self.source.graph))


if __name__ == '__main__':
    unittest.main()
