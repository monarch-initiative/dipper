import csv
import re
import logging
import gzip
import io
import sys
import os

import yaml

from dipper.sources.ZFIN import ZFIN
from dipper.sources.WormBase import WormBase
from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Reference import Reference
from dipper.models.Model import Model
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)
# get gene annotation from current.geneontology.com,
# which is the last official release (but not the bleeding edge)
GOGA = 'http://current.geneontology.org/annotations'
FTPEBI = 'ftp://ftp.uniprot.org/pub/databases/'     # best for North America
UPCRKB = 'uniprot/current_release/knowledgebase/'


class GeneOntology(Source):
    """
    This is the parser for the
    [Gene Ontology Annotations](http://www.geneontology.org),
    from which we process gene-process/function/subcellular
    location associations.

    We generate the GO graph to include the following information:
    * genes
    * gene-process
    * gene-function
    * gene-location

    We process only a subset of the organisms:

    Status: IN PROGRESS / INCOMPLETE

    """
    gaf_columns = [  # GAF2.1 files contain the following columns:
        'DB',
        'DB_Object_ID',
        'DB_Object_Symbol',
        'Qualifier',
        'GO_ID',
        'DB:Reference',
        'Evidence Code',
        'With (or) From',
        'Aspect',
        'DB_Object_Name',
        'DB_Object_Synonym',
        'DB_Object_Type',
        'Taxon and Interacting taxon',
        'Date',
        'Assigned_By',
        'Annotation_Extension',
        'Gene_Product_Form_ID'
    ]

    files = {
        '9615': {  # Canis lupus familiaris
            'file': 'goa_dog.gaf.gz',
            'url': GOGA + '/goa_dog.gaf.gz',
            'columnns': gaf_columns
        },
        '7227': {  # Drosophila melanogaster
            'file': 'fb.gaf.gz',
            'url': GOGA + '/fb.gaf.gz',
            'columnns': gaf_columns
        },
        '7955': {  # Danio rerio
            'file': 'zfin.gaf.gz',
            'url': GOGA + '/zfin.gaf.gz',
            'columnns': gaf_columns
        },
        '10090': {  # Mus musculus
            'file': 'mgi.gaf.gz',
            'url': GOGA + '/mgi.gaf.gz',
            'columnns': gaf_columns
        },
        '10116': {  # Rattus norvegicus
            'file': 'rgd.gaf.gz',
            'url': GOGA + '/rgd.gaf.gz',
            'columnns': gaf_columns
        },
        '6239': {  # Caenorhabditis elegans
            'file': 'wb.gaf.gz',
            'url': GOGA + '/wb.gaf.gz',
            'columnns': gaf_columns
        },
        '9823': {  # Sus scrofa
            'file': 'goa_pig.gaf.gz',
            'url': GOGA + '/goa_pig.gaf.gz',
            'columnns': gaf_columns
        },
        '9031': {  # Gallus gallus
            'file': 'goa_chicken.gaf.gz',
            'url': GOGA + '/goa_chicken.gaf.gz',
            'columnns': gaf_columns
        },
        '9606': {  # Homo sapiens
            'file': 'goa_human.gaf.gz',
            'url': GOGA + '/goa_human.gaf.gz',
            'columnns': gaf_columns
        },
        '9913': {  # Bos taurus
            'file': 'goa_cow.gaf.gz',
            'url': GOGA + '/goa_cow.gaf.gz',
            'columnns': gaf_columns
        },
        '559292': {  # Saccharomyces cerevisiae   4932
            'file': 'sgd.gaf.gz',
            'url': GOGA + '/sgd.gaf.gz',
            'columnns': gaf_columns
        },
        '4896': {  # Schizosaccharomyces pombe  (yeast)
            'file': 'pombase.gaf.gz',
            'url': GOGA + '/pombase.gaf.gz',
            'columnns': gaf_columns
        },
        '5782': {  # Dictyostelium (slime mold genus)
            'file': 'dictibase.gaf.gz',
            'url': GOGA + '/dictybase.gaf.gz',
            'columnns': gaf_columns
        },
        '5052':  {  # Aspergillus  (fungi)  http://www.aspergillusgenome.org/
            'file': 'aspgd.gaf.gz',
            'url': GOGA + '/aspgd.gaf.gz',
            'columnns': gaf_columns
        },

        # consider this after most others - should this be part of GO?
        # 'multispecies': {
        #   'file': 'gene_association.goa_uniprot.gz',
        #   'url': FTPEBI + 'GO/goa/UNIPROT/gene_association.goa_uniprot.gz'},

        'go-references': {
            'file': 'GO.references',
            # Quoth the header of this file: "This file is DEPRECATED.
            # Please see go-refs.json relative to this location"
            # (http://current.geneontology.org/metadata/go-refs.json)
            'url': 'http://www.geneontology.org/doc/GO.references'
        },
        'id-map': {  # 5GB mapping file takes 6 hours to DL ... maps UniProt to Ensembl
            'file': 'idmapping_selected.tab.gz',
            'url':  FTPEBI + UPCRKB + 'idmapping/idmapping_selected.tab.gz',
            # ftp://ftp.uniprot.org
            # /pub/databases/uniprot/current_release/knowledgebase/idmapping/README
            'columns': [
                'UniProtKB-AC',
                'UniProtKB-ID',
                'GeneID (EntrezGene)',
                'RefSeq',
                'GI',
                'PDB',
                'GO',
                'UniRef100',
                'UniRef90',
                'UniRef50',
                'UniParc',
                'PIR',
                'NCBI-taxon',
                'MIM',
                'UniGene',
                'PubMed',
                'EMBL',
                'EMBL-CDS',
                'Ensembl',
                'Ensembl_TRS',
                'Ensembl_PRO',
                'Additional PubMed'
            ]
        }
    }
    # consider moving the go-ref and id-map above to here in map_files
    map_files = {
        'eco_map': 'http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt',
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(
            graph_type,
            are_bnodes_skolemized,
            'go',
            ingest_title='Gene Ontology',
            ingest_url='http://www.geneontology.org',
            license_url=None,
            data_rights='http://geneontology.org/page/use-and-license'
            # file_handle=None
        )
        self.test_ids = []
        # note: dipper-etl defaults tax_ids to '9606'
        # note: sorting tax_ids for stable digest
        if tax_ids is not None and [] != set(tax_ids).difference(['9606']):
            LOG.info('Have %s given as taxon to ingest', str(tax_ids))
            self.tax_ids = sorted([str(x) for x in tax_ids])
            nottax = set(tax_ids) - set(self.files.keys())
            if nottax:
                LOG.error('Cant process taxon number(s):\t%s', str(nottax))
                self.tax_ids = list(set(self.tax_ids) - nottax)
        else:
            self.tax_ids = sorted(['9606', '10090', '7955'])

        LOG.info("Filtering to the following taxa: %s", self.tax_ids)

        # moving this from process_gaf() to avoid repeating this for each
        # file to be processed.
        if '7955' in self.tax_ids:
            self.zfin = ZFIN(self.graph_type, self.are_bnodes_skized)
        if '6239' in self.tax_ids:
            self.wbase = WormBase(self.graph_type, self.are_bnodes_skized)

        if 'gene' not in self.all_test_ids:
            LOG.warning("not configured with gene test ids.")
        else:
            self.test_ids = self.all_test_ids['gene']

        # build the id map for mapping uniprot ids to genes ... ONCE
        self.uniprot_entrez_id_map = self.get_uniprot_entrez_id_map()
        self.eco_map = self.get_eco_map(self.map_files['eco_map'])

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        if limit is not None:
            LOG.info("Only parsing first %s rows of each file", limit)
        LOG.info("Parsing files...")

        if self.test_only:
            self.test_mode = True

        for txid_num in list(set(self.files).intersection(self.tax_ids)):
            gaffile = '/'.join((self.rawdir, self.files[txid_num]['file']))
            self.process_gaf(gaffile, limit, self.uniprot_entrez_id_map, self.eco_map)

        LOG.info("Finished parsing.")

    def process_gaf(self, gaffile, limit, id_map=None, eco_map=None):

        if self.test_mode:
            graph = self.testgraph
        else:
            graph = self.graph

        model = Model(graph)
        geno = Genotype(graph)
        LOG.info("Processing Gene Associations from %s", gaffile)
        uniprot_hit = 0
        uniprot_miss = 0
        col = self.gaf_columns

        with gzip.open(gaffile, 'rb') as csvfile:
            reader = csv.reader(
                io.TextIOWrapper(csvfile, newline=""), delimiter='\t', quotechar='\"')
            for row in reader:
                # comments start with exclamation
                if row[0][0] == '!':
                    continue

                if len(row) != len(col):
                    LOG.error(
                        "Wrong number of columns %i, expected ... got:\n\t%s",
                        len(col), row)
                    exit(1)

                dbase = row[col.index('DB')].strip()
                gene_num = row[col.index('DB_Object_ID')].strip()
                gene_symbol = row[col.index('DB_Object_Symbol')].strip()
                qualifier = row[col.index('Qualifier')]
                go_id = row[col.index('GO_ID')].strip()
                ref = row[col.index('DB:Reference')].strip()
                eco_symbol = row[col.index('Evidence Code')].strip()
                with_or_from = row[col.index('With (or) From')]
                aspect = row[col.index('Aspect')].strip()
                gene_name = row[col.index('DB_Object_Name')]
                gene_synonym = row[col.index('DB_Object_Synonym')]
                # object_type = row[col.index('DB_Object_Type')].strip()
                taxon = row[col.index('Taxon and Interacting taxon')].strip()
                # date = row[col.index('Date')].strip()
                # assigned_by = row[col.index('Assigned_By')].strip()
                # annotation_extension = row[col.index('Annotation_Extension')]
                # gene_product_form_id = row[col.index('Gene_Product_Form_ID')]

                # test for required fields
                if '' in [row[:10], row[12]]:
                    LOG.error(
                        "Missing required part of annotation on row %i:\n%s",
                        reader.line_num, str(row[:-4]))
                    continue

                # (Don't) deal with qualifier NOT, contributes_to, colocalizes_with
                if re.search(r'NOT', qualifier):
                    continue

                if dbase in self.localtt:
                    dbase = self.localtt[dbase]
                uniprotid = None
                gene_id = None
                if dbase == 'UniProtKB':
                    if id_map is not None and gene_num in id_map:
                        gene_id = id_map[gene_num]
                        uniprotid = ':'.join((dbase, gene_num))
                        (dbase, gene_num) = gene_id.split(':')
                        uniprot_hit += 1
                    else:
                        # LOG.warning(
                        #   "UniProt id %s is without a 1:1 mapping to entrez/ensembl",
                        #    gene_num)
                        uniprot_miss += 1
                        continue
                else:
                    gene_num = gene_num.split(':')[-1]  # last
                    gene_id = ':'.join((dbase, gene_num))

                if self.test_mode and gene_id[:9] != 'NCBIGene:' and\
                        gene_num not in self.test_ids:
                    continue

                model.addClassToGraph(gene_id, gene_symbol,
                                      class_category=blv.Gene.value)
                if gene_name != '':
                    model.addDescription(gene_id, gene_name,
                                         subject_category=blv.Gene.value)
                if gene_synonym != '':
                    for syn in re.split(r'\|', gene_synonym):
                        syn = syn.strip()
                        if syn[:10] == 'UniProtKB:':
                            model.addTriple(
                                gene_id, self.globaltt['has gene product'], syn,
                                subject_category=blv.Gene.value,
                                object_category=blv.Gene.value)
                        elif re.fullmatch(graph.curie_regexp, syn) is not None:
                            LOG.warning(
                                'possible curie "%s" as a literal synomym for %s',
                                syn, gene_id)
                            model.addSynonym(gene_id, syn,
                                             class_category=blv.Gene.value,
                                             synonym_type_category=blv.Gene.value)
                        else:
                            model.addSynonym(gene_id, syn,
                                             class_category=blv.Gene.value)

                for txid in taxon.split('|'):
                    tax_curie = re.sub(r'taxon:', 'NCBITaxon:', txid)
                    geno.addTaxon(tax_curie, gene_id, genopart_category=blv.Gene.value)

                assoc = Assoc(graph, self.name,
                              subject_category=blv.Gene.value,
                              object_category=blv.OntologyClass.value)
                assoc.set_subject(gene_id)
                assoc.set_object(go_id)

                try:
                    eco_id = eco_map[eco_symbol]
                    assoc.add_evidence(eco_id)
                except KeyError:
                    LOG.error("Evidence code (%s) not mapped", eco_symbol)

                refs = re.split(r'\|', ref)
                for ref in refs:
                    ref = ref.strip()
                    if ref != '':
                        prefix = ref.split(':')[0]  # sidestep 'MGI:MGI:'
                        if prefix in self.localtt:
                            prefix = self.localtt[prefix]
                        ref = ':'.join((prefix, ref.split(':')[-1]))
                        refg = Reference(graph, ref)
                        if prefix == 'PMID':
                            ref_type = self.globaltt['journal article']
                            refg.setType(ref_type)
                        refg.addRefToGraph()
                        assoc.add_source(ref)

                # TODO add the source of the annotations from assigned by?

                rel = self.resolve(aspect, mandatory=False)
                if rel is not None and aspect == rel:
                    if aspect == 'F' and re.search(r'contributes_to', qualifier):
                        assoc.set_relationship(self.globaltt['contributes to'])
                    else:
                        LOG.error(
                            "Aspect: %s with qualifier: %s  is not recognized",
                            aspect, qualifier)
                elif rel is not None:
                    assoc.set_relationship(rel)
                    assoc.add_association_to_graph()
                else:
                    LOG.warning("No predicate for association \n%s\n", str(assoc))

                if uniprotid is not None:
                    assoc.set_description('Mapped from ' + uniprotid)
                # object_type should be one of:
                # protein_complex; protein; transcript; ncRNA; rRNA; tRNA;
                # snRNA; snoRNA; any subtype of ncRNA in the Sequence Ontology.
                # If the precise product type is unknown,
                # gene_product should be used
                ########################################################################

                # Derive G2P Associations from IMP annotations
                # in version 2.1 Pipe will indicate 'OR'
                # and Comma will indicate 'AND'.
                # in version 2.0, multiple values are separated by pipes
                # where the pipe has been used to mean 'AND'
                if eco_symbol == 'IMP' and with_or_from != '':
                    withitems = with_or_from.split('|')
                    phenotypeid = go_id + 'PHENOTYPE'
                    # create phenotype associations
                    for itm in withitems:
                        if itm == '' or re.match(
                                r'(UniProtKB|WBPhenotype|InterPro|HGNC)', itm):
                            LOG.warning(
                                "Skipping  %s from or with %s", uniprotid, itm)
                            continue
                        itm = re.sub(r'MGI\:MGI\:', 'MGI:', itm)
                        itm = re.sub(r'WB:', 'WormBase:', itm)

                        # for worms and fish, they might give a RNAi or MORPH
                        # in these cases make a reagent-targeted gene
                        if re.search('MRPHLNO|CRISPR|TALEN', itm):
                            targeted_gene_id = self.zfin.make_targeted_gene_id(
                                gene_id, itm)
                            geno.addReagentTargetedGene(itm, gene_id, targeted_gene_id)
                            # TODO PYLINT why is this needed?
                            # Redefinition of assoc type from
                            # dipper.models.assoc.Association.Assoc to
                            # dipper.models.assoc.G2PAssoc.G2PAssoc
                            assoc = G2PAssoc(
                                graph, self.name, targeted_gene_id, phenotypeid)
                        elif re.search(r'WBRNAi', itm):
                            targeted_gene_id = self.wbase.make_reagent_targeted_gene_id(
                                gene_id, itm)
                            geno.addReagentTargetedGene(itm, gene_id, targeted_gene_id)
                            assoc = G2PAssoc(
                                graph, self.name, targeted_gene_id, phenotypeid)
                        else:
                            assoc = G2PAssoc(graph, self.name, itm, phenotypeid)
                        for ref in refs:
                            ref = ref.strip()
                            if ref != '':
                                prefix = ref.split(':')[0]
                                if prefix in self.localtt:
                                    prefix = self.localtt[prefix]
                                ref = ':'.join((prefix, ref.split(':')[-1]))
                                assoc.add_source(ref)
                                # experimental phenotypic evidence
                                assoc.add_evidence(
                                    self.globaltt['experimental phenotypic evidence'])
                        assoc.add_association_to_graph()
                        # TODO should the G2PAssoc be the evidence for the GO assoc?

                if not self.test_mode and limit is not None and \
                        reader.line_num > limit:
                    break
            uniprot_tot = (uniprot_hit + uniprot_miss)
            uniprot_per = 0.0
            if uniprot_tot != 0:
                uniprot_per = 100.0 * uniprot_hit / uniprot_tot
            LOG.info(
                "Uniprot: %.2f%% of %i benefited from the 1/4 day id mapping download",
                uniprot_per, uniprot_tot)

    def get_uniprot_entrez_id_map(self):
        src_key = 'id-map'
        taxon_digest = GraphUtils.digest_id(str(self.tax_ids))
        id_map = {}
        smallfile = '/'.join((self.rawdir, 'id_map_' + taxon_digest + '.yaml'))
        bigfile = '/'.join((self.rawdir, self.files[src_key]['file']))

        # if processed smallfile exists and is newer than bigfile then use it instesd
        if os.path.isfile(smallfile) and \
                os.path.getctime(smallfile) > os.path.getctime(bigfile):
            LOG.info("Using the cheap mapping file %s", smallfile)
            with open(smallfile, 'r') as yamlreader:
                id_map = yaml.safe_load(yamlreader)
        else:
            LOG.info(
                "Expensive Mapping from Uniprot IDs to Entrez/ENSEMBL gene ids for %s",
                self.tax_ids)
            self.fetch_from_url(self.files[src_key]['url'], bigfile)
            col = self.files[src_key]['columns']
            ummapped_uniprot = 0
            with gzip.open(bigfile, 'rb') as csvfile:
                csv.field_size_limit(sys.maxsize)
                reader = csv.reader(  # warning this file is over 10GB unzipped
                    io.TextIOWrapper(csvfile, newline=""),
                    delimiter='\t', quotechar='\"')
                for row in reader:
                    uniprotkb_ac = row[col.index('UniProtKB-AC')].strip()
                    # uniprotkb_id = row[col.index('UniProtKB-ID')]
                    geneid = row[col.index('GeneID (EntrezGene)')].strip()
                    # refseq = row[col.index('RefSeq')]
                    # gi = row[col.index('GI')]
                    # pdb = row[col.index('PDB')]
                    # go = row[col.index('GO')]
                    # uniref100 = row[col.index('UniRef100')]
                    # unifref90 = row[col.index('UniRef90')]
                    # uniref50 = row[col.index('UniRef50')]
                    # uniparc = row[col.index('UniParc')]
                    # pir = row[col.index('PIR')]
                    ncbitaxon = row[col.index('NCBI-taxon')].strip()
                    # mim = row[col.index('MIM')]
                    # unigene = row[col.index('UniGene')]
                    # pubmed = row[col.index('PubMed')]
                    # embl = row[col.index('EMBL')]
                    # embl_cds = row[col.index('EMBL-CDS')]
                    ensembl = row[col.index('Ensembl')].strip()
                    # ensembl_trs = row[col.index('Ensembl_TRS')]
                    # ensembl_pro = row[col.index('Ensembl_PRO')]
                    # other_pubmed = row[col.index('Additional PubMed')]

                    if ncbitaxon not in self.tax_ids:
                        continue

                    # neither empty nor a list
                    if geneid != '' and ';' not in geneid:
                        id_map[uniprotkb_ac] = 'NCBIGene:' + geneid
                    elif ensembl != '' and ';' not in ensembl:
                        id_map[uniprotkb_ac] = 'ENSEMBL:' + ensembl
                    else:
                        ummapped_uniprot += 1

            LOG.info("Writing id_map out as %s", smallfile)
            with open(smallfile, 'w') as yamlwriter:
                yaml.dump(id_map, yamlwriter)
            LOG.warning('Did not find 1:1 gene IDs for %i uniprots', ummapped_uniprot)
        LOG.info(
            "Acquired %i 1:1 uniprot to [entrez|ensembl] mappings", len(id_map.keys()))

        return id_map

    def getTestSuite(self):
        import unittest
        from tests.test_geneontology import GeneOntologyTestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(GeneOntologyTestCase)

        return test_suite
