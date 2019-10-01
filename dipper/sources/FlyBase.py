import logging
import re
import csv
import gzip
import os
from ftplib import FTP
from typing import Dict, Tuple, List

from dipper.sources.PostgreSQLSource import PostgreSQLSource
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model
from dipper.models.Reference import Reference
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class FlyBase(PostgreSQLSource):
    """
    This is the [Drosophila Genetics](http://www.flybase.org/) resource,
    from which we process genotype and phenotype data about the fruit fly.

    Here, we connect to their public database and download preprocessed files

    Queries from the relational db
    1. allele-phenotype data: ../../sources/sql/fb/allele_phenotype.sql
    2. gene dbxrefs: ../../resources/sql/fb/gene_xref.sql

    Downloads:
    1. allele_human_disease_model_data_fb_*.tsv.gz - models of disease
    2. species.ab.gz - species prefix mappings
    3. fbal_to_fbgn_fb_*.tsv.gz - allele to gene
    4. fbrf_pmid_pmcid_doi_fb_*.tsv.gz -  flybase ref to pmid

    We connect using the
    [Direct Chado Access](http://gmod.org/wiki/
    Public_Chado_Databases#Direct_Chado_Access)

    When running the whole set,
    it performs best by dumping raw triples using the flag ```--format nt```.

    Note that this script underwent a major revision after commit bd5f555
    in which genotypes, stocks, and environments were removed

    """
    FLYFTP = 'ftp.flybase.net'
    FLYFILES = '/releases/current/precomputed_files'

    queries = {
        'allele_phenotype': {
            'query': '../../resources/sql/fb/allele_phenotype.sql',
            'file': 'allele_phenotype.tsv',
            'columns': [
                'allele_id',
                'pheno_desc',
                'pheno_type',
                'pub_id',
                'pub_title',
                'pmid_id',
            ]
        },
        'gene_xref': {
            'query': '../../resources/sql/fb/gene_xref.sql',
            'file': 'gene_xref.tsv',
            'columns': [
                'gene_id',
                'xref_id',
                'xref_source',
            ]
        }
    }

    files = {
        'disease_model': {
            'file': 'disease_model_annotations.tsv.gz',
            'url':  r'human_disease/disease_model_annotations_.*\.tsv\.gz',
            'columns': [
                'FBgn ID',
                'Gene symbol',
                'HGNC ID',
                'DO qualifier',
                'DO ID',
                'DO term',
                'Allele used in model (FBal ID)',
                'Allele used in model (symbol)',
                'Based on orthology with (HGNC ID)',
                'Based on orthology with (symbol)',
                'Evidence/interacting alleles',
                'Reference (FBrf ID)'
            ]
        },
        'species_map': {
            'file': 'species.ab.gz',
            'url': r'species/species\.ab\.gz',
            'columns': [
                'internal_id',
                'taxgroup',
                'abbreviation',
                'genus',
                'species name',
                'common name',
                'comment',
                'ncbi-taxon-id',
            ]
        },
        'allele_gene': {
            'file': 'fbal_to_fbgn_fb.tsv.gz',
            'url': r'alleles/fbal_to_fbgn_fb_.*\.tsv\.gz',
            'columns': [
                'AlleleID',
                'AlleleSymbol',
                'GeneID',
                'GeneSymbol',
            ]
        },
        'ref_pubmed': {
            'file': 'fbrf_pmid_pmcid_doi_fb.tsv.gz',
            'url': r'references/fbrf_pmid_pmcid_doi_fb_.*\.tsv\.gz',
            'columns': [
                'FBrf',
                'PMID',
                'PMCID',
                'DOI',
                'pub_type',
                'miniref',
                'pmid_added',
            ]
        }
    }

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skolemized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='flybase',
            ingest_title='FlyBase',
            ingest_url='http://www.flybase.org/',
            ingest_logo='source-flybase.png',
            license_url=None,
            data_rights='https://wiki.flybase.org/wiki/FlyBase_Wiki:General_disclaimer',
            file_handle=None
            )

    def fetch(self, is_dl_forced=False):
        """
        Fetch flat files and sql queries

        :param is_dl_forced: force download
        :return: None

        """
        # create the connection details for Flybase
        cxn = {
            'host': 'chado.flybase.org', 'database': 'flybase', 'port': 5432,
            'user': 'flybase', 'password': 'no password'}

        self.dataset.set_ingest_source(
            ''.join(('jdbc:postgresql://', cxn['host'], ':', str(cxn['port']),
                     '/', cxn['database'])), is_object_literal=True)

        # Get data from remote db
        # Each query takes 2 minutes or so
        for query_map in self.queries.values():
            query_fh = open(os.path.join(
                os.path.dirname(__file__), query_map['query']), 'r')
            query = query_fh.read()
            self.fetch_query_from_pgdb(
                query_map['file'], query, None, cxn)

        # Get flat files
        ftp = FTP(FlyBase.FLYFTP)
        ftp.login("anonymous", "info@monarchinitiative.org")

        for src_key, file in self.files.items():
            filename = self._resolve_filename(file['url'], ftp)
            # prepend ftp:// since this gets added to dataset rdf model
            self.files[src_key]['url'] = "ftp://" + self.FLYFTP + filename

        ftp.close()

        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        """
        Parse flybase files and add to graph

        :param limit: number of rows to process
        :return: None

        """
        if limit is not None:
            LOG.info("Only parsing first %d rows of each file", limit)
        LOG.info("Parsing files...")

        self._process_allele_phenotype(limit)
        self._process_allele_gene(limit)
        self._process_disease_model(limit)
        self._process_gene_xref(limit)

        LOG.info("Finished parsing.")
        LOG.info("Loaded %d nodes", len(self.graph))

    def _process_allele_phenotype(self, limit):
        """
        Make allele to phenotype associations using derived_pheno_class
        and derived_pheno_manifest cvterm in the flybase db, an example entry is:

        FBal0257663    @FBcv0000351:lethal@ | @FBcv0000308:female limited@,
                       with @FBal0130657:Scer\GAL4<up>dome-PG14</up>@

        The first term is the phenotype, and all follow up terms are qualifiers,
        self.globaltt['has_qualifier'])

        Our previous approach was to use the genotype id associated with
        FBal0257663/FBal0130657 , however, this required us to create blank
        nodes and was considered unnecessarily granular

        Note that sometimes identifiers do not exist for a term, eg
        @:heat sensitive | tetracycline conditional@

        derived_pheno_class - FBcv terms, these are phenotypes
        derived_pheno_manifest -  Anatomy terms FBbt, we currently
        make phenotype IRI equivalents that end up in UPheno, but
        this is being developed and updated, see
        https://github.com/monarch-initiative/dipper/issues/770

        Adds triples to self.graph

        :param limit: number of rows to process
        :return: None

        """
        model = Model(self.graph)
        src_key = 'allele_phenotype'
        raw = '/'.join((self.rawdir, self.queries[src_key]['file']))
        LOG.info("processing allele phenotype associations")
        col = self.queries[src_key]['columns']

        transgenic_alleles = self._get_foreign_transgenic_alleles()

        # flybase terms - terms we prefix with FlyBase:
        fly_prefixes = ['FBal', 'FBti', 'FBab', 'FBba', 'FBtp']

        # a alphanumeric id followed by a colon then
        # any character but a colon bordered by @s
        term_regex = re.compile(r'@([\w]*):([^:@]*)@')
        id_regex = re.compile(r'([a-zA-Z]+)(\d+)')

        with open(raw, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            row = next(reader)  # headers
            self.check_fileheader(col, row)

            for row in reader:
                allele_id = row[col.index('allele_id')]
                pheno_desc = row[col.index('pheno_desc')]
                pheno_type = row[col.index('pheno_type')]
                pub_id = row[col.index('pub_id')]
                pub_title = row[col.index('pub_title')]
                pmid_id = row[col.index('pmid_id')]

                # Don't get phenotypes for transgenic alleles
                if allele_id in transgenic_alleles:
                    continue

                allele_curie = 'FlyBase:' + allele_id

                terms = re.findall(term_regex, pheno_desc)

                if not terms:
                    LOG.warning('Could not parse pheno description: %s',
                                pheno_desc)
                    continue

                term_ids, term_labels = zip(*terms)
                id_match = re.match(id_regex, term_ids[0])

                if id_match is not None:
                    prefix, reference = id_match.group(1, 2)
                else:
                    raise ValueError("Could not parse id {}".format(term_ids[0]))

                # derived_pheno_class should all start with a FBcv term
                if pheno_type == 'derived_pheno_class' and prefix != 'FBcv':
                    LOG.warning('derived_pheno_class does not '
                                'start with FBcv: %s', pheno_desc)
                    continue

                # Create phenotype curie
                if pheno_type == 'derived_pheno_class':
                    phenotype_curie = prefix + ':' + reference
                elif pheno_type == 'derived_pheno_manifest':
                    if prefix == 'GO':
                        phenotype_curie = prefix + ':' + reference + 'PHENOTYPE'
                        phenotype_label = term_labels[0] + ' phenotype'
                        model.addClassToGraph(phenotype_curie, phenotype_label,
                                              blv.terms.PhenotypicFeature.value)
                    else:
                        phenotype_curie = 'OBO:' + prefix + reference + 'PHENOTYPE'
                else:
                    raise ValueError("Unexpected phenotype type: {}".
                                     format(pheno_type))

                if pmid_id:
                    ref_curie = 'PMID:' + pmid_id
                else:
                    ref_curie = 'FlyBase:' + pub_id
                    reference = Reference(self.graph, ref_curie)
                    reference.setTitle(pub_title)
                    reference.addRefToGraph()

                assoc = G2PAssoc(self.graph, self.name, allele_curie,
                                 phenotype_curie, self.globaltt['has phenotype'])
                assoc.add_source(ref_curie)
                # Associations need to be disambiguated via their qualifiers
                # see http://flybase.org/reports/FBal0207398 as an example
                assoc.set_association_id(assoc.make_association_id(
                    self.name,
                    allele_curie,
                    self.globaltt['has phenotype'],
                    phenotype_curie,
                    term_ids[1:]
                ))
                assoc.add_association_to_graph()
                assoc_id = assoc.get_association_id()

                # add the rest as qualifiers
                for term in term_ids[1:]:
                    if term:
                        # FBal, GO, FBti, FBab ...
                        id_match = re.match(id_regex, term)
                        if id_match is not None:
                            prefix, reference = id_match.group(1, 2)
                            if prefix in fly_prefixes:
                                term_curie = 'FlyBase:' + term
                            else:
                                term_curie = prefix + ':' + reference
                        else:
                            raise ValueError("Could not parse id {}".format(term))

                    else:
                        # There is not an id for a term,
                        # eg @:heat sensitive | tetracycline conditional@
                        continue

                    self.graph.addTriple(
                        assoc_id, self.globaltt['has_qualifier'], term_curie,
                        subject_category=blv.terms.InformationContentEntity.value,
                        object_category=blv.terms.InformationContentEntity.value)

                if limit is not None and reader.line_num > limit:
                    break

    def _species_to_ncbi_tax(self) -> Dict[str, Tuple[str, str]]:
        """
        Generates a dictionary of flybase species prefixes
        and a tuple of taxon group and ncbi taxon curies,
        eg "Hsap":("non-drosophilid eukaryote", "NCBITaxon:9606")

        Note that this file is missing some prefixes as of 2019_02

        :return: Dict with flybase prefixes as keys and a tuple of
                      taxon group and NCBI taxon curies as values

        """
        species_map = {}
        src_key = 'species_map'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("creating map of species prefixes and ncbi taxon curies")

        added_prefixes = {
            'P': ("drosophilid", self.globaltt['Drosophila melanogaster']),
            'Drer': ("non-drosophilid eukaryote", self.globaltt['Danio rerio'])
        }

        col = self.files[src_key]['columns']

        with gzip.open(raw, 'rt') as tsvfile:
            # Delimiter is ' | ' but csv requires 1 char
            reader = csv.reader(tsvfile, delimiter='|')
            # skip first lines
            next(reader)

            row = next(reader)  # headers
            row = [field.strip() for field in row]
            row[0] = row[0][2:]  # uncomment

            self.check_fileheader(col, row)

            next(reader)  # blank line

            for row in reader:
                row = [field.strip() for field in row]
                prefix = row[col.index('abbreviation')]
                tax_group = row[col.index('taxgroup')]
                taxon = row[col.index('ncbi-taxon-id')]
                taxon = taxon.replace("taxon", "NCBITaxon")

                species_map[prefix] = (tax_group, taxon)

            # Add in hard coded prefixes
            for prefix, taxon in added_prefixes.items():
                if prefix not in species_map:
                    species_map[prefix] = taxon

            return species_map

    def _flyref_to_pmid(self) -> Dict[str, str]:
        """
        Generates a dictionary of flybase reference and PMID curie mappings;
        "FBrf0241315":"PMID:30328653"

        :return: Dict with FBrf ids as keys and PMID curies as values
        """
        pub_map = {}
        src_key = 'ref_pubmed'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("creating map of flybase ref ids and pmids")

        col = self.files[src_key]['columns']

        # JR - I've set encoding to latin-1 to fix the UnicodeDecodeError that happens
        # when the default encoding (utf-8) is used. This possibly will break if/when
        # the encoding of this file upstream at Flybase is changed to utf-8. If so,
        # trying setting encoding='utf-8' below
        with gzip.open(raw, 'rt', encoding='latin-1') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            # skip first four lines
            for _ in range(0, 2):
                # file name, note
                next(reader)

            row = next(reader)  # headers
            row[0] = row[0][1:]  # uncomment

            self.check_fileheader(col, row)

            for row in reader:
                # File ends with a final comment
                if ''.join(row).startswith('#'):
                    continue

                pmid = row[col.index('PMID')]
                fly_ref = row[col.index('FBrf')]

                pub_map[fly_ref] = 'PMID:' + pmid

            return pub_map

    def _get_foreign_transgenic_alleles(self) -> List[str]:
        """
        Generates a list of transgenic alleles

        :return: List of FBal ids
        """
        transgenic_alleles = []
        species_map = self._species_to_ncbi_tax()
        src_key = 'allele_gene'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("Getting list of transgenic alleles")

        col = self.files[src_key]['columns']

        with gzip.open(raw, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            # skip first line, version info
            next(reader)
            row = next(reader)  # headers
            # header line starts with a hash and tab ??
            row = row[1:]

            self.check_fileheader(col, row)

            for row in reader:
                allele_id = row[col.index('AlleleID')]
                allele_label = row[col.index('AlleleSymbol')]

                # Add Allele and taxon, get anything that's not drosophila
                allele_prefix = re.findall(r'^(\w*)\\', allele_label)

                if len(allele_prefix) == 1:
                    try:
                        if species_map[allele_prefix[0]][0] != 'drosophilid':
                            transgenic_alleles.append(allele_id)
                    except KeyError:
                        # prefix is not in species prefix file
                        transgenic_alleles.append(allele_id)

        return transgenic_alleles

    def _process_gene_xref(self, limit):
        """
        Make eq axioms between flybase gene ids and ncbi gene and hgnc

        Note that there are a lot of genes in flybase from other organisms
        we make the eq axioms so that they clique merge in our large graph
        (for example, human genes should merge with HGNC)

        Adds triples to self.graph

        :param limit: number of rows to process
        :return: None

        """
        model = Model(self.graph)
        src_key = 'gene_xref'
        raw = '/'.join((self.rawdir, self.queries[src_key]['file']))
        LOG.info("processing gene xrefs")

        col = self.queries[src_key]['columns']

        with open(raw, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            row = next(reader)  # headers
            self.check_fileheader(col, row)

            for row in reader:
                gene_id = row[col.index('gene_id')]
                xref_id = row[col.index('xref_id')]
                xref_source = row[col.index('xref_source')]

                gene_curie = 'FlyBase:' + gene_id
                xref_prefix = None
                if xref_source == 'EntrezGene':
                    xref_prefix = 'NCBIGene'
                elif xref_source == 'HGNC':
                    xref_prefix = 'HGNC'
                xref_curie = xref_prefix + ':' + xref_id

                model.addEquivalentClass(gene_curie, xref_curie,
                                         subject_category=blv.terms.Gene.value,
                                         object_category=blv.terms.Gene.value)

                if limit is not None and reader.line_num > limit:
                    break

    def _process_allele_gene(self, limit):
        """
        Make associations between an allele and a gene
        Adds triples to self.graph

        Approach is to use the label nomenclature and species
        map to determine taxon.  Foreign Transgenes are filtered out.

        :param limit: number of rows to process
        :return: None

        """
        geno = Genotype(self.graph)
        species_map = self._species_to_ncbi_tax()
        src_key = 'allele_gene'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("processing allele to gene")

        col = self.files[src_key]['columns']

        with gzip.open(raw, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            # skip first line, version info
            next(reader)
            row = next(reader)  # headers
            # header line starts with a hash and tab ??
            row = row[1:]

            self.check_fileheader(col, row)

            for row in reader:
                allele_id = row[col.index('AlleleID')]
                allele_label = row[col.index('AlleleSymbol')]
                gene_id = row[col.index('GeneID')]
                gene_label = row[col.index('GeneSymbol')]

                allele_curie = 'FlyBase:' + allele_id
                gene_curie = 'FlyBase:' + gene_id

                # Add Allele and taxon, skip anything that's not drosophila
                allele_prefix = re.findall(r'^(\w*)\\', allele_label)

                if len(allele_prefix) == 1:
                    try:
                        if species_map[allele_prefix[0]][0] == 'drosophilid':
                            geno.addAllele(allele_curie, allele_label)
                            geno.addTaxon(
                                species_map[allele_prefix[0]][1],
                                allele_curie,
                                genopart_category=blv.terms.SequenceVariant.value)
                        else:
                            # If it's a foreign transgenic allele, skip
                            continue
                    except KeyError:
                        LOG.info("%s not in species prefix file", allele_prefix[0])
                        continue

                elif not allele_prefix:
                    geno.addAllele(allele_curie, allele_label)
                    geno.addTaxon(self.globaltt['Drosophila melanogaster'],
                                  allele_curie,
                                  genopart_category=blv.terms.SequenceVariant.value)
                else:
                    raise ValueError("Did not correctly parse allele label {}"
                                     .format(allele_label))
                # Process genes
                gene_prefix = re.findall(r'^(\w*)\\', gene_label)

                if len(gene_prefix) == 1:
                    try:
                        geno.addTaxon(species_map[gene_prefix[0]][1], gene_curie,
                                      genopart_category=blv.terms.Gene.value)

                        if species_map[gene_prefix[0]][0] == 'drosophilid':
                            geno.addGene(gene_curie, gene_label)
                        else:
                            # Don't create labels for non drosophila genes
                            geno.addGene(gene_curie)

                    except KeyError:
                        LOG.info("%s not in species prefix file", gene_prefix[0])
                        geno.addGene(gene_curie)

                elif not gene_prefix:
                    geno.addGene(gene_curie, gene_label)
                    geno.addTaxon(self.globaltt['Drosophila melanogaster'],
                                  allele_curie,
                                  genopart_category=blv.terms.SequenceVariant.value)
                else:
                    raise ValueError("Did not correct parse gene label {}"
                                     .format(gene_label))

                # Connect allele and gene with geno.addAffectedLocus()
                if allele_prefix and gene_prefix:
                    if allele_prefix[0] == gene_prefix[0]:
                        geno.addAffectedLocus(allele_curie, gene_curie)
                    else:
                        raise ValueError(
                            "Found allele and gene with different "
                            "prefixes: {}, {}".format(allele_id, gene_id))
                elif not allele_prefix and gene_prefix:
                    raise ValueError(
                        "Found allele and gene with different "
                        "prefixes: {}, {}".format(allele_id, gene_id))
                else:
                    # Both are melanogaster
                    geno.addAffectedLocus(allele_curie, gene_curie)

                if limit is not None and reader.line_num > limit:
                    break

    def _process_disease_model(self, limit):
        """
        Make associations between a disease and fly alleles
        Adds triples to self.graph

        Pulls FBrf to pmid eqs from _flyref_to_pmid
        for 2019_02 maps all but FlyBase:FBrf0211649

        Right now only model_of is processed from DOID_qualifier
        As of 2019_02 release there are also:
        ameliorates
        exacerbates
        DOES NOT ameliorate
        DOES NOT exacerbate
        DOES NOT model
        DOID_qualifier

        :param limit: number of rows to process
        :return: None

        """
        graph = self.graph
        pub_map = self._flyref_to_pmid()
        src_key = 'disease_model'
        raw = '/'.join((self.rawdir, self.files[src_key]['file']))
        LOG.info("processing disease models")
        col = self.files[src_key]['columns']

        transgenic_alleles = self._get_foreign_transgenic_alleles()

        with gzip.open(raw, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            # skip first four lines
            for _ in range(0, 4):
                # file name, generated datetime, db info, blank line
                next(reader)

            row = next(reader)  # headers
            row[0] = row[0][3:]  # uncomment

            self.check_fileheader(col, row)

            for row in reader:
                # File ends with a blank line and a final comment
                if ''.join(row).startswith('#') or not ''.join(row).strip():
                    continue

                allele_id = row[col.index('Allele used in model (FBal ID)')]
                flybase_ref = row[col.index('Reference (FBrf ID)')]
                evidence_or_allele = row[col.index('Evidence/interacting alleles')]
                doid_id = row[col.index('DO ID')]
                doid_qualifier = row[col.index('DO qualifier')]

                # Skip any foreign transgenic alleles
                if allele_id in transgenic_alleles:
                    continue

                # Rows without an allele id may be "potential models through orthology"
                # We skip these because we can already make that inference
                # http://flybase.org/reports/FBgn0000022
                if not allele_id:
                    continue

                allele_curie = 'FlyBase:' + allele_id
                if doid_qualifier == 'model of':
                    relation = self.globaltt['is model of']
                else:
                    # amelorates, exacerbates, and DOES NOT *
                    continue

                assoc = G2PAssoc(graph, self.name, allele_curie, doid_id, relation,
                                 phenotype_category=blv.terms.Disease.value)
                if flybase_ref != '':
                    ref_curie = None
                    try:
                        ref_curie = pub_map[flybase_ref]
                    except KeyError:
                        ref_curie = 'FlyBase:' + flybase_ref

                    assoc.add_source(ref_curie)
                if evidence_or_allele == 'inferred from mutant phenotype':
                    evidence_id = self.globaltt['mutant phenotype evidence']
                    assoc.add_evidence(evidence_id)
                else:
                    assoc.set_description(evidence_or_allele)

                assoc.add_association_to_graph()

                if limit is not None and reader.line_num > limit:
                    break

    @staticmethod
    def _resolve_filename(filename: str, ftp: FTP) -> str:
        """
        Resolve a file name from ftp server given a regex
        :return: str, file path on ftp server
        """

        # Represent file path as a list of directories
        dir_path = FlyBase.FLYFILES.split('/')
        # Also split on the filename to get prepending dirs
        file_path = filename.split('/')
        file_regex = file_path.pop()
        working_dir = "/".join(dir_path + file_path)

        LOG.info('Looking for remote files in %s', working_dir)

        ftp.cwd(working_dir)
        remote_files = ftp.nlst()

        files_to_download = [dnload for dnload in remote_files
                             if re.match(file_regex, dnload)]

        if len(files_to_download) > 1:
            raise ValueError("Could not resolve filename from regex, "
                             "too many matches for {}, matched: {}"
                             .format(file_regex, files_to_download))
        if not files_to_download:
            raise ValueError("Could not resolve filename from regex, "
                             "no matches for {}".format(file_regex))

        return working_dir + '/' + files_to_download[0]
