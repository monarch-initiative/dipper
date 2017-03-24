from dipper.sources.Source import Source
from dipper.models.Genotype import Genotype
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
from dipper.utils.DipperUtil import DipperUtil
from rdflib import RDFS, BNode
from dipper import config
import requests
import logging
import csv
import re
import os


logger = logging.getLogger(__name__)


class UDP(Source):
    """
    The National Institutes of Health (NIH) Undiagnosed Diseases Program (UDP)
    is part of the Undiagnosed Disease Network (UDN),
    an NIH Common Fund initiative that focuses on the most puzzling medical cases
    referred to the NIH Clinical Center in Bethesda, Maryland.
    from https://www.genome.gov/27544402/the-undiagnosed-diseases-program/

    Data is available by request for access via the NHGRI collaboration server:
    https://udplims-collab.nhgri.nih.gov/api

    Note the fetcher requires credentials for the UDP collaboration server
    Credentials are added via a config file, config.json, in the following format
    {
      "dbauth" : {
        "udp": {
          "user": "foo"
          "password": "bar"
      }
    }
    See dipper/config.py for more information

    Output of fetcher:
    udp_variants.tsv
    'Patient', 'Family', 'Chr', 'Build', 'Chromosome Position',
    'Reference Allele', 'Variant Allele', 'Parent of origin',
    'Allele Type', 'Mutation Type', 'Gene', 'Transcript', 'Original Amino Acid',
    'Variant Amino Acid', 'Amino Acid Change', 'Segregates with',
    'Position', 'Exon', 'Inheritance model', 'Zygosity', 'dbSNP ID', '1K Frequency',
    'Number of Alleles'

    udp_phenotypes.tsv
    'Patient', 'HPID', 'Present'

    The script also utilizes two mapping files
    udp_gene_map.tsv -  generated from scripts/fetch-gene-ids.py,
                        gene symbols from udp_variants
    udp_chr_rs.tsv - rsid(s) per coordinate greped from hg19 dbsnp file,
                     then disambiguated with eutils, see scripts/dbsnp/dbsnp.py

    """
    files = {
        'patient_phenotypes': {
            'file': 'udp_phenotypes.tsv'
        },
        'patient_variants': {
            'file': 'udp_variants.tsv'
        }
    }
    map_files = {
        'patient_ids': '../../resources/udp/patient_ids.yaml',
        'dbsnp_map': '../../resources/udp/udp_chr_rs.tsv',
        'gene_coord_map': '../../resources/udp/gene_coordinates.tsv'
    }

    UDP_SERVER = 'https://udplims-collab.nhgri.nih.gov/api'

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'udp')
        self.dataset = Dataset(
            'udp', 'UDP', 'https://rarediseases.info.nih.gov/')
        if graph_type != 'rdf_graph':
            raise ValueError("UDP requires a rdf_graph")

    def fetch(self, is_dl_forced=True):
        """
        Fetches data from udp collaboration server,
        see top level comments for class for more information
        :return:
        """

        username = config.get_config()['dbauth']['udp']['user']
        password = config.get_config()['dbauth']['udp']['password']

        credentials = (username, password)

        # Get patient map file:
        patient_id_map = self.open_and_parse_yaml(self.map_files['patient_ids'])
        udp_internal_ids = patient_id_map.keys()

        phenotype_fields = ['Patient', 'HPID', 'Present']

        # Get phenotype ids for each patient
        phenotype_params = {'method': 'search_subjects',
                            'subject_type': 'Phenotype',
                            'search_mode': 'DEEP',
                            'fields': 'Patient',
                            'conditions': 'equals',
                            'values': ','.join(udp_internal_ids),
                            'user_fields': ','.join(phenotype_fields)
                            }

        prioritized_variants = ['Patient', 'Gene', 'Chromosome Position',
                                'Variant Allele', 'Transcript']

        prioritized_params = {'method': 'search_subjects',
                              'subject_type': 'Variant Prioritization',
                              'search_mode': 'DEEP',
                              'fields': 'Patient',
                              'conditions': 'equals',
                              'values': ','.join(udp_internal_ids),
                              'user_fields': ','.join(prioritized_variants),
                              'format': 'json'}

        variant_fields = ['Patient', 'Family', 'Chr', 'Build', 'Chromosome Position',
                          'Reference Allele', 'Variant Allele', 'Parent of origin',
                          'Allele Type', 'Mutation Type', 'Gene', 'Transcript', 'Original Amino Acid',
                          'Variant Amino Acid', 'Amino Acid Change', 'Segregates with',
                          'Position', 'Exon', 'Inheritance model', 'Zygosity', 'dbSNP ID', '1K Frequency',
                          'Number of Alleles']

        variant_params = {'method': 'search_subjects',
                          'subject_type': 'Exome Analysis Results',
                          'search_mode': 'DEEP',
                          'fields': 'Patient',
                          'conditions': 'equals',
                          'user_fields': ','.join(variant_fields),
                          'format': 'json'}

        pheno_file = open('/'.join((self.rawdir, self.files['patient_phenotypes']['file'])), 'w')
        variant_file = open('/'.join((self.rawdir, self.files['patient_variants']['file'])), 'w')

        pheno_file.write('{0}\n'.format('\t'.join(phenotype_fields)))
        variant_file.write('{0}\n'.format('\t'.join(variant_fields)))

        variant_gene = self._fetch_data_from_udp(
            udp_internal_ids, prioritized_params, prioritized_variants, credentials)

        variant_gene_map = dict()
        for line in variant_gene:
            variant_gene_map.setdefault(line[0], []).append(
                # Try to make a unique value based on gene-pos-variantAlele-transcript
                # TODO make this a dict for readability purposes
                "{0}-{1}-{2}-{3}".format(line[1], line[2], line[3], line[4]))

        variant_info = self._fetch_data_from_udp(
            udp_internal_ids, variant_params, variant_fields, credentials)

        for line in variant_info:
            variant = "{0}-{1}-{2}-{3}".format(line[10], line[4], line[6], line[11])
            if variant in variant_gene_map[line[0]]:
                line[0] = patient_id_map[line[0]]
                line[4] = re.sub(r'\.0$', '', line[4])
                variant_file.write('{0}\n'.format('\t'.join(line)))

        phenotype_info = self._fetch_data_from_udp(
            udp_internal_ids, phenotype_params, phenotype_fields, credentials)

        for line in phenotype_info:
            line[0] = patient_id_map[line[0]]
            pheno_file.write('{0}\n'.format('\t'.join(line)))

        variant_file.close()
        pheno_file.close()

        return

    def parse(self, limit=None):
        """
        Override Source.parse()
        Args:
            :param limit (int, optional) limit the number of rows processed
        Returns:
            :return None
        """
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        phenotype_file = '/'.join((self.rawdir, self.files['patient_phenotypes']['file']))
        variant_file = '/'.join((self.rawdir, self.files['patient_variants']['file']))

        phenotype_file_handler = open(phenotype_file, 'r')
        variant_file_handler = open(variant_file, 'r')

        self._parse_patient_phenotypes(phenotype_file_handler, limit)
        self._parse_patient_variants(variant_file_handler)

        return

    def _fetch_data_from_udp(self, patient_ids, params, fields, credentials):
        data = []
        for patient in patient_ids:
            params['values'] = patient
            result_count = 1
            params['start'] = 0
            params['limit'] = 100
            while result_count > 0:

                logger.debug("processing {0} lines starting at {1} for patient {2}"
                             .format(params['limit'], params['start'], patient))

                req = requests.get(self.UDP_SERVER, params=params, auth=credentials)
                results = req.json()
                result_count = len(results['Subjects'])
                for res in results['Subjects']:
                    line = [res[field] for field in fields]
                    data.append(line)

                params['start'] = params['start'] + params['limit']

        return data

    def _parse_patient_variants(self, file):
        """
        :param file: file handler
        :return:
        """
        patient_var_map = self._convert_variant_file_to_dict(file)
        gene_coordinate_map = self._parse_gene_coordinates(self.map_files['gene_coord_map'])
        rs_map = self._parse_rs_map_file(self.map_files['dbsnp_map'])

        genotype = Genotype(self.graph)
        model = Model(self.graph)

        self._add_variant_gene_relationship(patient_var_map, gene_coordinate_map)

        for patient in patient_var_map:
            patient_curie = ':{0}'.format(patient)
            # make intrinsic genotype for each patient
            intrinsic_geno_bnode = self.make_id("{0}-intrinsic-genotype".format(patient), "_")
            genotype_label = "{0} genotype".format(patient)
            genotype.addGenotype(intrinsic_geno_bnode,
                                 genotype_label,
                                 genotype.genoparts['intrinsic_genotype'])

            self.graph.addTriple(patient_curie,
                                 genotype.object_properties['has_genotype'],
                                 intrinsic_geno_bnode)
            for variant_id, variant in patient_var_map[patient].items():
                build = variant['build']
                chromosome = variant['chromosome']
                position = variant['position']
                reference_allele = variant['reference_allele']
                variant_allele = variant['variant_allele']
                genes_of_interest = variant['genes_of_interest']
                rs_id = variant['rs_id']

                variant_label = ''
                variant_bnode\
                    = self.make_id("{0}".format(variant_id), "_")

                # maybe should have these look like the elif statements below
                if position and reference_allele and variant_allele:
                    variant_label = \
                        self._build_variant_label(build, chromosome,
                                                  position, reference_allele,
                                                  variant_allele, genes_of_interest)
                elif not position and reference_allele and variant_allele \
                        and len(genes_of_interest) == 1:

                        variant_label = \
                            self._build_variant_label(build, chromosome,
                                                      position, reference_allele,
                                                      variant_allele, genes_of_interest)
                elif position and (not reference_allele or not variant_allele) \
                        and len(genes_of_interest) == 1:

                        variant_label = "{0}{1}({2}):g.{3}".format(
                            build, chromosome, genes_of_interest[0], position)
                elif len(genes_of_interest) == 1:
                    variant_label = 'variant of interest' \
                                    ' in {0} gene of patient' \
                                    ' {1}'.format(genes_of_interest[0], patient)
                else:
                    variant_label = 'variant of interest' \
                                    ' in patient {0}'.format(patient)

                genotype.addSequenceAlteration(variant_bnode, None)
                # check if it we have built the label
                # in _add_variant_gene_relationship()
                labels = self.graph.objects(
                    BNode(re.sub(r'^_:', '', variant_bnode, 1)), RDFS['label'])

                label_list = list(labels)

                if len(label_list) == 0:
                    model.addLabel(variant_bnode, variant_label)

                self.graph.addTriple(variant_bnode,
                                     model.object_properties['in_taxon'],
                                     'NCBITaxon:9606')
                self.graph.addTriple(intrinsic_geno_bnode,
                                     genotype.object_properties['has_alternate_part'],
                                     variant_bnode)
                if rs_id:
                    dbsnp_curie = 'dbSNP:{0}'.format(rs_id)
                    model.addSameIndividual(variant_bnode, dbsnp_curie)

        self._add_variant_sameas_relationships(patient_var_map, rs_map)
        return

    def _add_variant_gene_relationship(self, patient_var_map, gene_coordinate_map):
        """
        Right now it is unclear the best approach on how to connect
        variants to genes.  In most cases has_affected_locus/GENO:0000418
        is accurate; however, there are cases where a variant is in the intron
        on one gene and is purported to causally affect another gene down or
        upstream.  In these cases we must first disambiguate which gene
        is the affected locus, and which gene(s) are predicated to be
        causully influenced by (RO:0002566)

        UPDATE 8-30: In the latest dataset we no longer have 1-many mappings
        between variants and genes, but leaving this here in case we see
        these in the future

        The logic followed here is:
        if mutation type contains downstream/upstream and more than one
        gene of interest, investigate coordinates of all genes to
        see if we can disambiguate which genes are which
        :return: None
        """
        genotype = Genotype(self.graph)
        dipper_util = DipperUtil()
        model = Model(self.graph)
        # Note this could be compressed in someway to remove one level of for looping
        for patient in patient_var_map:
            for variant_id, variant in patient_var_map[patient].items():
                variant_bnode = self.make_id("{0}".format(variant_id), "_")
                genes_of_interest = variant['genes_of_interest']
                if len(genes_of_interest) == 1:
                    # Assume variant is variant allele of gene
                    gene = genes_of_interest[0]
                    gene_id = dipper_util.get_ncbi_id_from_symbol(gene)
                    self._add_gene_to_graph(
                        gene, variant_bnode, gene_id,
                        genotype.object_properties['has_affected_locus'])

                elif re.search(r'upstream|downstream',
                               variant['type'],
                               flags=re.I):
                    # Attempt to disambiguate
                    ref_gene = []
                    up_down_gene = []
                    unmatched_genes = []
                    for gene in variant['genes_of_interest']:
                        if gene_id \
                                and gene_id != '' \
                                and gene_id in gene_coordinate_map:
                            if gene_coordinate_map[gene_id]['start'] \
                                    <= variant['position']\
                                    <= gene_coordinate_map[gene_id]['end']:
                                gene_info = {
                                    'symbol': gene,
                                    'strand': gene_coordinate_map[gene_id]['strand']
                                }
                                ref_gene.append(gene_info)
                            else:
                                up_down_gene.append(gene)
                        else:
                            unmatched_genes.append(gene)
                    if len(ref_gene) == 1:
                        self._add_gene_to_graph(
                            ref_gene[0]['symbol'], variant_bnode, gene_id,
                            genotype.object_properties['has_affected_locus'])

                        # update label with gene
                        gene_list = [ref_gene[0]['symbol']]  # build label expects list
                        variant_label = self._build_variant_label(
                            variant['build'], variant['chromosome'],
                            variant['position'], variant['reference_allele'],
                            variant['variant_allele'], gene_list)
                        model.addLabel(variant_bnode, variant_label)

                    # In some cases there are multiple instances
                    # of same gene from dupe rows in the source
                    # Credit http://stackoverflow.com/a/3844832
                    elif len(ref_gene) > 0 and ref_gene[1:] == ref_gene[:-1]:
                        self._add_gene_to_graph(
                            ref_gene[0]['symbol'], variant_bnode, gene_id,
                            genotype.object_properties['has_affected_locus'])

                        # build label function expects list
                        gene_list = [ref_gene[0]['symbol']]
                        variant_label = self._build_variant_label(
                            variant['build'], variant['chromosome'],
                            variant['position'], variant['reference_allele'],
                            variant['variant_allele'], gene_list)
                        model.addLabel(variant_bnode, variant_label)

                    # Check if reference genes are on different strands
                    elif len(ref_gene) == 2:
                        strands = [st['strand'] for st in ref_gene]
                        if "minus" in strands and "plus" in strands:
                            for r_gene in ref_gene:
                                self._add_gene_to_graph(
                                    r_gene['symbol'], variant_bnode, gene_id,
                                    genotype.object_properties['has_affected_locus'])
                        else:
                            logger.warn("unable to map intron variant"
                                        " to gene coordinates: {0}"
                                        .format(variant))
                            for r_gene in ref_gene:
                                self._add_gene_to_graph(
                                    r_gene['symbol'], variant_bnode, gene_id,
                                    model.object_properties['causally_influences'])
                    elif re.search(r'intron', variant['type'], flags=re.I):
                        logger.warn("unable to map intron variant"
                                    " to gene coordinates: {0}"
                                    .format(variant))
                    for neighbor in up_down_gene:
                        self._add_gene_to_graph(
                            neighbor, variant_bnode, gene_id,
                            model.object_properties['causally_influences'])
                    # Unmatched genes are likely because we cannot map to an NCBIGene
                    # or we do not have coordinate information
                    for unmatched_gene in unmatched_genes:
                        self._add_gene_to_graph(
                            unmatched_gene, variant_bnode, gene_id,
                            model.object_properties['causally_influences'])

        return

    @staticmethod
    def _convert_variant_file_to_dict(file):
        """
        Converts tsv to dicts with this structure
        {
            'patient_1': {
                'variant-id': {
                    'build': 'hg19'
                    'chromosome': 'chr7',
                    'reference_allele': 'A',
                    'variant_allele': 'G',
                    'position': '1234',
                    'rs_id' : 'RS1234',
                    'type': 'SNV',
                    'genes_of_interest' : [SHH, BRCA1]
                }
            }
        }
        TODO the above structure can be simplified as it results
             in duplicated variant dicts, the join-by-id approach
             below would be an improvement:
        {
            'patient_1': {
                'variants': [variant-id1,variant-id2, 3, 4]
            }
            'variants : [
               'variant-id1': {
                    'build': hg19
                    'chromosome': 'chr7',
                    'reference_allele': 'A',
                    'variant_allele': 'G',
                    'position': '1234'
                    'rs_id' : 'RS1234',
                    'type': 'SNV",
                    'genes_of_interest' : [SHH, BRCA1]
                }
            ]
        }

        If any part of the core variant information is missing
        (build, chr, bp change(s), the line number will be used
        to make the variant unique

        Variant id will be used downstream to form blank nodes (checksumed)

        See docstring for _add_variant_gene_relationship for explanation
        on why there is a one to many mapping between variants and genes

        Values are normalized with these rules:
        1. Basepairs are upper case
        2. HG19 -> hg19
        3. X -> chrX
        :return: dict
        """
        patient_variant_map = {}
        line_num = 0
        reader = csv.reader(file, delimiter="\t")
        for row in reader:

            (patient, family, chromosome, build, position,
             reference_allele, variant_allele, parent_of_origin,
             allele_type, mutation_type, gene_symbol, transcript,
             reference_aa, variant_aa, aa_change, segregates_with,
             locus, exon, inheritance_model, zygosity, dbSNP_ID, frequency,
             num_of_alleles) = row

            if patient == 'Patient':  # skip header
                line_num += 1
                continue

            if patient not in patient_variant_map:
                patient_variant_map[patient] = {}

            formatted_chr = re.sub(r'^CHR', 'chr', chromosome, flags=re.I)

            if re.fullmatch(r'[XY]|[0-9]{1,2}', chromosome, flags=re.I):
                formatted_chr = "chr{0}".format(chromosome.upper())

            formatted_build = re.sub(r'^HG', 'hg', build, flags=re.I)
            ref_base = reference_allele.upper()
            var_base = variant_allele.upper()
            rs_id = ''

            # Catch misformatted data
            if re.search(r'LEFT FLANK|NM_|EXON', ref_base):
                ref_base = ''

            if re.search(r'LEFT FLANK|NM_|EXON', var_base):
                var_base = ''

            if dbSNP_ID != '':
                match = re.fullmatch(r'^(rs\d+).*', dbSNP_ID)
                if match:
                    rs_id = match.group(1)

            # Format variant object
            variant_info = [formatted_chr, formatted_build, position,
                            ref_base, var_base]

            if '' in variant_info:
                filt_list = [info for info in variant_info if info != '']
                variant_id = str(line_num) + '-' + '-'.join(filt_list)
            else:
                variant_id = '-'.join(variant_info)

            if variant_id in patient_variant_map[patient]:
                patient_variant_map[patient][variant_id]['genes_of_interest'].append(gene_symbol)
            else:
                patient_variant_map[patient][variant_id] = {
                    'build': formatted_build,
                    'position': position,
                    'chromosome': formatted_chr,
                    'reference_allele': ref_base,
                    'variant_allele': var_base,
                    'type': mutation_type,
                    'rs_id': ''
                }
                if rs_id:
                    patient_variant_map[patient][variant_id]['rs_id'] = rs_id

                patient_variant_map[patient][variant_id]['genes_of_interest'] = [gene_symbol]

            line_num += 1

        return patient_variant_map

    def _parse_patient_phenotypes(self, file, limit=None):
        """
        :param file: file handler
        :param limit: limit rows processed
        :return:
        """
        model = Model(self.graph)
        line_counter = 0
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            (patient_id, hpo_curie, present) = row
            patient_curie = ':{0}'.format(patient_id)
            if patient_id == 'Patient':  # skip header
                line_counter += 1
                continue

            model.addPerson(patient_curie, patient_id)

            self.graph.addTriple(patient_curie,
                                 model.object_properties['has_phenotype'],
                                 "DOID:4")
            if present == 'yes':
                self.graph.addTriple(patient_curie,
                                     model.object_properties['has_phenotype'],
                                     hpo_curie)

            line_counter += 1
            if not self.testMode and limit is not None \
                    and line_counter >= limit:
                break

    @staticmethod
    def _parse_gene_coordinates(file):
        """
        :param file: file path
        :param limit: limit (int, optional) limit the number of rows processed
        :return: dict
        """
        id_map = {}
        if os.path.exists(os.path.join(os.path.dirname(__file__), file)):
            with open(os.path.join(os.path.dirname(__file__), file)) as tsvfile:
                reader = csv.reader(tsvfile, delimiter="\t")
                for row in reader:
                    (gene_curie, start, end, strand, build) = row
                    id_map[gene_curie] = {
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'build': build
                    }
        return id_map

    @staticmethod
    def _parse_rs_map_file(file):
        """
        Parses rsID mapping file from dbSNP
        Outputs dict where keys are coordinates in the format
        {chromsome}-{position}

        {
            chr1-1234: [
                {
                    'type': 'snp'
                    'rs_id': 'rs1234'
                    'alleles': 'A/G/T'
                }
            ]
        }

        :param file: file path
        :param limit: limit (int, optional) limit the number of rows processed
        :return: dict
        """
        rs_map = {}
        if os.path.exists(os.path.join(os.path.dirname(__file__), file)):
            with open(os.path.join(os.path.dirname(__file__), file)) as tsvfile:
                reader = csv.reader(tsvfile, delimiter="\t")
                for row in reader:
                    (chromosome, position, rs_id, var_type, alleles) = row
                    map_key = "chr{0}-{1}".format(chromosome, position)
                    rs_info = {
                        'type': var_type,
                        'rs_id': rs_id,
                        'alleles': alleles
                    }
                    if map_key in rs_map:
                        rs_map[map_key].append(rs_info)
                    else:
                        rs_map[map_key] = [rs_info]
        return rs_map

    @staticmethod
    def _build_variant_label(build, chromosome, position,
                             reference_allele, variant_allele,
                             gene_symbols=None):
        """
        Function to build HGVS variant labels
        :param build: {str} build id
        :param chromosome: {str} chromosome
        :param position: {str} variation position as string or int
        :param reference_allele: {str} single letter ref bp
        :param variant_allele: {str} single letter bp change
        :param gene_symbol: {str} gene symbol (hgvs)
        :return: {str} variant label
        """
        variant_label = ''
        prefix = ''
        if gene_symbols and len(gene_symbols) == 1 and gene_symbols[0]:
            prefix = "{0}{1}({2})".format(build, chromosome, gene_symbols[0])
        else:
            prefix = "{0}{1}".format(build, chromosome)
        if reference_allele == '-':
            variant_label = "{0}:g.{1}ins{2}".format(
                            prefix, position, variant_allele)
        elif variant_allele == '-':
            variant_label = "{0}:g.{1}del{2}".format(
                            prefix, position, reference_allele)
        else:
            variant_label = "{0}:g.{1}{2}>{3}".format(
                prefix, position, reference_allele, variant_allele)
        return variant_label

    def _add_gene_to_graph(self, gene, variant_bnode, gene_id, relation):
        """
        :param gene:
        :param variant_bnode:
        :return:
        """
        model = Model(self.graph)
        if gene_id:
            self.graph.addTriple(variant_bnode, relation, gene_id)
        elif gene:
            logger.info("gene {0} not mapped to NCBI gene,"
                        " making blank node".format(gene))
            gene_bnode = self.make_id("{0}".format(gene), "_")
            model.addIndividualToGraph(gene_bnode, gene)
            self.graph.addTriple(variant_bnode, relation, gene_bnode)

    def _add_variant_sameas_relationships(self, patient_var_map, rs_map):
        """
        Adds same as relationships between udp variant bnodes and dbsnp ids
        :param patient_var_map:
        :param rs_map:
        :return:
        """
        model = Model(self.graph)
        for patient in patient_var_map:
            for variant_id, variant in patient_var_map[patient].items():
                variant_bnode = self.make_id("{0}".format(variant_id), "_")
                build = variant['build']
                chromosome = variant['chromosome']
                position = variant['position']
                reference_allele = variant['reference_allele']
                variant_allele = variant['variant_allele']
                if build and chromosome and position\
                        and reference_allele and variant_allele:
                    if re.fullmatch(r'[ATCG]', reference_allele)\
                            and re.fullmatch(r'[ATCG]', variant_allele):
                        # variation is snp
                        rs_id = self._get_rs_id(variant, rs_map, 'snp')
                        if rs_id:
                            dbsnp_curie = 'dbSNP:rs{0}'.format(rs_id)
                            model.addSameIndividual(variant_bnode, dbsnp_curie)

                    elif re.fullmatch(r'\-', reference_allele)\
                            or re.fullmatch(r'\-', variant_allele):
                        rs_id = self._get_rs_id(variant, rs_map, 'indel')
                        if rs_id is not None:
                            dbsnp_curie = 'dbSNP:rs{0}'.format(rs_id)
                            model.addSameIndividual(variant_bnode, dbsnp_curie)
                    else:
                        rs_id = self.\
                            _get_rs_id(variant, rs_map, 'indel')
                        if rs_id is not None:
                            dbsnp_curie = 'dbSNP:rs{0}'.format(rs_id)
                            model.addSameIndividual(variant_bnode, dbsnp_curie)
        return

    @staticmethod
    def _get_rs_id(variant, rs_map, variant_type):
        """
        Given a variant dict, return unambiguous RS ID
        TODO
        Some sequence alterations appear to have mappings to dbsnp's notation
        for example,
        reference allele: TTTTTTTTTTTTTT
        variant allele:   TTTTTTTTTTTTTTT
        Is theoretically the same as -/T, we should clarify with UDP and then add
        functionality to map this notation to the more common -/T
        :param variant:
        :param rs_map:
        :param type: snp or indel
        :return:
        """
        rs_id = None
        if variant_type == 'snp':
            variant_key = "{0}-{1}".format(variant['chromosome'], variant['position'])
            if variant_key in rs_map:
                snp_candidates = [rs_dict for rs_dict in rs_map[variant_key] if rs_dict['type'] == 'snp']
                if len(snp_candidates) == 1:
                    rs_id = snp_candidates[0]["rs_id"]
        elif variant_type == 'indel':
            rs_candidates = []
            variant_key = "{0}-{1}".format(variant['chromosome'], variant['position'])

            if variant_key in rs_map:
                snp_candidates = [rs_dict for rs_dict in rs_map[variant_key] if rs_dict['type'] == 'in-del']
                for candidate in snp_candidates:
                    alleles = candidate['alleles'].split('/')
                    if variant['reference_allele'] in alleles \
                            and variant['variant_allele'] in alleles:
                        rs_candidates.append(candidate['rs_id'])

                if len(rs_candidates) == 1:
                    rs_id = rs_candidates[0]
                elif len(rs_candidates) > 1:
                    logger.info("ambiguous rs mapping for:"
                                " {0}\n candidate ids: {1}".format(variant, rs_candidates))
                else:
                    logger.info("rs at coordinate but no match found"
                                " for variant {0}\n candidate ids: {1}".format
                                (variant, rs_map[variant_key]))
        else:
            logger.warn("type: {0} unsupported".format(variant_type))
        return rs_id
