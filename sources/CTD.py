import csv
import re
import gzip
from sources.Source import Source
from models.Dataset import Dataset
from models.Chem2DiseaseAssoc import Chem2DiseaseAssoc


class CTD(Source):

    files = {
        'interactions': {'file': 'CTD_chemicals_diseases.tsv.gz',
                         'url': 'http://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz'}
    }

    def __init__(self):
        Source.__init__(self, 'ctd')
        self.load_bindings()
        self.dataset = Dataset('ctd', 'CTD', 'http://ctdbase.org')

    def fetch(self, is_dl_forced):
        """
        :return: None
        """
        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        """
        Parses version and interaction information from CTD
        :param limit limit the number of rows processed
        :return:None
        """
        if limit is not None:
            print("Only parsing first", limit, "rows")

        print("Parsing files...")
        self._parse_interactions_file(limit, self.files['interactions']['file'])

        print("Done parsing files.")

        return

    def _parse_interactions_file(self, limit, file):
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
                    # only get direct associations
                    if row[5] != '':
                        # Process data here
                        self._process_interactions(row)
                    row_count += 1
                    if limit is not None and row_count >= limit:
                        break

    def _process_interactions(self, row):
        """
        :param row
        :return:None
        """
        self._check_list_len(row, 10)
        (chem_name, chem_id, cas_rn, disease_name, disease_id, direct_evidence,
         inferred_gene_symbol, inference_score, omim_ids, pubmed_ids) = row
        evidence_pattern = re.compile('^therapeutic|marker\/mechanism$')
        dual_evidence = re.compile('^marker\/mechanism\|therapeutic$')

        if re.match(evidence_pattern, direct_evidence):
            # TODO check if id/node already exists
            assoc_id = self.make_id('ctd' + chem_id + disease_id + direct_evidence)
            reference_list = self._process_pubmed_ids(pubmed_ids)
            evidence_code = self._set_evidence_code(direct_evidence)
            chem_mesh_id = 'MESH:'+chem_id

            assoc = Chem2DiseaseAssoc(assoc_id, chem_mesh_id, disease_id,
                                      reference_list, evidence_code)
            assoc.loadObjectProperties(self.graph)
            assoc.addAssociationNodeToGraph(self.graph)

        return

    def _process_pubmed_ids(self, pubmed_ids):
        """
        :param pubmed_ids -  string representing publication
                           ids seperated by a | symbol
        :return: list of ids with the PUBMED prefix
        """
        id_list = pubmed_ids.split('|')
        for (i, val) in enumerate(id_list):
            id_list[i] = 'PMID:'+val
        return id_list

    def _set_evidence_code(self, evidence):
        """
        :param evidence
        :return: ECO evidence code
        """
        ECO_MAP = {
            'therapeutic': 'ECO:0000269',
            'marker/mechanism': 'ECO:0000306'
        }
        return ECO_MAP[evidence]