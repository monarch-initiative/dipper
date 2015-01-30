import csv
import re
import gzip
from sources.Source import Source
from models.Dataset import Dataset


class CTD(Source):

    files = {
        'interactions': {'file': 'CTD_chemicals_diseases.tsv.gz',
                         'url': 'http://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz'}
    }
    fetchdate = ''

    def __init__(self):
        Source.__init__(self, 'ctd')
        self.load_bindings()
        self.dataset = Dataset('ctd', 'CTD', 'http://ctdbase.org')

    def fetch(self):
        """
        :return: None
        """
        fetchdate = self.fetch_files(self.files)
        self.fetchdate = fetchdate
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
        self._get_interactions(limit, self.files['interactions']['file'])

        print("Done parsing files.")

        return

    def _get_interactions(self, limit, file):
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
                        self._set_version(match.group(1))
                        is_versioned = True
                elif re.match('^#', ' '.join(row)):
                    next
                else:
                    # only get direct associations
                    if row[5] != '':
                        # Process data here
                        foo=1
                    row_count += 1
                    if limit is not None and row_count >= limit:
                        break

    def _set_version(self, version):
        """
        Replace whitespace and colons with dashes
        and set the dataset version
        :param version: string containing version ID
        :return:None
        """
        version = re.sub(r'\s|:', '-', version)
        self.dataset.setVersion(self.fetchdate, version)
        return

    def _add_evidence_associations(self, subject, object):
        """
        :param subject
        :param object
        :return:None
        """

        return
