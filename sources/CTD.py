import csv, os, datetime, re, gzip
from datetime import datetime
from sources.Source import Source
from models.Dataset import Dataset


class CTD(Source):

    files = {
        'interactions': {'file': 'CTD_chemicals_diseases.tsv.gz',
                         'url': 'http://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz'}
    }
    date = ''

    def __init__(self):
        Source.__init__(self, 'ctd')
        self.load_bindings()
        self.dataset = Dataset('foo', 'bar', 'baz')

    def fetch(self):
        """
        :return: None
        """
        filedate = self.fetch_files(self.files)
        self.date = filedate
        return

    def parse(self, limit=None):
        """
        :return:None
        """
        version_pattern = re.compile('^# Report created: (.+)$')
        is_versioned = False
        file_path = '/'.join((self.rawdir, self.files['interactions']['file']))
        with gzip.open(file_path, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                # Scan the header lines until we get the version
                # There is no official version sp we are using
                # the upload timestamp instead
                if is_versioned is False:
                    match = re.match(version_pattern,' '.join(row))
                    if match:
                       self.set_version(match.group(1))
                       is_versioned = True

        return

    def set_version(self,version):
        return






