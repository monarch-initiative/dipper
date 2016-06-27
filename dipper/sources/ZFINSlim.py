from intermine.webservice import Service
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.sources.Source import Source
from dipper.sources.ZFIN import ZFIN
from dipper.models.Dataset import Dataset

import logging
import datetime

logger = logging.getLogger(__name__)


class ZFINSlim(Source):
    """
    zfin mgi model only containing Gene to phenotype associations
    Using the file here: https://zfin.org/downloads/phenoGeneCleanData_fish.txt
    """
    files = {
        'g2p_clean': {
            'file': 'phenoGeneCleanData_fish.txt.txt',
            'url': 'https://zfin.org/downloads/phenoGeneCleanData_fish.txt'
        },
        'zpmap': {
            'file': 'zp-mapping.txt',
            'url': 'http://compbio.charite.de/hudson/job/zp-owl-new/lastSuccessfulBuild/artifact/zp.annot_sourceinfo'
        }
    }

    def __init__(self):
        super().__init__('zfin-slim')
        self.dataset = Dataset(
            'zfin_slim', 'ZFINSlim', 'http://zfin.org/')

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        self.load_bindings()
        zfin_parser = ZFIN()
        zp_file = '/'.join((self.rawdir, self.files['zpmap']['file']))
        zp_map = zfin_parser._load_zp_mappings(zp_file)

