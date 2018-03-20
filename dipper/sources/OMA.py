from .OrthoXML import OrthoXML
from dipper.models.Dataset import Dataset
import logging

__author__ = 'Adrian Altenhoff'

logger = logging.getLogger(__name__)


class OMA(OrthoXML):
    BENCHMARK_BASE = "https://omabrowser.org/ReferenceProteomes"
    files = {'oma_hogs': {
        'file': 'OMA_GETHOGs-2_2017-04.orthoxml.gz',
        'url': BENCHMARK_BASE + "/OMA_GETHOGs-2_2017-04.orthoxml.gz"},
    }

    def __init__(self, graph_type, are_bnodes_skolemized, tax_ids=None):
        super().__init__(graph_type, are_bnodes_skolemized, 'OMA', tax_ids=tax_ids)
        self.dataset = Dataset(
            'OMA_HOGs', 'Ortholgous MAtrix Hierarchical Orthologous Groups',
            'https://omabrowser.org/',
            license_url="https://creativecommons.org/licenses/by-sa/2.5/")

    def getTestSuite(self):
        import unittest
        from tests.test_oma import OMATestCase

        test_suite = unittest.TestLoader().loadTestsFromTestCase(
            OMATestCase)

        return test_suite
