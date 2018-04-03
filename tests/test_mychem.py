import json
import unittest
from dipper.sources.MyChem import MyChem
from dipper.graph.RDFGraph import RDFGraph

RESOURCE = ""

class TestGwasSNPModel(unittest.TestCase):
    """
    Test the modelling of a  SNP to trait association
    from sample GWAS catalog data
    """

    def setUp(self):
        self.source = MyChem('rdf_graph', True)
        self.source.graph = RDFGraph(True) # Reset graph
        self.test_data = ""

    def tearDown(self):
        self.source = None
