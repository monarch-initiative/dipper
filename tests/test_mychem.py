import json
import unittest
import logging
from dipper.sources.MyChem import MyChem
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.TestUtils import TestUtils
import os

logging.basicConfig()
logging.getLogger().setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

TESTDATA = os.path.join(os.path.dirname(__file__), 'resources/mychem/mychem.json')


class TestMyChemParser(unittest.TestCase):

    def setUp(self):
        self.test_util = TestUtils()
        self.source = MyChem('rdf_graph', True)

        # Replaces source.fetch()
        data_fh = open(TESTDATA, 'r')
        self.test_data = json.load(data_fh)
        data_fh.close()
        self.source.drugbank_targets.append(self.test_data[0])
        self.source.drugcentral_interactors.append(self.test_data[0])

    def tearDown(self):
        self.source = None

    def test_parse(self):
        self.source.graph = RDFGraph(True)  # Reset graph
        self.assertTrue(len(list(self.source.graph)) == 0)

        self.source.parse()

        triples = """
        UNII:46U771ERWK RO:0002606 SNOMED:386761002 ;
            rdfs:subClassOf CHEBI:23367 .

        SNOMED:386761002 rdfs:label "Local anesthesia" ;
            rdfs:subClassOf DOID:4 .
        """

        # dbg
        logger.debug(
            "Reference graph: %s", self.source.graph.serialize(format="turtle")
            .decode("utf-8")
        )
        self.assertTrue(self.test_util.test_graph_equality(
            triples, self.source.graph))
