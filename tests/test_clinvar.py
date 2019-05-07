#!/usr/bin/env python3
import unittest
import logging
import os
from unittest.mock import patch
from dipper.graph.RDFGraph import RDFGraph
from dipper.utils.rdf2dot import rdf2dot
from dipper.utils.TestUtils import TestUtils
from dipper.sources.ClinVarXML_alpha import parse as clinvar_parse

logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)
LOG = logging.getLogger(__name__)

TEST_PATH = os.path.join(os.path.dirname(__file__), 'resources/clinvar')
NT_PATH  = TEST_PATH + "/nt/"
DOT_PATH = TEST_PATH + "/dot/"
XML_PATH = TEST_PATH + "/input/"
TTL_PATH = TEST_PATH + "/output/"

# IDs for the files in resources/clinvar/{RCV}.xml.gz, reused for output {RCV}.nt
RCVS = [
    'RCV000112698',
    'RCV000162061',
    'RCV000175394',
    'RCV000416376',
    'RCV000498447',
    'RCV000763295',
    'RCV000087646'
]

MAP_FILE = 'gene_condition_source_id'


class ClinVarTestCase(unittest.TestCase):

    def test_parse(self):
        for rcv in RCVS:
            output_nt = rcv + '.nt'
            input_xml = rcv + '.xml.gz'
            reference_ttl = TTL_PATH + rcv + '.ttl'
            with self.subTest(rcv=rcv):

                mock_args = [
                    "test_clinvar.py",
                    "--inputdir", XML_PATH,
                    "--filename", input_xml,
                    "--mapfile",  MAP_FILE,
                    "--destination", NT_PATH,
                    "--output", output_nt
                ]

                patch('sys.argv', mock_args).start()
                clinvar_parse()
                query_graph = RDFGraph()
                query_graph.bind_all_namespaces()
                query_graph.parse(NT_PATH + output_nt, format='nt')

                with open(reference_ttl, 'r') as ref_fh:
                    ref_graph = "\n".join(ref_fh.readlines())

                # debug
                LOG.debug(
                    "Reference graph: %s",
                    query_graph.serialize(format="turtle").decode("utf-8"))

                # Convert output from ClinVar parse to dot then png
                dot_file_path = DOT_PATH + rcv + ".dot"
                with open(dot_file_path, 'w') as dot_file:
                    rdf2dot(query_graph, dot_file)

                self.assertTrue(TestUtils.test_graph_equality(ref_graph, query_graph))

                #self.assertTrue(True)
