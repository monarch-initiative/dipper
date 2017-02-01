import logging
import requests
from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.assoc.Association import Assoc
from dipper.models.Evidence import Evidence
from dipper.models.Provenance import Provenance
from dipper.models.Model import Model
from dipper import curie_map
from pathlib import Path
import json

logger = logging.getLogger(__name__)


class MyDrug(Source):
    """
    Drugs and Compounds stored in the BioThings database
    """
    MY_DRUG_API = 'http://c.biothings.io/v1/query'

    files = {
        'aeolus': {
            'file': 'aeolus.json'
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'mydrug')

        self.dataset = Dataset(
            'MyDrug', 'Drugs and Compounds in BioThings',
            'http://c.biothings.io/')

    def fetch(self, is_dl_forced=False):
        """
        Note there is a unpublished mydrug client that works like this:
        from mydrug import MyDrugInfo
        md = MyDrugInfo()
        r = list(md.query('_exists_:aeolus', fetch_all=True))

        :param is_dl_forced: boolean, force download
        :return:
        """
        dir_path = Path(self.rawdir)
        aeolus_file = dir_path / self.files['aeolus']['file']
        if self.checkIfRemoteIsNewer(aeolus_file):
            aeolis_fh = aeolus_file.open('w')
            aeolis_fh.write("[\n")
            params = {
                'q': '_exists_:aeolus',
                'from': 0,
                'rows': 10
            }
            result_count = params['rows']

            while params['from'] < result_count:
                solr_request = requests.get(self.MY_DRUG_API, params=params)
                response = solr_request.json()
                for index, doc in enumerate(response['hits']):
                    if params['from'] == 0 and index == 0:
                        aeolis_fh.write("{}".format(json.dumps(doc)))
                    else:
                        aeolis_fh.write(",\n{}".format(json.dumps(doc)))
                if params['from'] % 500 == 0:
                    logger.info("Fetched {} documents".format(params['from']))
                result_count = response['total']
                params['from'] += params['rows']

            aeolis_fh.write("\n]")
            aeolis_fh.close()

    def parse(self, limit=None, or_limit=1):
        """
        Parse mydrug files
        :param limit: int limit json docs processed
        :param or_limit: int odds ratio limit
        :return: None
        """
        dir_path = Path(self.rawdir)
        aeolus_file = dir_path / self.files['aeolus']['file']
        aeolus_fh = aeolus_file.open('r')
        count = 0
        for line in aeolus_fh.readlines():
            if limit is not None and count >= limit:
                break
            line = line.rstrip("\n,")
            if line != '[' and line != ']':
                self._parse_aeolus_data(document=json.loads(line),
                                        or_limit=or_limit)
                count += 1
            if count % 500 == 0:
                    logger.info("Processed {} documents".format(count))

        aeolus_fh.close()
        return

    def _parse_aeolus_data(self, document, or_limit=None):
        model = Model(self.graph)

        rxcui_curie = "RXCUI:{}".format(document['aeolus']['rxcui'])
        uni_curie = "UNII:{}".format(document['aeolus']['unii'])
        model.addLabel(rxcui_curie, document['aeolus']['drug_name'])
        model.addLabel(uni_curie, document['aeolus']['drug_name'])

        model.addSameIndividual(rxcui_curie, uni_curie)
        self.graph.addTriple(rxcui_curie,
                             model.annotation_properties['inchi_key'],
                             document['unii']['inchikey'],
                             object_is_literal=True)

        if or_limit is not None:
            outcomes = (outcome for outcome in document['aeolus']['outcomes']
                        if 'ror' in outcome and outcome['ror'] >= or_limit)
        else:
            outcomes = (outcome for outcome in document['aeolus']['outcomes'])

        for outcome in outcomes:
            drug2outcome_assoc = Assoc(self.graph, self.name)

            meddra_curie = "MEDDRA:{}".format(outcome['code'])
            model.addLabel(meddra_curie, outcome['name'])

            drug2outcome_assoc.sub = rxcui_curie
            drug2outcome_assoc.obj = meddra_curie
            drug2outcome_assoc.rel = Assoc.object_properties['causes_or_contributes']
            drug2outcome_assoc.description = \
                "A proportional reporting ratio or odds " \
                "ratio greater than or equal to {} in the " \
                "AEOLUS data was the significance cut-off " \
                "used for creating drug-outcome associations".format(or_limit)
            drug2outcome_assoc.add_association_to_graph()
            drug2outcome_assoc.add_predicate_object(
                Assoc.annotation_properties['probabalistic_quantifier'],
                outcome['ror'], 'Literal')

            self._add_outcome_evidence(drug2outcome_assoc.assoc_id, outcome)
            self._add_outcome_provenance(drug2outcome_assoc.assoc_id, outcome)

    def _add_outcome_provenance(self, association, outcome):
        """
        :param association: str association curie
        :param outcome: dict (json)
        :return: None
        """
        provenance = Provenance(self.graph)
        base = curie_map.get_base()

        provenance.add_agent_to_graph(base, 'Monarch Initiative')
        self.graph.addTriple(
            association, provenance.object_properties['asserted_by'], base)

    def _add_outcome_evidence(self, association, outcome):
        """
        :param association: str association curie
        :param outcome: dict (json)
        :return: None
        """
        evidence = Evidence(self.graph, association)
        source = {
            'curie': "DOI:10.5061/dryad.8q0s4/1",
            'label': "Data from: A curated and standardized adverse "
                     "drug event resource to accelerate drug safety research",
            'type': 'IAO:0000100'
        }
        reference = {
            'curie': "PMID:27193236",
            'label': None,
            'type': "IAO:0000311"
        }
        evidence_curie = self.make_id("{0}{1}{2}".format(
            association, outcome['id'], self.name
        ))
        evidence_type = "ECO:0000180"
        evidence.add_supporting_evidence(evidence_curie, evidence_type)

        evidence.add_supporting_publication(
            evidence_curie, reference['curie'], reference['label'],
            reference['type'])

        evidence.add_source(
            evidence_curie, source['curie'], source['label'], source['type'])

        count_bnode = self.make_id(
            "{0}{1}{2}".format(evidence_curie,
                               outcome['case_count'], self.name), prefix="_")
        pr_ratio_bnode = self.make_id(
            "{0}{1}{2}{3}".format(evidence_curie, outcome['prr'], self.name, 'prr'),
            prefix="_")
        odds_ratio_bnode = self.make_id(
            "{0}{1}{2}{3}".format(evidence_curie, outcome['ror'], self.name, 'ror'),
            prefix="_")

        evidence.add_data_individual(
            count_bnode, ind_type=evidence.data_types['count'])
        evidence.add_data_individual(
            pr_ratio_bnode,
            ind_type=evidence.data_types['proportional_reporting_ratio'])
        evidence.add_data_individual(
            odds_ratio_bnode, ind_type=evidence.data_types['odds_ratio'])

        value_map = {
            count_bnode: outcome['case_count'],
            pr_ratio_bnode: outcome['prr'],
            odds_ratio_bnode: outcome['ror']
        }
        evidence.add_supporting_data(evidence_curie, value_map)
        return

    # Override
    def checkIfRemoteIsNewer(self, localfile):
        """
        Need to figure out how biothings records releases,
        for now if the file exists we will assume it is
        a fully downloaded cache
        :param localfile: str file path
        :return: boolean True if remote file is newer else False
        """
        is_remote_newer = False
        if localfile.exists() \
                and localfile.stat().st_size > 0:
            logger.info("File exists locally, using cache")
        else:
            is_remote_newer = True
            logger.info("No cache file, fetching entries")
        return is_remote_newer
