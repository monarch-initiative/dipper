from dipper.sources.Source import Source
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
import logging
from SPARQLWrapper import SPARQLWrapper, JSON
import requests

__author__ = 'timputman'

logger = logging.getLogger(__name__)


class MyChem(Source):
    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'mychem')
        self.dataset = Dataset(
            'mychem', 'MYCHEM', 'https://mychem.info/', None,
            None)

        self.global_terms = Source.open_and_parse_yaml('../../translationtable/global_terms.yaml')
        self.inchikeys = MyChem.chunks(l=MyChem.get_inchikeys(), n=10)
        self.drugbank_targets = list()
        self.drugcentral_interactors = list()

    def fetch(self, is_dl_forced=False):
        self.fetch_from_mychem()
        return

    def parse(self, limit=None):
        try:
            for index, record in enumerate(self.drugbank_targets):
                record_data = {
                    'drugbank_id': 'DrugBank:{}'.format(record['drugbank']['drugbank_id']),
                    'unii': None,
                    'targets': [],

                }
                if 'unii' in record.keys():
                    if isinstance(record['unii'], dict):
                        record_data['unii'] = 'UNII:{}'.format(record['unii']['unii'])
                else:
                    continue
                if 'targets' in record['drugbank'].keys():
                    targets = record['drugbank']['targets']
                    targets = MyChem.return_target_list(targets)
                    for targ in targets:
                        uniprot = MyChem.check_uniprot(targ)
                        if uniprot is not None:
                            actions = MyChem.format_actions(targ)
                            for act in actions:
                                record_data['targets'].append({
                                    'uniprot': uniprot,
                                    'action': 'MONARCH:{}'.format(act),
                                    'name': targ['name']
                                })
                else:
                    continue
                self.make_triples(source='drugbank', package=record_data)
            for index, record in enumerate(self.drugcentral_interactors):
                record_data = {
                    'interactions': [],
                    'indications': [],
                }
                if 'unii' in record.keys():
                    if isinstance(record['unii'], dict):
                        record_data['unii'] = 'UNII:{}'.format(record['unii']['unii'])
                else:
                    continue
                if 'bioactivity' in record['drugcentral'].keys():
                    bioactivity = record['drugcentral']['bioactivity']
                    bioactivity = MyChem.return_target_list(bioactivity)
                    for bioact in bioactivity:
                        target_class = None
                        target_name = None
                        if 'target_class' in bioact.keys():
                            target_class = bioact['target_class']
                        if 'target' in bioact.keys():
                            target_name = bioact['target']
                        if 'uniprot_id' in bioact.keys():
                            uniprot_ids = bioact['uniprot_id'].split("|")
                            for up in uniprot_ids:
                                record_data['interactions'].append({
                                    'uniprot': 'UniProtKB:{}'.format(up),
                                    'target_class': target_class,
                                    'target_name': target_name,
                                })
                        else:
                            continue
                elif 'drug_use' in record['drugcentral'].keys():
                    drug_use = record['drugcentral']['drug_use']
                    drug_use = MyChem.return_target_list(drug_use)
                    for dr_record in drug_use:
                        if dr_record['relation'] != 'contraindication' and 'snomed_id' in dr_record.keys():
                            record_data['indications'].append({
                                'snomed_id': 'SNOMED:{}'.format(dr_record['snomed_id']),
                                'snomed_name': dr_record['snomed_name'],
                            })

                else:
                    continue
                self.make_triples(source='drugcentral', package=record_data)
        except Exception as e:
            print(e)
        return

    def make_triples(self, source, package):
        model = Model(self.graph)
        if source == 'drugbank':
            for target in package['targets']:
                model.addTriple(subject_id=package['unii'],predicate_id=target['action'],obj=target['uniprot'])
                model.addLabel(subject_id=target['uniprot'], label=target['name'])
                model.addTriple(subject_id=target['uniprot'],
                                predicate_id=Model.object_properties['subclass_of'],
                                obj='SO:0000104')
                model.addTriple(subject_id=package['drugbank_id'],
                                predicate_id=Model.object_properties['equivalent_class'],
                                obj=package['unii'])
                model.addTriple(subject_id=target['action'],
                                predicate_id='rdfs:subPropertyOf',
                                obj='RO:0002436')
                model.addTriple(subject_id=package['unii'],
                                predicate_id=Model.object_properties['subclass_of'],
                                obj='CHEBI:23367')
        if source == 'drugcentral':
            for indication in package['indications']:
                model.addTriple(subject_id=package['unii'], predicate_id='RO:0002606', obj=indication['snomed_id'])
                model.addTriple(subject_id=package['unii'],
                                predicate_id=Model.object_properties['subclass_of'],
                                obj='CHEBI:23367')
                model.addTriple(subject_id=indication['snomed_id'],
                                predicate_id=Model.object_properties['subclass_of'],
                                obj='DOID:4')
                model.addLabel(subject_id=indication['snomed_id'], label=indication['snomed_name'])
            for interaction in package['interactions']:
                model.addTriple(subject_id=package['unii'], predicate_id='RO:0002436', obj=interaction['uniprot'])
                # model.addLabel(subject_id=interaction['uniprot'], label='Protein_{}'.format(interaction['uniprot']))
                model.addLabel(subject_id=interaction['uniprot'], label=interaction['target_name'])
                model.addTriple(subject_id=package['unii'],
                                predicate_id=Model.object_properties['subclass_of'],
                                obj='CHEBI:23367')
                model.addDescription(subject_id=interaction['uniprot'], description=interaction['target_class'])
                model.addTriple(subject_id=interaction['uniprot'],
                                predicate_id=Model.object_properties['subclass_of'],
                                obj='SO:0000104')


        return

    def fetch_from_mychem(self):
        count = 0
        for k in self.inchikeys:
            count +=1
            if count < 10:
                ids = ",".join(k)
                fields = 'drugbank.targets,drugbank.drugbank_id,unii.unii,drugcentral.drug_use,drugcentral.bioactivity'
                records = MyChem.get_drug_record(ids=ids, fields=fields)
                for record in records:
                    if 'drugbank' in record.keys():
                        self.drugbank_targets.append(record)
                    if 'drugcentral' in record.keys():
                        self.drugcentral_interactors.append(record)
            return

    @staticmethod
    def add_relation(results, relation):
        triples = list()
        for index, result in enumerate(results['results']['bindings']):
            result.update({'relation': relation})
            triples.append(result)
        return triples

    @staticmethod
    def execute_query(query):
        endpoint = SPARQLWrapper('https://query.wikidata.org/sparql')
        endpoint.setQuery(query)
        endpoint.setReturnFormat(JSON)
        return endpoint.query().convert()

    @staticmethod
    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    @staticmethod
    def get_inchikeys():
        all_ids = list()
        query = '''
        SELECT ?drug ?inchi WHERE {
        ?drug wdt:P235 ?inchi;
              wdt:P652 ?unii.
        }
        '''
        r = MyChem.execute_query(query=query)
        for identifier in r['results']['bindings']:
            all_ids.append(identifier['inchi']['value'])
        return all_ids

    @staticmethod
    def get_drug_record(ids, fields):
        url = 'http://mychem.info/v1/drug'
        params = {
            'ids': ids,
            'fields': fields
        }

        r = requests.post(url=url, params=params)
        return r.json()

    @staticmethod
    def format_actions(target_dict):
        final = list()

        def replace_characters(action):
            if ' ' in action:
                action = action.split()
                return "_".join(action)
            elif '/' in action:
                action = action.split('/')
                return "_".join(action)
            else:
                return action

        if 'actions' in target_dict.keys():
            if isinstance(target_dict['actions'], str):
                final.append(replace_characters(target_dict['actions']))
            if isinstance(target_dict['actions'], list):
                for act in target_dict['actions']:
                    final.append(replace_characters(act))
        return final

    @staticmethod
    def check_uniprot(target_dict):
        if 'uniprot' in target_dict.keys():
            return 'UniProtKB:{}'.format(target_dict['uniprot'])

    @staticmethod
    def return_target_list(targ_in):
        targ_list = list()
        if isinstance(targ_in, dict):
            targ_list.append(targ_in)
        if isinstance(targ_in, list):
            targ_list = targ_list + targ_in
        return targ_list

