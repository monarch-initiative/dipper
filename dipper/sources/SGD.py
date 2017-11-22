from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.models.Provenance import Provenance
from dipper.models.Dataset import Dataset
from dipper.models.Reference import Reference
from ontobio.ontol_factory import OntologyFactory
import logging
from pprint import pprint
import pandas as pd

__author__ = 'timputman'

logger = logging.getLogger(__name__)


class SGD(Source):
    """
    Ingest of Saccharomyces Genome Database (SGD) phenotype associations

    """
    SGD_BASE = 'https://downloads.yeastgenome.org/curation/literature/'
    files = {
        'sgd_phenotype': {
            'file': 'phenotype_data.tab',
            'url': SGD_BASE + 'phenotype_data.tab'},
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'sgd')
        self.dataset = Dataset(
            'sgd', 'SGD', 'https://www.yeastgenome.org/', None,
            None)

        self.global_terms = Source.open_and_parse_yaml('../../translationtable/global_terms.yaml')

    def fetch(self, is_dl_forced=False):
        """
        Override Source.fetch()
        Fetches resources from rat_genome_database using the rat_genome_database ftp site
        Args:
            :param is_dl_forced (bool): Force download
        Returns:
            :return None
        """
        self.get_files(is_dl_forced)
        return

    def parse(self, limit=None):
        """
        Override Source.parse()
        Args:
            :param limit (int, optional) limit the number of rows processed
        Returns:
            :return None
        """
        if limit is not None:
            logger.info("Only parsing first %d rows", limit)

        sgd_file = '/'.join((self.rawdir, self.files['sgd_phenotype']['file']))
        columns = ['Feature Name', 'Feature Type', 'Gene Name', 'SGDID', 'Reference', 'Experiment Type', 'Mutant Type',
                   'Allele', 'Strain Background', 'Phenotype', 'Chemical', 'Condition', 'Details', 'Reporter']
        sgd_df = pd.read_csv(sgd_file, sep='\t', names=columns)
        records = sgd_df.to_dict(orient='records')

        # load apo for term mapping
        ofactory = OntologyFactory()
        apo_ont = ofactory.create("apo")
        apo_nodes = apo_ont.nodes()
        # dict schema { 'term': 'apo_id' }
        apo_term_id = dict()
        for node in apo_nodes:
            apo_term_id[apo_ont.label(node)] = node

        for index, assoc in enumerate(records):
            if limit is not None and index > limit:
                break
            if index < 10:
                # removed description and mapp Experiment Type to apo term
                experiment_type = assoc['Experiment Type'].split('(')[0]
                assoc['experiment_type'] = apo_term_id[experiment_type]
                sgd_phenotype = assoc['Phenotype']
                pheno_obj = {
                    'entity': {
                        'term': None,
                        'apo_id': None
                    },
                    'quality': {
                        'term': None,
                        'apo_id': None
                    },
                    'has_quality': False  # False = phenotype was descriptive and don't bother looking for a quality
                }
                phenotype = assoc['Phenotype']
                if ':' in phenotype:
                    pheno_obj['has_quality'] = True
                    ent_qual = sgd_phenotype.split(': ')
                    entity = ent_qual[0]
                    quality = ent_qual[1]
                    pheno_obj['entity']['term'] = entity
                    pheno_obj['entity']['apo_id'] = apo_term_id[entity]
                    pheno_obj['quality']['term'] = quality
                    pheno_obj['quality']['apo_id'] = apo_term_id[quality]
                else:
                    pheno_obj['entity']['term'] = phenotype
                    pheno_obj['entity']['apo_id'] = apo_term_id[phenotype]
                assoc['pheno_obj'] = pheno_obj
                self.make_association(assoc)

        return

    def make_association(self, record):
        """
        contstruct the association
        :param record:
        :return: modeled association of  genotype to mammalian phenotype
        """
        model = Model(self.graph)

        # define the triple
        gene = 'SGD:{}'.format(record['SGDID'])
        relation = Model.object_properties['has_phenotype']  # has phenotype
        phenotype = record['pheno_obj']['entity']['apo_id']
        g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=phenotype, pred=relation)
        g2p_assoc.make_association_id(definedby='http://identifiers.org/SGD',
                                      subject=gene,
                                      predicate=relation,
                                      object=phenotype)
        g2p_assoc.set_association_id()
        if record['pheno_obj']['has_quality']:
            # add quality to association
            quality = record['pheno_obj']['quality']['apo_id']
            has_quality = Model.object_properties['has_quality']
            g2p_assoc.add_predicate_object(predicate=has_quality, object_node=quality) # this is where it is breaking

        # # # add the references
        references = record['Reference']
        references = references.replace(' ', '')
        references = references.split('|')
        #
        # # # created RGDRef prefix in curie map to route to proper reference URL in RGD
        # references = [x.replace('RGD', 'RGDRef') if 'PMID' not in x else x for x in references]
        if len(references) > 0:
            # make first ref in list the source
            g2p_assoc.add_source(identifier=references[0])
            ref_model = Reference(
                self.graph, references[0],
                Reference.ref_types['publication']
            )
            ref_model.addRefToGraph()
        #
        # # if len(references) > 1:
            # create equivalent source for any other refs in list
            for ref in references[1:]:
                ref = ref.replace('SGD_REF', 'SGD')
                model.addSameIndividual(sub=references[0], obj=ref)
        #
        g2p_assoc.add_evidence(record['experiment_type'])

        try:
            g2p_assoc.add_association_to_graph()
        except Exception as e:
            print(e)
        return


example_record = {'Allele': ' ',
                  'Chemical': ' ',
                  'Condition': ' ',
                  'Details': ' ',
                  'Experiment Type': 'classical genetics',
                  'Feature Name': 'IMI1',
                  'Feature Type': 'not in systematic sequence of S288C',
                  'Gene Name': 'IMI1',
                  'Mutant Type': None,
                  'Phenotype': 'mitochondrial genome maintenance: abnormal',
                  'Reference': 'PMID: 26091838|SGD_REF: S000180603',
                  'Reporter': ' ',
                  'SGDID': 'S000149345',
                  'Strain Background': 'W303',
                  'pheno_obj': {'entity': {'apo_id': 'APO:0000105',
                                           'term': 'mitochondrial genome maintenance'},
                                'has_quality': True,
                                'quality': {'apo_id': 'APO:0000002', 'term': 'abnormal'}}}
