from dipper.sources.Source import Source
from dipper.models.assoc.Association import Assoc
from dipper.models.Model import Model
from dipper.models.Dataset import Dataset
from dipper.models.Reference import Reference
from ontobio.ontol_factory import OntologyFactory
import logging
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
        self.apo_term_id = SGD.make_apo_map()

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
        for index, assoc in enumerate(records):
            if isinstance(assoc['Gene Name'], str):
                if limit is not None and index > limit:
                    break
                self.make_association(assoc)

        return

    def make_association(self, record):
        """
        contstruct the association
        :param record:
        :return: modeled association of  genotype to mammalian phenotype
        """
        # prep record
        # remove description and mapp Experiment Type to apo term
        experiment_type = record['Experiment Type'].split('(')[0]
        experiment_type = experiment_type.split(',')
        record['experiment_type'] = list()
        for exp_type in experiment_type:
            exp_type = exp_type.lstrip().rstrip()
            record['experiment_type'].append(
                {
                    'id': self.apo_term_id[exp_type],
                    'term': exp_type,
                })
        sgd_phenotype = record['Phenotype']
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
        phenotype = record['Phenotype']
        if ':' in phenotype:
            pheno_obj['has_quality'] = True
            ent_qual = sgd_phenotype.split(': ')
            entity = ent_qual[0]
            quality = ent_qual[1]
            pheno_obj['entity']['term'] = entity
            pheno_obj['entity']['apo_id'] = self.apo_term_id[entity]
            pheno_obj['quality']['term'] = quality
            pheno_obj['quality']['apo_id'] = self.apo_term_id[quality]
        else:
            pheno_obj['entity']['term'] = phenotype
            pheno_obj['entity']['apo_id'] = self.apo_term_id[phenotype]
        record['pheno_obj'] = pheno_obj

        # begin modeling
        model = Model(self.graph)

        # define the triple
        gene = 'SGD:{}'.format(record['SGDID'])
        relation = Model.object_properties['has_phenotype']  # has phenotype

        if record['pheno_obj']['has_quality']:
            pheno_label = '{0}:{1}'.format(
                record['pheno_obj']['entity']['term'],
                record['pheno_obj']['quality']['term'])
            pheno_id = 'MONARCH:{0}{1}'.format(
                record['pheno_obj']['entity']['apo_id'].replace(':', '_'),
                record['pheno_obj']['quality']['apo_id'].replace(':', '_')
            )
            g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=pheno_id, pred=relation)
        else:
            pheno_label = record['pheno_obj']['entity']['term']
            pheno_id = record['pheno_obj']['entity']['apo_id']
            g2p_assoc = Assoc(self.graph, self.name, sub=gene, obj=pheno_id, pred=relation)
            assoc_id = g2p_assoc.make_association_id(definedby='yeastgenome.org', subject=gene, predicate=relation,
                                                     object=pheno_id)
            g2p_assoc.set_association_id(assoc_id=assoc_id)

        # add to graph to mint assoc id
        g2p_assoc.add_association_to_graph()

        model.addLabel(subject_id=gene, label=record['Gene Name'])

        # add the association triple
        model.addTriple(subject_id=gene, predicate_id=relation, obj=pheno_id)

        # make pheno subclass of UPHENO:0001001
        model.addTriple(subject_id=pheno_id, predicate_id=Model.object_properties['subclass_of'], obj='UPHENO:0001001')

        # label nodes
        # pheno label
        model.addLabel(subject_id=pheno_id, label=pheno_label)


        # add the descripiton: all the unmodeled data in a '|' delimited list
        description = [
            'genomic_background: {}'.format(record['Strain Background']),
            'allele: {}'.format(record['Allele']),
            'chemical: {}'.format(record['Chemical']),
            'condition: {}'.format(record['Condition']),
            'details: {}'.format(record['Details']),
            'feature_name: {}'.format(record['Feature Name']),
            'gene_name: {}'.format(record['Gene Name']),
            'mutant_type: {}'.format(record['Mutant Type']),
            'reporter: {}'.format(record['Reporter']),
        ]
        g2p_assoc.description = " | ".join(description)

        # add the references
        references = record['Reference']
        references = references.replace(' ', '')
        references = references.split('|')

        #  created RGDRef prefix in curie map to route to proper reference URL in RGD
        if len(references) > 0:
            # make first ref in list the source
            g2p_assoc.add_source(identifier=references[0])
            ref_model = Reference(
                self.graph, references[0],
                Reference.ref_types['publication']
            )
            ref_model.addRefToGraph()

        if len(references) > 1:
            # create equivalent source for any other refs in list
            for ref in references[1:]:
                model.addSameIndividual(sub=references[0], obj=ref)

        # add experiment type as evidence
        for exp_type in record['experiment_type']:
            g2p_assoc.add_evidence(exp_type['id'])
            model.addLabel(subject_id=exp_type['id'], label=exp_type['term'])

        try:
            g2p_assoc.add_association_to_graph()
        except Exception as e:
            print(e)
        return

    @staticmethod
    def make_apo_map():
        # load apo for term mapping
        ofactory = OntologyFactory()
        apo_ont = ofactory.create("apo")
        apo_nodes = apo_ont.nodes()
        # dict schema { 'term': 'apo_id' }
        apo_term_id = dict()
        for node in apo_nodes:
            label = apo_ont.label(node)
            apo_term_id[label] = node
        return apo_term_id
