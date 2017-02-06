from dipper.models.Reference import Reference
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.sources.Source import Source
from dipper.sources.ZFIN import ZFIN
from dipper.models.Dataset import Dataset
from dipper.models.Model import Model
import csv
import logging

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

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'zfin_slim')
        self.dataset = Dataset(
            'zfin_slim', 'ZFINSlim', 'http://zfin.org/')

    def fetch(self, is_dl_forced=False):
        self.get_files(is_dl_forced)

    def parse(self, limit=None):
        zfin_parser = ZFIN(self.graph_type, self.are_bnodes_skized)
        model = Model(self.graph)
        zp_file = '/'.join((self.rawdir, self.files['zpmap']['file']))
        g2p_file = '/'.join((self.rawdir, self.files['g2p_clean']['file']))
        zfin_parser.zp_map = zfin_parser._load_zp_mappings(zp_file)

        with open(g2p_file, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:

                (internal_id, symbol, gene_id, subterm1_id, subterm1_label,
                 pc_rel_id, pc_rel_label, superterm1_id, superterm1_label,
                 quality_id, quality_name, modifier, subterm2_id,
                 subterm2_label, pc_rel2_id, pc_rel2_id, superterm2_id,
                 superterm2_label, fish_id, fish_label, start_stage, end_stage,
                 environment, pub_id, figure_id, unknown_field) = row

                zp_id = zfin_parser._map_sextuple_to_phenotype(
                    superterm1_id, subterm1_id, quality_id, superterm2_id,
                    subterm2_id, modifier)

                gene_curie = "ZFIN:{0}".format(gene_id)
                model.makeLeader(gene_curie)
                pub_curie = "ZFIN:{0}".format(pub_id)
                if zp_id:
                    assoc = G2PAssoc(self.graph, self.name, gene_curie, zp_id)
                    if pub_id:
                        reference = Reference(self.graph, pub_curie,
                                              Reference.ref_types['document'])
                        reference.addRefToGraph()
                        assoc.add_source(pub_curie)

                    assoc.add_evidence('ECO:0000059')
                    assoc.add_association_to_graph()



