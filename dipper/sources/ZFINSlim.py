import csv
import logging

from dipper.models.Reference import Reference
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.sources.Source import Source
from dipper.sources.ZFIN import ZFIN
from dipper.models.Model import Model


LOG = logging.getLogger(__name__)
# note: currently no log issued


class ZFINSlim(Source):
    """
    zfin mgi model only containing Gene to phenotype associations
    Using the file here: https://zfin.org/downloads/phenoGeneCleanData_fish.txt
    """
    files = {
        'g2p_clean': {
            'file': 'phenoGeneCleanData_fish.txt',
            'url': 'https://zfin.org/downloads/phenoGeneCleanData_fish.txt',
            # https://zfin.org/downloads#  header Documentation is burried in UI crap
        },
        'zpmap': {
            'file': 'zp-mapping-2019.txt',
            'url': 'http://purl.obolibrary.org/obo/zp/src/curation/id_map_zfin.tsv'
                   # ^^ Nico's updated mapping, May 2019
        }
    }

    def __init__(self, graph_type, are_bnodes_skolemized, skip_stats=False):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            name='zfinslim',
            ingest_title='Simplified ZFIN',
            ingest_url='https://zfin.org/',
            ingest_logo="https://github.com/monarch-initiative/monarch-ui/blob/master/public/img/sources/source-zfin.png",
            license_url=None,
            data_rights='http://zfin.org/warranty.html',
            # file_handle=None
        )
        self.dataset.set_citation(
            'https://wiki.zfin.org/display/general/ZFIN+db+information')

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
                (internal_id,
                 symbol,
                 gene_id,
                 subterm1_id,
                 subterm1_label,
                 pc_rel_id,
                 pc_rel_label,
                 superterm1_id,
                 superterm1_label,
                 quality_id,
                 quality_name,
                 modifier,
                 subterm2_id,
                 subterm2_label,
                 pc_rel2_id,
                 pc_rel2_label,
                 superterm2_id,
                 superterm2_label,
                 fish_id,
                 fish_label,
                 start_stage,
                 end_stage,
                 environment,
                 pub_id,
                 figure_id
                ) = row

                if modifier != "abnormal":
                    LOG.warning("skipping phenotype with modifier != abnormal: " + modifier)
                    continue

                zp_id = zfin_parser._map_octuple_to_phenotype(subterm1_id,
                                                              pc_rel_id,
                                                              superterm1_id,
                                                              quality_id,
                                                              subterm2_id,
                                                              pc_rel2_id,
                                                              superterm2_id,
                                                              modifier)

                gene_curie = "ZFIN:{0}".format(gene_id)
                model.makeLeader(gene_curie)
                pub_curie = "ZFIN:{0}".format(pub_id)
                if zp_id:
                    assoc = G2PAssoc(self.graph, self.name, gene_curie, zp_id)
                    if pub_id:
                        reference = Reference(self.graph, pub_curie,
                                              self.globaltt['document'])
                        reference.addRefToGraph()
                        assoc.add_source(pub_curie)

                    assoc.add_evidence(
                        self.globaltt['experimental phenotypic evidence'])
                    assoc.add_association_to_graph()
