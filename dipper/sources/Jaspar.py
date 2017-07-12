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
import cvs

logger = logging.getLogger(__name__)


class Jaspar(Source):
    """
    Putative upstream binding site motifs
    
    """

    UCSC = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/'
    NCBI = 'ftp://ftp.ncbi.nih.gov'
    files = {
        'file': 'GENE_JASPAR_bucket_count.tab',
        'urls' = {
             # the remote files to create 'GENE_JASPAR_bucket_count.tab' 
             'url': 'http://jaspar.genereg.net/html/DOWNLOAD/bed_files/MA*.bed',
             '1k': UCSC + '/upstream1000.fa.gz',
             '2k': UCSC + '/upstream2000.fa.gz',
             '5k': UCSC + '/upstream5000.fa.gz',
             'rsng':  NCBI + '/gene/DATA/gene2refseq.gz' # 500M
        }     
    }

    def __init__(self, graph_type, are_bnodes_skolemized):
        super().__init__(graph_type, are_bnodes_skolemized, 'jaspar')

        self.dataset = Dataset(
            'Jaspar', 'Putative binding site motifs',
            'http://jaspar.genereg.net/')

    def fetch(self, is_dl_forced=False):
        """
        Unaware of any update schedule for this dataset,
        Readme is four years old but bed files refreshed a month ago.
        Since they are still mapping to hg19
        I am not expecting urgent bleeding edge update behavior.
        also their webserver is a bit picky, chokes on wget requests,
        accepts curl requests
        
        :param is_dl_forced: boolean, force download
        :return:
        """
        return
        

    def parse(self, limit=None):
        """
        Parse jaspar file
        :param limit: int limit rows processed

        :return: None
        """
        dir_path = Path(self.rawdir)
        jaspar_file = dir_path / self.files['file']
        with open(jaspar_file, 'r', encoding="utf8") as csvfile:
            jaspar_fh = csv.reader(csvfile, delimiter='\t')
        count = 0
        model = Model(self.graph)
        for line in jaspar_fh.readlines():
            if limit is not None and count >= limit:
                break
            (ncbigene, motif, bucket, howmany) = line
            count += 1

            motif_label = ncbigene + ' ' + motif + ' ' + bucket
            bnode = '_:' + digest_id(motif_label)

            # <ncbigene><SO:five_prime_flanking_region><bnode>
            self.graph.addTriple(
                ncbigene,
                resolve('SO:five_prime_flanking_region', local_tt),
                bnode,
                object_is_literal=False)
                
            
            
            spo(bnode, 'rdf:type', 'SO:TF_binding_site') 
            spo(bnode, 'rdfs:label', motif_label) 
            spo(bnode, 'rdf:value', count) 
            spo(bnode, 'GENO:has_extent', bucket) 
            spo(bnode, 'OIO:has_xbxref', motif)

            self.graph.addTriple(bnode,
                             model.annotation_properties['inchi_key'],
                             document['unii']['inchikey'],
                             object_is_literal=True)
           

            
            if count % 50000 == 0:
                    logger.info("Processed {} gene_motif_buckets".format(count))

        jaspar_fh.close()
        return



    def _parse_jaspar_data(self, document, or_limit=None):
        model = Model(self.graph)

        rxcui_curie = "RXCUI:{}".format(document['jaspar']['rxcui'])
        uni_curie = "UNII:{}".format(document['jaspar']['unii'])
        model.addLabel(rxcui_curie, document['jaspar']['drug_name'])
        model.addLabel(uni_curie, document['jaspar']['drug_name'])

        model.addSameIndividual(rxcui_curie, uni_curie)
        
        self.graph.addTriple(rxcui_curie,
                             model.annotation_properties['inchi_key'],
                             document['unii']['inchikey'],
                             object_is_literal=True)



        for outcome in outcomes:
            drug2outcome_assoc = Assoc(self.graph, self.name)

            meddra_curie = "MEDDRA:{}".format(outcome['code'])
            model.addLabel(meddra_curie, outcome['name'])

            drug2outcome_assoc.sub = rxcui_curie
            drug2outcome_assoc.obj = meddra_curie
            drug2outcome_assoc.rel = \
                Assoc.object_properties['causes_or_contributes']
            drug2outcome_assoc.description = \
                "A proportional reporting ratio or odds " \
                "ratio greater than or equal to {} in the " \
                "jaspar data was the significance cut-off " \
                "used for creating drug-outcome associations".format(or_limit)
            drug2outcome_assoc.add_association_to_graph()
            drug2outcome_assoc.add_predicate_object(
                Assoc.annotation_properties['probabalistic_quantifier'],
                outcome['ror'], 'Literal')

            self._add_outcome_evidence(drug2outcome_assoc.assoc_id, outcome)
            self._add_outcome_provenance(drug2outcome_assoc.assoc_id, outcome)



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
