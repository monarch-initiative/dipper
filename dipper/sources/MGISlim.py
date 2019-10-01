import logging
import datetime

from intermine.webservice import Service
from dipper.models.assoc.G2PAssoc import G2PAssoc
from dipper.sources.Source import Source
from dipper.models.Reference import Reference
from dipper.models.BiolinkVocabulary import BioLinkVocabulary as blv

LOG = logging.getLogger(__name__)


class MGISlim(Source):
    """
    slim mgi model only containing Gene to phenotype associations
    Uses mousemine: http://www.mousemine.org/mousemine/begin.do
    python lib api  http://intermine.org/intermine-ws-python/

    """

    def __init__(self,
                 graph_type,
                 are_bnodes_skolemized,
                 data_release_version=None):
        super().__init__(
            graph_type=graph_type,
            are_bnodes_skized=are_bnodes_skolemized,
            data_release_version=data_release_version,
            name='mgi_slim',
            ingest_description='Simplified Mouse Genome Informatics',
            ingest_url='http://www.mousemine.org/mousemine/service',
            ingest_logo="source-mgi.png",
            license_url='http://www.informatics.jax.org/mgihome/other/copyright.shtml',
            # data_rights=None,
            # file_handle=None
        )

        # "NCBITaxon:
        self.txid = self.globaltt['Mus musculus'][10:]

    def fetch(self, is_dl_forced=False):
        return

    def parse(self, limit=None):

        count = 0
        for num in range(10, 100):
            fuzzy_gene = "MGI:{0}*".format(num)
            gene = "MGI:{0}".format(num)
            service = Service("http://www.mousemine.org/mousemine/service")
            logging.getLogger('Model').setLevel(logging.ERROR)
            logging.getLogger('JSONIterator').setLevel(logging.ERROR)
            query = service.new_query("OntologyAnnotation")
            query.add_constraint("subject", "SequenceFeature")
            query.add_constraint("ontologyTerm", "MPTerm")
            query.add_view(
                "subject.primaryIdentifier", "subject.symbol",
                "subject.sequenceOntologyTerm.name", "ontologyTerm.identifier",
                "ontologyTerm.name", "evidence.publications.pubMedId",
                "evidence.comments.type", "evidence.comments.description"
            )
            query.add_sort_order("OntologyAnnotation.ontologyTerm.name", "ASC")
            query.add_constraint("subject.organism.taxonId", "=", self.txid, code="A")
            query.add_constraint("subject", "LOOKUP", fuzzy_gene, code="B")
            query.add_constraint(
                "subject.primaryIdentifier", "CONTAINS", gene, code="C")
            query.outerjoin("evidence.comments")

            for row in query.rows():
                mgi_curie = row["subject.primaryIdentifier"]
                mp_curie = row["ontologyTerm.identifier"]
                pub_curie = "PMID:{0}".format(row["evidence.publications.pubMedId"])
                assoc = G2PAssoc(self.graph, self.name, mgi_curie, mp_curie,
                                 entity_category=blv.terms.Gene.value)
                if row["evidence.publications.pubMedId"]:
                    reference = Reference(
                        self.graph, pub_curie, self.globaltt['journal article'])
                    reference.addRefToGraph()
                    assoc.add_source(pub_curie)

                assoc.add_evidence(self.globaltt['experimental phenotypic evidence'])
                assoc.add_association_to_graph()

            if not count % 10 and count != 0:
                count_from = count - 10
                LOG.info(
                    "%s processed ids from MGI:%i* to MGI:%i*",
                    datetime.datetime.now(), count_from, count)

            count += 1
            if limit and count >= limit:
                break

        return
