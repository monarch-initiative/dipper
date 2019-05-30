from enum import Enum, auto

class BioLinkVocabulary(Enum):
    # need an IRI to write these out as triples, so let's use these:
    # https://github.com/monarch-initiative/scigraph-docker/blob/master/conf/monarch-load-global.yaml#L1
    # probably will need to and should refactor this
    publication = 'http://purl.obolibrary.org/obo/IAO_0000310'
    disease = 'http://purl.obolibrary.org/obo/MONDO_0000001'
    gene = 'http://purl.obolibrary.org/obo/SO_0000704'
    phenotypic_feature = 'http://purl.obolibrary.org/obo/UPHENO_0001001'
