from enum import Enum, auto


class BioLinkVocabulary(Enum):
    publication = "biolink:publication"
    disease = "biolink:disease"
    gene = "biolink:gene"
    namedThing = "biolink:NamedThing"
    category = "biolink:category"
