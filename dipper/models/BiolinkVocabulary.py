from enum import Enum


class BioLinkVocabulary(Enum):
    publication = "biolink:publication"
    disease = "biolink:disease"
    gene = "biolink:gene"
    namedThing = "biolink:NamedThing"
    category = "biolink:category"
