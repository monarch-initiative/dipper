from enum import Enum


class BioLinkVocabulary(Enum):
    category = "biolink:category"  # this is a node property, not a slot
    Publication = "biolink:Publication"
    Disease = "biolink:Disease"
    Gene = "biolink:Gene"
    Genotype = "biolink:Genotype"
    NamedThing = "biolink:NamedThing"
    PopulationOfIndividualOrganisms = "biolink:PopulationOfIndividualOrganisms"
    SequenceVariant = "biolink:SequenceVariant"  # synonym for allele
    InformationContentEntity = "biolink:InformationContentEntity"
    Zygosity = "biolink:Zygosity"
    PhenotypicFeature = "biolink:PhenotypicFeature"
    Publications = "biolink:Publications"
    BiologicalSex = "biolink:BiologicalSex"