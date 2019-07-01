from enum import Enum


class BioLinkVocabulary(Enum):
    category = "biolink:category" # this is a node property, not a class
    publication = "biolink:Publication"
    disease = "biolink:Disease"
    gene = "biolink:Gene"
    genotype = "biolink:Genotype"
    namedThing = "biolink:NamedThing"
    population_of_individual_organisms = "biolink:PopulationOfIndividualOrganisms"
