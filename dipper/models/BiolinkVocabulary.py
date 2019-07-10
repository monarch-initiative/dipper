from enum import Enum

"""
This is a class to support categorizing everything we ingest using biolink 
categories.

If you need to add an Enum for something that isn't already here,
consult this YAML file:
https://github.com/biolink/biolink-model/blob/master/biolink-model.yaml
"""

class BioLinkVocabulary(Enum):
    category = "biolink:category"  # this is a node property, not a class
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
    OrganismTaxon = "biolink:OrganismTaxon"
    MolecularEntity = "biolink:MolecularEntity"
    GenomicEntity = "biolink:GenomicEntity"  # we are using this for chromosome
    Case = "biolink:Case"
    OntologyClass = "biolink:OntologyClass"
    GeneFamily = "biolink:GeneFamily"
    Transcript = "biolink:Transcript"
    Protein = "biolink:Protein"
    ChemicalSubstance = "biolink:ChemicalSubstance"
    GenomicSequenceLocalization = "biolink:GenomicSequenceLocalization"
    IndividualOrganism = "biolink:IndividualOrganism"
    GeneGrouping = "biolink:GeneGrouping"
    Genome = "biolink:Genome"
    GenomeBuild = "biolink:GenomeBuild"
    BiologicalProcess = "biolink:BiologicalProcess"
    LifeStage = "biolink:LifeStage"
    Environment = "biolink:Environment"
    EvidenceType = "biolink:EvidenceType"
    AnatomicalEntity = "biolink:AnatomicalEntity"
    Provider = "biolink:Provider"
    Procedure = "biolink:Procedure"
    Association = "biolink:Association"
