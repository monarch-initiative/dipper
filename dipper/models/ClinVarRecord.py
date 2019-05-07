"""
https://www.ncbi.nlm.nih.gov/clinvar/docs/details/
Object mapping to XML schema
"""

from typing import List, Optional, Union


class Gene:
    """
    ClinVar Gene
    Intentionally leaves out label/symbol, this should come from HGNC
    """
    def __init__(self,
                 id: Union[str,int,None],
                 association_to_allele: str):
        self.id = id
        self.association_to_allele = association_to_allele


class Allele:
    """
    ClinVar Allele
    Alleles can have 0 to many genes

    These are called alleles and variants on the ClinVar UI,
    and variant, single nucleotide variant, etc in the XML

    id: allele id
    label: label
    variant_type: single nucleotide variant
    genes: gene(s) in which the variant is found
    synonyms: eg HGVC
    dbsnp: dbSNP curies
    """

    def __init__(self,
                 id: str,
                 label: Optional[str] = None,
                 variant_type: Optional[str] = None,
                 genes: Optional[List[Gene]] = None,
                 synonyms: Optional[List[str]] = None,
                 dbsnps: Optional[List[str]] = None):
        self.id = id
        self.label = label
        self.variant_type = variant_type
        self.genes = genes if genes is not None else []
        self.synonyms = synonyms if synonyms is not None else []
        self.dbsnps = dbsnps if dbsnps is not None else []


class Genovar:
    """
    Sequence feature entity that is linked to a disease
    in an Reference ClinVar Record, it is either a
    Variant or Genotype
    """
    def __init__(self,
                 id: str,
                 label: Optional[str] = None,
                 variant_type: Optional[str] = None):
        self.id = id
        self.label = label
        self.variant_type = variant_type


class Variant(Genovar):
    """
    ClinVar variant, variants can have one or more alleles

    These are called variants and alleles on the ClinVar UI, and
    Variants in the XML
    """
    def __init__(self,
                 id: str,
                 label: Optional[str] = None,
                 alleles: Optional[List[Allele]] = None,
                 variant_type: Optional[str] = None):
        super().__init__(id, label, variant_type)
        self.alleles = alleles if alleles is not None else []


class Genotype(Genovar):
    """
    ClinVar genotype
    Example: Compound Heterozygote, Diplotype

    These are called variants on the ClinVar UI, and a GenotypeSet
    in the XML
    """
    def __init__(self,
                 id: str,
                 label: Optional[str] = None,
                 variants: Optional[List[Variant]] = None,
                 variant_type: Optional[str] = None):
        super().__init__(id, label, variant_type)
        self.variants = variants if variants is not None else []


class Condition:
    """
    ClinVar condition
    """
    def __init__(self,
                 id: str,
                 label: Optional[str] = None,
                 database: Optional[str] = None,
                 medgen_id: Optional[str] = None):
        self.id = id
        self.label = label
        self.database = database
        self.medgen_id = medgen_id


class ClinVarRecord:
    """
    Reference ClinVar Record (RCV)
    id: RCV id
    accession: RCV accession (eg RCV000123456)
    created: Created date
    updated: Updated date
    genovar: the variant or genotype associated with the condition(s)
    condition: The condition(s) for which this allele set was interpreted, with
               links to databases with defining information about that condition.
    """
    def __init__(self,
                 id: str,
                 accession: str,
                 created: str,
                 updated: str,
                 genovar: Genovar,
                 conditions: Optional[List[Condition]] = None):
        self.id = id
        self.accession = accession
        self.created = created
        self.updated = updated
        self.genovar = genovar
        self.conditions = conditions if conditions is not None else []
