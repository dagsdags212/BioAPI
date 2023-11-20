from enum import Enum


"""
A complete list of valid Entrez databases. More details can be found on the official NCBI website:
>> https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
"""
class EntrezDB(Enum):
    BIOPROJECT          = "bioproject"
    BIOSAMPLE           = "biosample"
    BOOKS               = "books"
    CONSERVED_DOMAINS   = "cdd"
    DBGAP               = "gap"
    DBVAR               = "dbvar"
    GENE                = "gene"
    GENOME              = "genome"
    GEO_DATASETS        = "gds"
    GEO_PROFILES        = "geoprofiles"
    HOMOLOGENE          = "homologene"
    MESH                = "mesh"
    NCBI_TOOLKIT        = "toolkit"
    NLM_CATALOG         = "nlmcatalog"
    NUCLEOTIDE          = "nuccore"
    POPSET              = "popset"
    PROBE               = "probe"
    PROTEIN             = "protein"
    PROTEIN_CLUSTERS    = "proteinclusters"
    PUBCHEM_BIOASSAY    = "pcassay"
    PUBCHEM_COMPOUND    = "pccompound"
    PUBCHEM_SUBSTANCE   = "pcsubstance"
    PUBMED              = "pubmed"
    PUBMED_CENTRAL      = "pmc"
    SNP                 = "snp"
    SRA                 = "sra"
    STRUCTURE           = "structure"
    TAXONOMY            = "taxonomy"

    def __repr__(self) -> str:
        return self.value

    def __str__(self) -> str:
        return self.name
