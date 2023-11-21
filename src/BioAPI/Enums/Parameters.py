from enum import Enum, EnumMeta


class MyEnum(EnumMeta):
    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.value

class FileFormat(MyEnum):
    """Supported file formats for parsing files dervide from the EFetch Entrez utility"""
    GENBANK = "genbank"

class Rettype(MyEnum):
    """Supported return types (rettype)"""
    FASTA   = "fasta"
    GENBANK = "gb"

class Retmode(MyEnum):
    """Supported return modes (retmode)"""
    TEXT    = "text"


