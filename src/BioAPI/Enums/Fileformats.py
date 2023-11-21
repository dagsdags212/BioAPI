from enum import Enum


class Fileformat:
    """A list of file formats capable of being parsed by the SeqIO module"""
    ABI             = "abi"
    ABI_TRIM        = "abi-trim"
    ACE             = "ace"
    CIF_ATOM        = "cif-atom"
    CIF_SEQRES      = "cif-seqres"
    CLUSTAL         = "clustal"
    EMBL            = "embl"
    FASTA           = "fasta"
    FASTA_2LINE     = "fasta-2line"
    FASTQ           = "fastq"
    FASTQ_SOLEXA    = "fastq-solexa"
    FASTQ_ILLUMINA  = "fastq-illumina"
    GCK             = "gck"
    GENBANK         = "genbank"
    IG              = "ig"
    IMGT            = "imgt"
    NEXUS           = "nexus"
    PDB_SEQRES      = "pdb-seqres"
    PDB_ATOM        = "pdb-atom"
    PHD             = "phd"
    PHYLIP          = "phylip"
    PIR             = "pir"
    SEQXML          = "seqxml"
    SFF             = "sff"
    SFF_TRIM        = "sff-trim"
    SNAPGENE        = "snapgene"
    STOCKHOLM       = "stockholm"
    TAB             = "tab"
    QUAL            = "qual"
    UNIPROT_XML     = "uniprot-xml"
    XDNA            = "xdna"

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.value
