import os
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from .Enums import *


class Credentials:
    def __init__(self, email: str = None, api_key: str = None, tool: str = None) -> None:
        self.email = email if email else os.environ.get("EMAIL")
        self.api_key = api_key if api_key else os.environ.get("API_KEY")
        self.tool = tool if tool else None


class EntrezApi:
    def __init__(self, email: str = None, api_key: str = None, tool: str = None) -> None:
        """Stores user credentials upon instantiation, if given by the user."""
        self.credentials = Credentials(email=email, api_key=api_key, tool=tool)

    def _send_credentials(self) -> None:
        """Sends user credentials to the Entrez server. Must be done for every request."""
        Entrez.email = self.credentials.email
        Entrez.api_key = self.credentials.api_key
        Entrez.tool = self.credentials.tool

    def fetch_by_id(
        self,
        accession: str,
        db: EntrezDB,
        rettype: Rettype = Rettype.GENBANK,
        retmode: Retmode = Retmode.TEXT) -> SeqRecord:
        """
        Sends a request to the Entrez server to fetch a record from the given accession id.

        Parameters
        ==========
        >> accession (str) : an alphanumeric string that maps to the record
        >> db        (Enum): a valid Entrez databae
        >> rettype   (Enum): a valid return type (default = gb)
        >> retmode   (Enum): a valid return mode (deefault = text)

        Returns
        =======
        A SeqRecord object containing the record data
        """
        self._send_credentials()
        handle = Entrez.efetch(db=db, id=accession, rettype=rettype, retmode=retmode)
        fformats = {
            "fasta": "fasta",
            "gb": "genbank",
        }
        record = SeqIO.read(handle, fformats.get(rettype, "gb"))
        handle.close()
        return record
