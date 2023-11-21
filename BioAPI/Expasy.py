import os, sys
from typing import Dict
from pathlib import Path
import requests
from Bio import ExPASy
from Bio import SwissProt


class ExpasyApi:
    def __init__(self) -> None:
        self.records: Dict[str, SwissProt.Record] = {}

    def fetch_from_id(self, accession: str) -> SwissProt.Record:
        """
        Retrieves a protein record from the ExPASy server from an accession id.

        Parameters
        ==========
        >> accession (str): an alphanumeric string that maps a record to the database

        Returns
        =======
        A SwissProt.Record object containing protein data
        """
        handle = ExPASy.get_sprot_raw(accession)
        record = SwissProt.read(handle)
        self.records[accession] = record
        handle.close()
        return record

    def save(self, accession: str, filename: str, dir: Path = None) -> Path:
        """
        Retrieves a record directly from the ExPASy RESTful API. Writes the record locally.

        Parameters
        ==========
        >> accession  (str)           : an alphanumeric string that maps to the record
        >> filename   (str)           : name for the record
        >> dir        (Path, optional): base path of the target directory; creates as `BioAPI_data` directory at
                                        the HOME/USERPROFILE path if not provided

        Returns
        =======
        (Path) A filepath for the saved record
        """
        BASE_URL = "https://rest.uniprot.org/uniprotkb/"
        url = BASE_URL + accession + ".txt"
        res = requests.get(url)
        res.raise_for_status()
        # set filepath to default (located in user's $HOME variable) if `dir` is not provided
        filename = Path(filename + ".txt")
        if dir:
            filepath = dir / filename
        else:
            HOME_ENV = "HOME" if sys.platform == "linux" else "USERPROFILE"
            HOME = Path(os.environ.get(HOME_ENV))
            default_dir = HOME / "BioAPI_data"
            if not os.path.isdir(default_dir):
                os.mkdir(default_dir)
            filepath = default_dir / filename

        write_handle = open(filepath, "w")
        for line in res.text.split("\n"):
            write_handle.write(line.rstrip() + "\n")
        write_handle.close()
        return filepath
