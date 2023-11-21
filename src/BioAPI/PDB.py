import os
import sys
from typing import Optional
from pathlib import Path
import requests
from Bio import PDB


class PDBApi:

    def __init__(self) -> None:
        self.records = {}

    @classmethod
    def save_record(self, accession: str, filename: Optional[str] = None, dir: Optional[Path] = None) -> Path:
        """
        Retrieves a PDB record from the RCSB server from an accession id.

        Parameters
        ==========
        >> accession (str)          : an alphanumeric string that maps a record to the database
        >> filename  (Optional, str): name for the PDB record
        >> dir       (Optional, str): base path of the target directory; creates a `BioAPI_data` directory at
                                      the HOME/USERPROFILE path if not provided

        Returns
        =======
        A Path object pointing to the saved PDB record
        """
        assert len(accession) >= 4, "`accession` must be of length 4 or greater"
        url = f"https://files.rcsb.org/download/{accession}.pdb"
        res = requests.get(url)
        res.raise_for_status()
        # set filepath to default (located in user's $HOME variable) if `dir` is not provided
        filename = Path(f"{filename}.pdb") if filename else Path(f"{accession}.pdb")
        if dir:
            filepath = dir / filename
        else:
            HOME_ENV = "HOME" if sys.platform == "linux" else "USERPROFILE"
            HOME = Path(os.environ.get(HOME_ENV))
            default_dir = HOME / "BioAPI_data"
            if not os.path.isdir(default_dir):
                os.mkdir(default_dir)
            filepath = default_dir / filename
            print(f"Parameter `dir` not provided. Saving file to {default_dir}")
        # append filepath to records attribute
        self.records[accession] = filepath
        # write response text to filepath
        write_handle = open(filepath, "w")
        for line in res.text.split("\n"):
            write_handle.write(line.rstrip() + "\n")
        write_handle.close()
        print(f"Saved `{filename}` at {filepath}")
        return filepath

