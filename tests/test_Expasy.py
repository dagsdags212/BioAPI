import unittest
import sys, os, subprocess
from pathlib import Path
from Bio import SwissProt
from BioAPI import *


class TestExpasyApiClass(unittest.TestCase):

    def setUp(self):
        self.app = ExpasyApi()
        self.acc = "G3ECR1"
        HOME_ENV = "HOME" if sys.platform == "linux" else "USERPROFILE"
        HOME = Path(os.environ.get(HOME_ENV))
        self.default_dir = HOME / "BioAPI_data"

    def test_fetch_by_id_return_type(self):
        """Expect `fetch_by_id` to return a SwissProt.Record object"""
        record = self.app.fetch_from_id(self.acc)
        self.assertIsInstance(record, SwissProt.Record)

    def test_save_without_provided_dir(self):
        """Expect that a default directory (`BioAPI_data`) will be created"""
        path = self.app.save(self.acc, self.acc+".txt")
        self.assertTrue(os.path.isdir(self.default_dir))

    def test_save_output(self):
        """Expect the first line to contain `ID` and the second line to contain `AC`"""
        path = self.app.save(self.acc, self.acc+".txt")
        self.assertIsInstance(path, Path)
        with open(path, "r") as fhandle:
            line1 = fhandle.readline()
            line2 = fhandle.readline()
        fhandle.close()
        self.assertIn("ID", line1)
        self.assertIn("AC", line2)
