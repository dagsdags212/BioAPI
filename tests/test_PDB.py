import unittest
import sys, os, subprocess
from pathlib import Path
from BioAPI import *


class TestPDBApiClass(unittest.TestCase):

    def setUp(self):
        self.app = PDBApi()
        self.acc = "4CUP"
        HOME_ENV = "HOME" if sys.platform == "linux" else "USERPROFILE"
        HOME = Path(os.environ.get(HOME_ENV))
        self.default_dir = HOME / "BioAPI_data"

    def test_save_record_return_type(self):
        """Expect `save_record` to return a Path object"""
        path = self.app.save_record(self.acc)
        self.assertIsInstance(path, Path)
        os.remove(path)

    def test_save_record_without_provided_dir(self):
        """Expect that a default directory (`BioAPI_data`) will be created"""
        path = self.app.save_record(self.acc)
        self.assertTrue(os.path.isdir(self.default_dir))
        os.remove(path)

    def test_save_record_output_file(self):
        """
        Checks for the presence of:
            > "HEADER" in first line
            > "TITLE" in the second line
            > "COMPND" in the third line
        """
        path = self.app.save_record(self.acc)
        self.assertIsInstance(path, Path)
        with open(path, "r") as fhandle:
            line1 = fhandle.readline()
            line2 = fhandle.readline()
            line3 = fhandle.readline()
        fhandle.close()
        self.assertIn("HEADER", line1)
        self.assertIn("TITLE", line2)
        self.assertIn("COMPND", line3)
        os.remove(path)
