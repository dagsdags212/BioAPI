import unittest
import os
from BioAPI.Fasta import Credentials


class TestCredentailsClass(unittest.TestCase):

    def setUp(self):
        email = os.environ.get("EMAIL")
        api_key = os.environ.get("API_KEY")
        tool = "BioAPI module"
        self.creds = Credentials(email, api_key, tool)

    def test_attributes(self):
        self.assertEqual(self.creds.tool, "BioAPI module")

if __name__ == "__main__":
    unittest.main()
