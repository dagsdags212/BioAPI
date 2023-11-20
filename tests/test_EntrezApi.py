import unittest
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BioAPI import Credentials, EntrezApi
from BioAPI.Enums import *

class TestEntrezApiClass(unittest.TestCase):

    def setUp(self):
        email, api_key = os.environ.get("EMAIL"), os.environ.get("API_KEY")
        self.creds = Credentials(email, api_key, "BioAPI library")
        self.app = EntrezApi(self.creds)

    def test_send_credentials(self):
        """Expect `_send_credentials` to return None"""
        self.assertIsNone(self.app._send_credentials())

    def test_fetch_by_id_return_type(self):
        """Expect `fetch_by_id` to return a SeqRecord object"""
        accession = "EU490707"
        db = EntrezDB.NUCLEOTIDE
        rettype = Rettype.GENBANK
        retmode = Retmode.TEXT
        record = self.app.fetch_by_id(accession, db, rettype, retmode)
        self.assertIsInstance(record, SeqRecord)
        self.assertIsInstance(record.seq, Seq)
