import logging
import unittest
import pandas

from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder

logger = logging.getLogger('ARGDrugTablePointfinderTest')


class ARGDrugTablePointfinderTest(unittest.TestCase):

    def setUp(self):
        self.arg_drug_table = ARGDrugTablePointfinder()

    def testNoneEntry(self):
        # Tests when the entry for the drug is "None"
        drug = self.arg_drug_table.get_drug("escherichia_coli", "parC", 57)

        # Specifically, we're interested in it not crashing.
        self.assertTrue(pandas.isna(drug))
