import logging
import unittest
import math

from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder

logger = logging.getLogger('ARGDrugTablePointfinderTest')


class ARGDrugTablePointfinderTest(unittest.TestCase):

    def setUp(self):
        self.arg_drug_table = ARGDrugTablePointfinder()

    def testNoneEntry(self):
        # Tests when the entry for the drug is "None"
        drug = self.arg_drug_table.get_drug("escherichia_coli", "parC", 57)

        # Specifically, we're interested in it not crashing, and the current behaviour is to return NaN.
        print("The drug is: " + str(drug))
        self.assertTrue(math.isnan(drug))
