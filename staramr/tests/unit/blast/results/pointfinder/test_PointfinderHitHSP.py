import unittest

import pandas as pd

from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.exceptions.InvalidPositionException import InvalidPositionException


class PointfinderHitHSPTest(unittest.TestCase):

    def testBuildPointfinderHitHSPSuccess(self):
        blast_record = pd.Series({'sstart': 1, 'send': 10, 'qstart': 1, 'qend': 10, 'sstrand': 'plus'})

        PointfinderHitHSP(file=None, blast_record=blast_record)

    def testBuildPointfinderHitHSPFailInvalidSubjectCoords(self):
        blast_record = pd.Series({'sstart': 10, 'send': 1, 'qstart': 1, 'qend': 10, 'sstrand': 'plus'})

        self.assertRaises(InvalidPositionException, PointfinderHitHSP, None, blast_record)

    def testBuildPointfinderHitHSPInvalidQueryCoords(self):
        blast_record = pd.Series({'sstart': 1, 'send': 10, 'qstart': 10, 'qend': 1, 'sstrand': 'plus'})

        self.assertRaises(InvalidPositionException, PointfinderHitHSP, None, blast_record)
