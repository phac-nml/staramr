import unittest

import pandas as pd

from staramr.results.AMRDetectionSummaryResistance import AMRDetectionSummaryResistance


class AMRDetectionSummaryResistanceTest(unittest.TestCase):

    def setUp(self):
        self.columns_resfinder = ('Isolate ID', 'Gene', 'Predicted Phenotype', '%Identity', '%Overlap',
                                  'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession')
        self.columns_pointfinder = ('Isolate ID', 'Gene', 'Predicted Phenotype', 'Type', 'Position', 'Mutation',
                                    '%Identity', '%Overlap', 'HSP Length/Total Length')

        # Resfinder tables
        self.resfinder_table_empty = pd.DataFrame([],
                                                  columns=self.columns_resfinder)

        self.resfinder_table = pd.DataFrame([
            ['file1', 'blaIMP-42', 'ampicillin, amoxi/clav, cefoxitin, ceftriaxone, meropenem', 99.73, 100.00,
             '741/741', 'blaIMP-42_1_AB753456', 1, 741, 'AB753456']
        ],
            columns=self.columns_resfinder)

        self.resfinder_table_duplicate_resistances = pd.DataFrame([
            ['file1', 'blaIMP-42', 'ampicillin', 99.73, 100.00,
             '741/741', 'blaIMP-42_1_AB753456', 1, 741, 'AB753456'],
            ['file1', 'blaCTX-M-55', 'ampicillin, ceftriaxone', 99.73, 100.00,
             '741/741', 'x', 1, 741, 'AB753456']
        ],
            columns=self.columns_resfinder)

        self.pointfinder_table = pd.DataFrame([
            ['file1', 'gyrA', 'ciprofloxacin I/R, nalidixic acid', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0,
             '2637/2637'],
        ],
            columns=self.columns_pointfinder)

        self.pointfinder_table_duplicate = pd.DataFrame([
            ['file1', 'gyrA', 'ampicillin, ceftriaxone, ciprofloxacin I/R', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96,
             100.0,
             '2637/2637'],
        ],
            columns=self.columns_pointfinder)

        self.files = ['file1']

    def testResfinder(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.files, self.resfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('ampicillin, amoxi/clav, cefoxitin, ceftriaxone, meropenem',
                         summary['Predicted Phenotype'].iloc[0], 'Phenotypes not equal')

    def testResfinderDuplicateResistances(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.files, self.resfinder_table_duplicate_resistances)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaCTX-M-55, blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('ampicillin, ceftriaxone',
                         summary['Predicted Phenotype'].iloc[0], 'Phenotypes not equal')

    def testPointfinder(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.files, self.resfinder_table_empty,
                                                              self.pointfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('ciprofloxacin I/R, nalidixic acid',
                         summary['Predicted Phenotype'].iloc[0], 'Phenotypes not equal')

    def testPointfinderResfinder(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.files, self.resfinder_table, self.pointfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, gyrA', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('ampicillin, amoxi/clav, cefoxitin, ceftriaxone, meropenem, ciprofloxacin I/R, nalidixic acid',
                         summary['Predicted Phenotype'].iloc[0], 'Phenotypes not equal')

    def testPointfinderResfinderDuplicate(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.files, self.resfinder_table_duplicate_resistances,
                                                              self.pointfinder_table_duplicate)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaCTX-M-55, blaIMP-42, gyrA', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('ampicillin, ceftriaxone, ciprofloxacin I/R',
                         summary['Predicted Phenotype'].iloc[0], 'Phenotypes not equal')
