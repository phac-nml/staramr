import unittest

import pandas

from staramr.results.AMRDetectionSummary import AMRDetectionSummary


class AMRDetectionSummaryTest(unittest.TestCase):

    def setUp(self):
        self.columns_resfinder = ('Isolate ID', 'Gene', '%Identity', '%Overlap',
                                  'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession')
        self.columns_pointfinder = ('Isolate ID', 'Gene', 'Type', 'Position', 'Mutation',
                                    '%Identity', '%Overlap', 'HSP Length/Total Length')

        # Resfinder tables
        self.resfinder_table_empty = pandas.DataFrame([],
                                                      columns=self.columns_resfinder)

        self.resfinder_table1 = pandas.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table1_files = ['file1']

        self.resfinder_table_mult_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_resistance_files = ['file1']

        self.resfinder_table_mult_gene_same_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_gene_same_resistance_files = ['file1']

        self.resfinder_table_mult_same_gene_same_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_same_gene_same_resistance_files = ['file1']

        self.resfinder_table_mult_file = pandas.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
            ['file2', 'blaIMP-42', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_file_files = ['file1', 'file2']

        # Pointfinder tables
        self.pointfinder_table = pandas.DataFrame([
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        self.pointfinder_table_files = ['file1']

        self.pointfinder_table_multiple_gene = pandas.DataFrame([
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
            ['file1', 'gyrAB', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        self.pointfinder_table_multiple_gene_files = ['file1']

    def testResfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table1_files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderSingleGeneWithNegativesNonIncluded(self):
        files = ['file1', 'file2']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderSingleGeneWithNegativesIncluded(self):
        files = ['file1', 'file2']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')

    def testResfinderSingleGeneWithTwoNegativesIncluded(self):
        files = ['file1', 'file2', 'file3']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(3, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')
        self.assertEqual('file3', summary.index[2], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[2], 'Negative gene not valid')

    def testResfinderSingleGeneWithTwoNegativesIncludedDifferentOrder(self):
        files = ['file2', 'file1', 'file3']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(3, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')
        self.assertEqual('file3', summary.index[2], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[2], 'Negative gene not valid')

    def testResfinderMultipleGeneResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_resistance_files,
                                                    self.resfinder_table_mult_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderMultipleGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_gene_same_resistance_files,
                                                    self.resfinder_table_mult_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderMultipleSameGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_same_gene_same_resistance_files,
                                                    self.resfinder_table_mult_same_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderMultipleFile(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_file_files,
                                                    self.resfinder_table_mult_file)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file2', summary.index[1], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[1], 'Genes not equal')

    def testPointfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.pointfinder_table_files, self.resfinder_table_empty,
                                                    self.pointfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPointfinderSingleMultipleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.pointfinder_table_multiple_gene_files,
                                                    self.resfinder_table_empty, self.pointfinder_table_multiple_gene)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA, gyrAB', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPointfinderSingleMultipleGeneSame(self):
        df = pandas.DataFrame([
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        files = ['file1']

        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, df)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA, gyrA', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPointfinderResfinderSingleGene(self):
        files = ['file1']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.pointfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, gyrA', summary['Genotype'].iloc[0], 'Genes not equal')
