import unittest

import pandas

from staramr.results.AMRDetectionSummary import AMRDetectionSummary

class AMRDetectionSummaryTest(unittest.TestCase):

    def setUp(self):
        self.columns_resfinder = ('FILE', 'GENE', 'RESFINDER_PHENOTYPE', '%IDENTITY', '%OVERLAP',
                                  'DB_SEQ_LENGTH/QUERY_HSP', 'CONTIG', 'START', 'END', 'ACCESSION')
        self.columns_pointfinder = ('FILE', 'GENE', 'RESFINDER_PHENOTYPE', 'CODON_POSITION', 'NUCLEOTIDE', 'AMINO_ACID',
                                    '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP')

        # Resfinder tables
        self.resfinder_table_empty = pandas.DataFrame([],
                                                      columns=self.columns_resfinder)

        self.resfinder_table1 = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table1_files = ['file1']

        self.resfinder_table_mult_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 'New resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_resistance_files = ['file1']

        self.resfinder_table_mult_gene_same_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_gene_same_resistance_files = ['file1']

        self.resfinder_table_mult_same_gene_same_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_same_gene_same_resistance_files = ['file1']

        self.resfinder_table_mult_file = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 'New resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
            ['file2', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_file_files = ['file1', 'file2']

        # Pointfinder tables
        self.pointfinder_table = pandas.DataFrame([
            ['file1', 'gyrA', 'pfResistance', 67, 'GCC -> CCC', 'A -> P', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        self.pointfinder_table_files = ['file1']

        self.pointfinder_table_multiple_gene = pandas.DataFrame([
            ['file1', 'gyrA', 'pfResistance', 67, 'GCC -> CCC', 'A -> P', 99.96, 100.0, '2637/2637'],
            ['file1', 'gyrAB', 'pfResistance2', 67, 'GCC -> CCC', 'A -> P', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        self.pointfinder_table_multiple_gene_files = ['file1']

    def testResfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table1_files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')

    def testResfinderSingleGeneWithNegativesNonIncluded(self):
        files = ['file1', 'file2']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')

    def testResfinderSingleGeneWithNegativesIncluded(self):
        files = ['file1', 'file2']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['GENE'].iloc[1], 'Negative gene not valid')
        self.assertEqual('Sensitive', summary['RESFINDER_PHENOTYPE'].iloc[1],
                         'Resfinder negative phenotype not equal')

    def testResfinderSingleGeneWithTwoNegativesIncluded(self):
        files = ['file1', 'file2', 'file3']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(3, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['GENE'].iloc[1], 'Negative gene not valid')
        self.assertEqual('Sensitive', summary['RESFINDER_PHENOTYPE'].iloc[1],
                         'Resfinder negative phenotype not equal')
        self.assertEqual('file3', summary.index[2], 'Negative file not included')
        self.assertEqual('None', summary['GENE'].iloc[2], 'Negative gene not valid')
        self.assertEqual('Sensitive', summary['RESFINDER_PHENOTYPE'].iloc[2],
                         'Resfinder negative phenotype not equal')

    def testResfinderSingleGeneWithTwoNegativesIncludedDifferentOrder(self):
        files = ['file2', 'file1', 'file3']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(3, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['GENE'].iloc[1], 'Negative gene not valid')
        self.assertEqual('Sensitive', summary['RESFINDER_PHENOTYPE'].iloc[1],
                         'Resfinder negative phenotype not equal')
        self.assertEqual('file3', summary.index[2], 'Negative file not included')
        self.assertEqual('None', summary['GENE'].iloc[2], 'Negative gene not valid')
        self.assertEqual('Sensitive', summary['RESFINDER_PHENOTYPE'].iloc[2],
                         'Resfinder negative phenotype not equal')

    def testResfinderMultipleGeneResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_resistance_files, self.resfinder_table_mult_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, New resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')

    def testResfinderMultipleGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_gene_same_resistance_files, self.resfinder_table_mult_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')

    def testResfinderMultipleSameGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_same_gene_same_resistance_files, self.resfinder_table_mult_same_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')

    def testResfinderMultipleFile(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_file_files, self.resfinder_table_mult_file)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, New resistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Resfinder phenotype not equal')
        self.assertEqual('file2', summary.index[1], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[1], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[1],
                         'Resfinder phenotype not equal')

    def testPointfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.pointfinder_table_files, self.resfinder_table_empty, self.pointfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('pfResistance', summary['RESFINDER_PHENOTYPE'].iloc[0], 'Pointfinder phenotype not equal')

    def testPointfinderSingleMultipleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.pointfinder_table_multiple_gene_files, self.resfinder_table_empty, self.pointfinder_table_multiple_gene)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA, gyrAB', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('pfResistance, pfResistance2', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Pointfinder phenotype not equal')

    def testPointfinderSingleMultipleGeneSame(self):
        df = pandas.DataFrame([
            ['file1', 'gyrA', 'pfResistance', 67, 'GCC -> CCC', 'A -> P', 99.96, 100.0, '2637/2637'],
            ['file1', 'gyrA', 'pfResistance', 67, 'GCC -> CCC', 'A -> P', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        files = ['file1']

        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, df)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA, gyrA', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('pfResistance, pfResistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Pointfinder phenotype not equal')

    def testPointfinderResfinderSingleGene(self):
        files = ['file1']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.pointfinder_table)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, gyrA', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, pfResistance', summary['RESFINDER_PHENOTYPE'].iloc[0],
                         'Pointfinder phenotype not equal')
