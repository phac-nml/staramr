import unittest

import pandas

from amr.results.AMRDetectionSummary import AMRDetectionSummary


class AMRDetectionSummaryTest(unittest.TestCase):

    def setUp(self):
        columns_resfinder = ('FILE', 'GENE', 'RESFINDER_PHENOTYPE', '%IDENTITY', '%OVERLAP',
                     'DB_SEQ_LENGTH/QUERY_HSP', 'CONTIG', 'START', 'END', 'ACCESSION')
        columns_pointfinder = ('FILE', 'GENE', 'CODON_POSITION', 'NUCLEOTIDE', 'AMINO_ACID',
                     '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP')

        self.pointfinder_table = pandas.DataFrame([
            ['file1', 'gyrA', 67, 'GCC -> CCC', 'A -> P', 99.96, 100.0, '2637/2637'],
        ],
            columns=columns_pointfinder)

        self.resfinder_table1 = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741, 'AB753456'],
        ],
            columns=columns_resfinder)

        self.resfinder_table_mult_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 'New resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=columns_resfinder)

        self.resfinder_table_mult_gene_same_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=columns_resfinder)

        self.resfinder_table_mult_same_gene_same_resistance = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=columns_resfinder)

        self.resfinder_table_mult_file = pandas.DataFrame([
            ['file1', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 'New resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
            ['file2', 'blaIMP-42', 'Beta-lactam resistance', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=columns_resfinder)


    def testResfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary['FILE'].iloc[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0], 'Resfinder phenotype not equal')


    def testResfinderMultipleGeneResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary['FILE'].iloc[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, New resistance', summary['RESFINDER_PHENOTYPE'].iloc[0], 'Resfinder phenotype not equal')


    def testResfinderMultipleGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary['FILE'].iloc[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0], 'Resfinder phenotype not equal')


    def testResfinderMultipleSameGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_same_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary['FILE'].iloc[0], 'File name not equal')
        self.assertEqual('blaIMP-42, blaIMP-42', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[0], 'Resfinder phenotype not equal')


    def testResfinderMultipleFile(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_file)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary['FILE'].iloc[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['GENE'].iloc[0], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance, New resistance', summary['RESFINDER_PHENOTYPE'].iloc[0], 'Resfinder phenotype not equal')
        self.assertEqual('file2', summary['FILE'].iloc[1], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['GENE'].iloc[1], 'Genes not equal')
        self.assertEqual('Beta-lactam resistance', summary['RESFINDER_PHENOTYPE'].iloc[1], 'Resfinder phenotype not equal')
