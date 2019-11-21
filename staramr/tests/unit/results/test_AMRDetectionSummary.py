import logging
import unittest

import pandas as pd
from os import path

from staramr.results.AMRDetectionSummary import AMRDetectionSummary
from staramr.results.AMRDetectionSummaryResistance import AMRDetectionSummaryResistance

logger = logging.getLogger('AMRDetectionSummaryTest')


class AMRDetectionSummaryTest(unittest.TestCase):

    def setUp(self):
        self.test_data_dir = path.join(path.dirname(__file__), '..', 'data')
        self.columns_resfinder = ('Isolate ID', 'Gene', '%Identity', '%Overlap',
                                  'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession')
        self.columns_pointfinder = ('Isolate ID', 'Gene', 'Type', 'Position', 'Mutation',
                                    '%Identity', '%Overlap', 'HSP Length/Total Length')
        self.columns_plasmidfinder = ('Isolate ID', 'Gene', '%Identity', '%Overlap',
                                      'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession')
        self.columns_quality_module = ('Isolate ID','Genome Length','N50 value','Number of Contigs Under 1000 bp','Quality Module','Quality Module Feedback')
        self.detailed_summary = ('Isolate ID', 'Gene', 'Predicted Phenotype', '%Identity', '%Overlap',
                                 'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession', 'Data Type')

        # Resfinder tables
        self.resfinder_table_empty = pd.DataFrame([],
                                                  columns=self.columns_resfinder)

        self.resfinder_table1 = pd.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table1_files = ['file1']

        self.resfinder_table_mult_resistance = pd.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_resistance_files = ['file1']

        self.resfinder_table_mult_gene_same_resistance = pd.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_gene_same_resistance_files = ['file1']

        self.resfinder_table_mult_same_gene_same_resistance = pd.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_same_gene_same_resistance_files = ['file1']

        self.resfinder_table_mult_file = pd.DataFrame([
            ['file1', 'blaIMP-42', 99.73, 100.00, '741/741', 'blaIMP-42_1_AB753456', 1, 741,
             'AB753456'],
            ['file1', 'newGene', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
            ['file2', 'blaIMP-42', 99.73, 100.00, '741/741', 'newGene', 1, 741,
             'AB753456'],
        ],
            columns=self.columns_resfinder)
        self.resfinder_table_mult_file_files = ['file1', 'file2']

        # Plasmidfinder tables
        self.plasmidfinder_table_empty = pd.DataFrame([],
                                                      columns=self.columns_plasmidfinder)

        self.plasmidfinder_table1 = pd.DataFrame([
            ['file4', 'IncFIB(S)', 100.00, 100.00, '643/643', 'ref|NC_003277.2|', 17653, 17011,
             'FN432031'],
        ],
            columns=self.columns_plasmidfinder)
        self.plasmidfinder_table1_files = ['file4']

        self.plasmidfinder_table_mult_resistance = pd.DataFrame([
            ['file4', 'IncFIB(S)', 100.00, 100.00, '643/643', 'ref|NC_003277.2|', 17653, 17011,
             'FN432031'],
            ['file4', 'IncFII(S)', 100.00, 100.00, '262/262', 'ref|NC_003277.2|', 1665, 1926,
             'CP000858'],
        ],
            columns=self.columns_plasmidfinder)
        self.plasmidfinder_table_mult_resistance_files = ['file4']

        self.plasmidfinder_table_mult_gene_same_resistance_files = ['file4']

        self.plasmidfinder_table_mult_gene_same_resistance = pd.DataFrame([
            ['file4', 'IncFIB(S)', 100.00, 100.00, '643/643', 'ref|NC_003277.2|', 17653, 17011,
             'FN432031'],
            ['file4', 'IncFII(S)', 100.00, 100.00, '262/262', 'ref|NC_003277.2|', 1665, 1926,
             'CP000858'],
        ],
            columns=self.columns_plasmidfinder)

        self.plasmidfinder_table_mult_same_gene_same_resistance_files = ['file4']

        self.plasmidfinder_table_mult_same_gene_same_resistance = pd.DataFrame([
            ['file4', 'IncFIB(S)', 100.00, 100.00, '643/643', 'ref|NC_003277.2|', 17653, 17011,
             'FN432031'],
            ['file4', 'IncFII(S)', 100.00, 100.00, '262/262', 'ref|NC_003277.2|', 1665, 1926,
             'CP000858'],
        ],
            columns=self.columns_plasmidfinder)

        self.plasmidfinder_table_mult_file_files = ['file4', 'file5']

        self.plasmidfinder_table_mult_file = pd.DataFrame([
            ['file4', 'IncFIB(S)', 100.00, 100.00, '643/643', 'ref|NC_003277.2|', 17653, 17011,
             'FN432031'],
            ['file4', 'IncFII(S)', 100.00, 100.00, '262/262', 'ref|NC_003277.2|', 1665, 1926,
             'CP000858'],
            ['file5', 'IncFIB(K)', 98.93, 100.00, '560/560', 'ref|NC_006856.1|', 118238, 117679,
             'JN233704'],
        ],
            columns=self.columns_plasmidfinder)

        # Pointfinder tables
        self.pointfinder_table = pd.DataFrame([
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        self.pointfinder_table_files = ['file1']

        self.pointfinder_table_multiple_gene = pd.DataFrame([
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
            ['file1', 'gyrAB', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        self.pointfinder_table_multiple_gene_files = ['file1']

        self.pointfinder_table_empty = pd.DataFrame([],
                                                    columns=self.columns_pointfinder)

        self.quality_module_table_single_file = pd.DataFrame([['file1',6000000,11000,0,'Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')
                
        self.quality_module_table_two_files = pd.DataFrame([['file1',6000000,11000,0,'Pass',''],
                                    ['file2','6000000','11000','0','Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')

        self.quality_module_table_three_files = pd.DataFrame([['file1',6000000,11000,0,'Pass',''],
                                    ['file2','6000000','11000','0','Pass',''],
                                    ['file3','6000000','11000','0','Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')

        self.quality_module_table_file_4 = pd.DataFrame([
                                    ['file4',6000000,11000,0,'Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')

        self.quality_module_table_file_4_and_5 = pd.DataFrame([['file4','6000000','11000','0','Pass',''],
                                    ['file5',6000000,11000,0,'Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')

        self.quality_module_table_file_4_and_5_and_6 = pd.DataFrame([['file4','6000000','11000','0','Pass',''],
                                    ['file5',6000000,11000,0,'Pass',''],
                                    ['file6',6000000,11000,0,'Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')



        # Detailed Summary Tables

        self.quality_module_table_file_1 = pd.DataFrame([['file1',6000000,11000,0,'Pass',''],
                                    ['file1','6000000','11000','0','Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')    

        self.quality_module_table_file_1_and_2_and_4_and_5= pd.DataFrame([['file1',6000000,11000,0,'Pass',''],
                                    ['file2',6000000,11000,0,'Pass',''],
                                    ['file4',6000000,11000,0,'Pass',''],
                                    ['file5',6000000,11000,0,'Pass',''],],
                                    columns=self.columns_quality_module).set_index('Isolate ID')      

        self.plasmidfinder_table_None = None

        self.detailed_summary_multi_files = ['file1', 'file2', 'file4', 'file5']

        

    def testResfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table1_files, self.resfinder_table1, self.quality_module_table_single_file,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()
        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')
        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderSingleGeneWithNegativesNonIncluded(self):
        files = ['file1', 'file2']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.quality_module_table_single_file, self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderSingleGeneWithNegativesIncluded(self):
        files = ['file1', 'file2']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.quality_module_table_two_files, self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file2', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')

    def testResfinderSingleGeneWithTwoNegativesIncluded(self):
        files = ['file1', 'file2', 'file3']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.quality_module_table_three_files, self.plasmidfinder_table_empty)

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
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.quality_module_table_three_files, self.plasmidfinder_table_empty)

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
                                                    self.resfinder_table_mult_resistance,
                                                    self.quality_module_table_single_file,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderMultipleGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_gene_same_resistance_files,
                                                    self.resfinder_table_mult_gene_same_resistance,
                                                    self.quality_module_table_single_file,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderMultipleSameGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_same_gene_same_resistance_files,
                                                    self.resfinder_table_mult_same_gene_same_resistance,
                                                    self.quality_module_table_single_file,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, blaIMP-42', summary['Genotype'].iloc[0], 'Genes not equal')

    def testResfinderMultipleFile(self):
        amr_detection_summary = AMRDetectionSummary(self.resfinder_table_mult_file_files,
                                                    self.resfinder_table_mult_file, self.quality_module_table_two_files, self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, newGene', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file2', summary.index[1], 'File name not equal')
        self.assertEqual('blaIMP-42', summary['Genotype'].iloc[1], 'Genes not equal')

    def testPlasmidfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.plasmidfinder_table1_files, self.resfinder_table_empty,self.quality_module_table_file_4,
                                                    self.plasmidfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')
        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S)', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPlasmidfinderSingleGeneWithNegativesNonIncluded(self):
        files = ['file4', 'file5']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, self.quality_module_table_file_4_and_5, self.plasmidfinder_table1)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S)', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPlasmidfinderSingleGeneWithNegativesIncluded(self):
        files = ['file4', 'file5']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, self.quality_module_table_file_4_and_5, self.plasmidfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S)', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file5', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')

    def testPlasmidfinderSingleGeneWithTwoNegativesIncluded(self):
        files = ['file4', 'file5', 'file6']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, self.quality_module_table_file_4_and_5_and_6, self.plasmidfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(3, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S)', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file5', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')
        self.assertEqual('file6', summary.index[2], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[2], 'Negative gene not valid')

    def testPlasmidfinderSingleGeneWithTwoNegativesIncludedDifferentOrder(self):
        files = ['file5', 'file4', 'file6']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, self.quality_module_table_file_4_and_5_and_6, self.plasmidfinder_table1)

        summary = amr_detection_summary.create_summary(True)

        self.assertEqual(3, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S)', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file5', summary.index[1], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[1], 'Negative gene not valid')
        self.assertEqual('file6', summary.index[2], 'Negative file not included')
        self.assertEqual('None', summary['Genotype'].iloc[2], 'Negative gene not valid')

    def testPlasmidfinderMultipleGeneResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.plasmidfinder_table_mult_resistance_files,
                                                    self.resfinder_table_empty,
                                                    self.quality_module_table_file_4,
                                                    self.plasmidfinder_table_mult_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S), IncFII(S)', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPlasmidfinderMultipleGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.plasmidfinder_table_mult_gene_same_resistance_files,
                                                    self.resfinder_table_empty,
                                                    self.quality_module_table_file_4,
                                                    self.plasmidfinder_table_mult_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S), IncFII(S)', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPlasmidfinderMultipleSameGeneSameResistance(self):
        amr_detection_summary = AMRDetectionSummary(self.plasmidfinder_table_mult_same_gene_same_resistance_files,
                                                    self.resfinder_table_empty,
                                                    self.quality_module_table_file_4,
                                                    self.plasmidfinder_table_mult_same_gene_same_resistance)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S), IncFII(S)', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPlasmidfinderMultipleFile(self):
        amr_detection_summary = AMRDetectionSummary(self.plasmidfinder_table_mult_file_files,
                                                    self.resfinder_table_empty, self.quality_module_table_file_4_and_5, self.plasmidfinder_table_mult_file)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(2, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file4', summary.index[0], 'File name not equal')
        self.assertEqual('IncFIB(S), IncFII(S)', summary['Genotype'].iloc[0], 'Genes not equal')
        self.assertEqual('file5', summary.index[1], 'File name not equal')
        self.assertEqual('IncFIB(K)', summary['Genotype'].iloc[1], 'Genes not equal')

    def testPointfinderSingleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.pointfinder_table_files, self.resfinder_table_empty,self.quality_module_table_single_file,
                                                    self.pointfinder_table, self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPointfinderSingleMultipleGene(self):
        amr_detection_summary = AMRDetectionSummary(self.pointfinder_table_multiple_gene_files,
                                                    self.resfinder_table_empty, self.quality_module_table_single_file, self.pointfinder_table_multiple_gene,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA, gyrAB', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPointfinderSingleMultipleGeneSame(self):
        df = pd.DataFrame([
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
            ['file1', 'gyrA', 'codon', 67, 'GCC -> CCC (A -> P)', 99.96, 100.0, '2637/2637'],
        ],
            columns=self.columns_pointfinder)
        files = ['file1']

        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table_empty, self.quality_module_table_single_file, df,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('gyrA, gyrA', summary['Genotype'].iloc[0], 'Genes not equal')

    def testPointfinderResfinderSingleGene(self):
        files = ['file1']
        amr_detection_summary = AMRDetectionSummary(files, self.resfinder_table1, self.quality_module_table_single_file, self.pointfinder_table,
                                                    self.plasmidfinder_table_empty)

        summary = amr_detection_summary.create_summary()

        self.assertEqual(1, len(summary.index), 'Invalid number of rows in results')

        self.assertEqual('file1', summary.index[0], 'File name not equal')
        self.assertEqual('blaIMP-42, gyrA', summary['Genotype'].iloc[0], 'Genes not equal')

    def testDetailedSummary_noPlasmid_noPoint(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table1.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(2, len(detailed_summary.index), 'Invalid number of rows, expected 2')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')

        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('blaIMP-42', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')

    def testDetailedSummary_noPoint(self):
        plasmid_table = self.plasmidfinder_table1
        plasmid_table['Isolate ID'] = 'file1'
        plasmid_table = plasmid_table.set_index('Isolate ID')

        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table1.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              plasmidfinder_dataframe=plasmid_table)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(2, len(detailed_summary.index), 'Invalid number of rows, expected 2')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('IncFIB(S)', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('blaIMP-42', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')

    def testDetailedSummary_noPlasmid(self):
        point_table = self.pointfinder_table
        point_table = point_table.set_index('Isolate ID')

        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table1.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              pointfinder_dataframe=point_table)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(3, len(detailed_summary.index), 'Invalid number of rows, expected 3')

        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[2], 'File name not equal')

        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')
        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[2], 'Incorrect Data Type')

        self.assertEqual('blaIMP-42', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('gyrA', detailed_summary['Gene'].iloc[1], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[2], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[2], 'Predicted Phenotype not equal')

    def testDetailedSummary_noRes(self):
        point_table = self.pointfinder_table
        point_table = point_table.set_index('Isolate ID')

        plasmid_table = self.plasmidfinder_table1
        plasmid_table['Isolate ID'] = 'file1'
        plasmid_table = plasmid_table.set_index('Isolate ID')

        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table_empty.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              pointfinder_dataframe=point_table,
                                                              plasmidfinder_dataframe=plasmid_table)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(2, len(detailed_summary.index), 'Invalid number of rows, expected 2')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('IncFIB(S)', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('gyrA', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')

    def testDetailedSummary_noRes_noPoint(self):
        plasmid_table = self.plasmidfinder_table1
        plasmid_table['Isolate ID'] = 'file1'
        plasmid_table = plasmid_table.set_index('Isolate ID')

        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table_empty.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              plasmidfinder_dataframe=plasmid_table)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(2, len(detailed_summary.index), 'Invalid number of rows, expected 2')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('IncFIB(S)', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('Sensitive', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')

    def testDetailedSummary_noRes_noPlasmid(self):
        point_table = self.pointfinder_table
        point_table = point_table.set_index('Isolate ID')

        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table_empty.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              pointfinder_dataframe=point_table)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(2, len(detailed_summary.index), 'Invalid number of rows, expected 2')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')

        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('gyrA', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')

    def testDetailedSummary_allFinders(self):
        point_table = self.pointfinder_table
        point_table = point_table.set_index('Isolate ID')

        plasmid_table = self.plasmidfinder_table1
        plasmid_table['Isolate ID'] = 'file1'
        plasmid_table = plasmid_table.set_index('Isolate ID')

        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table1.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              pointfinder_dataframe=point_table,
                                                              plasmidfinder_dataframe=plasmid_table)

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(3, len(detailed_summary.index), 'Invalid number of rows, expected 3')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[2], 'File name not equal')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[2], 'Incorrect Data Type')

        self.assertEqual('IncFIB(S)', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('blaIMP-42', detailed_summary['Gene'].iloc[1], 'Genes not equal')
        self.assertEqual('gyrA', detailed_summary['Gene'].iloc[2], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')
        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[2], 'Predicted Phenotype not equal')

    def testDetailedSummary_noFinders(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.resfinder_table1_files,
                                                              self.resfinder_table_empty.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1,
                                                              self.pointfinder_table_empty.set_index('Isolate ID'),
                                                              self.plasmidfinder_table_empty.set_index('Isolate ID'))

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(2, len(detailed_summary.index), 'Invalid number of rows, expected 2')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[1], 'File name not equal')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('None', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('', detailed_summary['Predicted Phenotype'].iloc[0], 'Predicted Phenotype not equal')
        self.assertEqual('Sensitive', detailed_summary['Predicted Phenotype'].iloc[1], 'Predicted Phenotype not equal')

    def testDetailedSummary_multiFiles(self):
        amr_detection_summary = AMRDetectionSummaryResistance(self.detailed_summary_multi_files,
                                                              self.resfinder_table_mult_file.set_index('Isolate ID'),
                                                              self.quality_module_table_file_1_and_2_and_4_and_5,
                                                              self.pointfinder_table_multiple_gene.set_index(
                                                                  'Isolate ID'),
                                                              self.plasmidfinder_table_mult_file.set_index(
                                                                  'Isolate ID'))

        detailed_summary = amr_detection_summary.create_detailed_summary()

        self.assertEqual(12, len(detailed_summary.index), 'Invalid number of rows, expected 12')
        self.assertEqual('file1', detailed_summary.index[0], 'File name not equal')
        self.assertEqual('file1', detailed_summary.index[4], 'File name not equal')

        self.assertEqual('file2', detailed_summary.index[5], 'File name not equal')
        self.assertEqual('file2', detailed_summary.index[6], 'File name not equal')

        self.assertEqual('file4', detailed_summary.index[7], 'File name not equal')
        self.assertEqual('file4', detailed_summary.index[9], 'File name not equal')

        self.assertEqual('file5', detailed_summary.index[10], 'File name not equal')
        self.assertEqual('file5', detailed_summary.index[11], 'File name not equal')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[0], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[1], 'Incorrect Data Type')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[5], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[6], 'Incorrect Data Type')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[8], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[9], 'Incorrect Data Type')

        self.assertEqual('Plasmid', detailed_summary['Data Type'].iloc[10], 'Incorrect Data Type')
        self.assertEqual('Resistance', detailed_summary['Data Type'].iloc[11], 'Incorrect Data Type')

        self.assertEqual('None', detailed_summary['Gene'].iloc[0], 'Genes not equal')
        self.assertEqual('blaIMP-42', detailed_summary['Gene'].iloc[1], 'Genes not equal')

        self.assertEqual('None', detailed_summary['Gene'].iloc[5], 'Genes not equal')
        self.assertEqual('blaIMP-42', detailed_summary['Gene'].iloc[6], 'Genes not equal')

        self.assertEqual('IncFIB(S)', detailed_summary['Gene'].iloc[7], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[9], 'Genes not equal')

        self.assertEqual('IncFIB(K)', detailed_summary['Gene'].iloc[10], 'Genes not equal')
        self.assertEqual('None', detailed_summary['Gene'].iloc[11], 'Genes not equal')

        self.assertEqual('Sensitive', detailed_summary['Predicted Phenotype'].iloc[9], 'Predicted Phenotype not equal')
        self.assertEqual('Sensitive', detailed_summary['Predicted Phenotype'].iloc[11], 'Predicted Phenotype not equal')
