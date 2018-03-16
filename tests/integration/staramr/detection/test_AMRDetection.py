import unittest
from os import path
import pandas

from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.detection.AMRDetection import AMRDetection

class AMRDetectionIT(unittest.TestCase):

    def setUp(self):
        self.resfinder_database_dir = path.join("databases", "resfinder")
        self.pointfinder_database_root_dir = path.join("databases", "pointfinder")

        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_database_dir)
        self.pointfinder_database = None
        self.blast_handler = BlastHandler(self.resfinder_database, 2, self.pointfinder_database)

        self.amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database)

        self.test_data_dir = path.join("tests", "integration", "data")

    def testResfinderBetaLactam2MutationsSuccess(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['GENE'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 99.73, places=2, msg='Wrong pid')

    def testResfinderBetaLactam2MutationsFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")]
        self.amr_detection.run_amr_detection(files, 99.8, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

    def testResfinderBetaLactamDelStartSuccess(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 91)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['GENE'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 100.00, places=2, msg='Wrong pid')

    def testResfinderBetaLactamDelStartFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 92)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

    def testPointfinderSalmonellaA67PSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['GENE'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P.fsa', msg='Wrong file')
        self.assertEqual(result['TYPE'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['POSITION'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['MUTATION'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%OVERLAP'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['DB_SEQ_LENGTH/QUERY_HSP'].iloc[0], '2637/2637', msg='Wrong lengths')

    def testPointfinderSalmonellaA67PFailPID(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99.97, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

    def testPointfinderSalmonellaA67TFail(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67T.fsa")]
        amr_detection.run_amr_detection(files, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

    def testPointfinderSalmonellaA67PReverseComplementSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67P-rc.fsa")]
        amr_detection.run_amr_detection(files, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['GENE'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P-rc.fsa', msg='Wrong file')
        self.assertEqual(result['TYPE'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['POSITION'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['MUTATION'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%OVERLAP'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['DB_SEQ_LENGTH/QUERY_HSP'].iloc[0], '2637/2637', msg='Wrong lengths')

    def testPointfinderSalmonella_16S_rrSD_C1065T_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "16S_rrsD-1T1065.fsa")]
        amr_detection.run_amr_detection(files, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['GENE'] == '16S_rrsD (C1065T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S_rrsD-1T1065.fsa', msg='Wrong file')
        self.assertEqual(result['TYPE'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['POSITION'].iloc[0], 1065, msg='Wrong codon position')
        self.assertEqual(result['MUTATION'].iloc[0], 'C -> T', msg='Wrong mutation')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 99.935, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%OVERLAP'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['DB_SEQ_LENGTH/QUERY_HSP'].iloc[0], '1544/1544', msg='Wrong lengths')

    def testResfinderExcludeNonMatches(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database, False)
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa"),
                 path.join(self.test_data_dir, "non-match.fsa")]
        amr_detection.run_amr_detection(files, 99, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        summary_results.loc['beta-lactam-blaIMP-42-mut-2.fsa']

    def testResfinderIncludeNonMatches(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database, True)
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa"),
                 path.join(self.test_data_dir, "non-match.fsa")]
        amr_detection.run_amr_detection(files, 99, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 2, 'Wrong number of rows in result')

        result_beta_lactam = summary_results.loc['beta-lactam-blaIMP-42-mut-2.fsa']
        self.assertTrue(isinstance(result_beta_lactam, pandas.Series), 'Wrong type of results returned')
        self.assertEqual(result_beta_lactam['GENE'], 'blaIMP-42', 'Wrong genotype')

        result_sensitive = summary_results.loc['non-match.fsa']
        self.assertTrue(isinstance(result_sensitive, pandas.Series), 'Wrong type of results returned')
        self.assertEqual(result_sensitive['GENE'], 'None', 'Wrong genotype')

    def testNonMatches(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database, True)
        files = [path.join(self.test_data_dir, "non-match.fsa")]
        amr_detection.run_amr_detection(files, 99, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        result_sensitive = summary_results.loc['non-match.fsa']
        self.assertTrue(isinstance(result_sensitive, pandas.Series), 'Wrong number of results detected')
        self.assertEqual(result_sensitive['GENE'], 'None', 'Wrong genotype')


if __name__ == '__main__':
    unittest.main()
