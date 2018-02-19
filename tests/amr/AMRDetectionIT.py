from os import path

import unittest

from amr.AMRDetection import AMRDetection
from amr.blast.BlastHandler import BlastHandler
from amr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from amr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase

class AMRDetectionIT(unittest.TestCase):

    def setUp(self):
        self.resfinder_database_dir = path.join("databases", "resfinder")
        self.pointfinder_database_root_dir = path.join("databases", "pointfinder")

        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_database_dir)
        pointfinder_database = None
        blast_handler = BlastHandler(self.resfinder_database, pointfinder_database, threads=2)

        self.amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        self.test_data_dir = path.join("tests", "data")


    def testResfinderBetaLactam2MutationsSuccess(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['GENE'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result['RESFINDER_PHENOTYPE'].iloc[0], 'Beta-lactam resistance', 'Wrong phenotype')
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
        self.assertEqual(result['RESFINDER_PHENOTYPE'].iloc[0], 'Beta-lactam resistance', 'Wrong phenotype')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 100.00, places=2, msg='Wrong pid')


    def testResfinderBetaLactamDelStartFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 92)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')


    def testPointfinderSalmonellaA67PSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, pointfinder_database, threads=2)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['GENE'] == 'gyrA']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 99.962, places=3, msg='Wrong pid')

if __name__ == '__main__':
    unittest.main()