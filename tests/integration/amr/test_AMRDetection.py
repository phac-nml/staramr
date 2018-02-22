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

        self.test_data_dir = path.join("tests", "integration", "data")


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

        result = pointfinder_results[pointfinder_results['GENE'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P.fsa', msg='Wrong file')
        self.assertEqual(result['RESFINDER_PHENOTYPE'].iloc[0], 'Quinolones', msg='Wrong phenotype')
        self.assertEqual(result['CODON_POSITION'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['NUCLEOTIDE'].iloc[0], 'GCC -> CCC', msg='Wrong nucleotide')
        self.assertEqual(result['AMINO_ACID'].iloc[0], 'A -> P', msg='Wrong amino acid')
        self.assertAlmostEqual(result['%IDENTITY'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%OVERLAP'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['DB_SEQ_LENGTH/QUERY_HSP'].iloc[0], '2637/2637', msg='Wrong lengths')


    def testPointfinderSalmonellaA67PFailPID(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, pointfinder_database, threads=2)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99.97, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')


    def testPointfinderSalmonellaA67TFail(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_database_root_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, pointfinder_database, threads=2)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database)

        files = [path.join(self.test_data_dir, "gyrA-A67T.fsa")]
        amr_detection.run_amr_detection(files, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

if __name__ == '__main__':
    unittest.main()