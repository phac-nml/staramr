import logging
import tempfile
import unittest
from os import path

from Bio import SeqIO

from staramr.blast.JobHandler import JobHandler
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager
from staramr.databases.resistance.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder
from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.detection.AMRDetectionResistance import AMRDetectionResistance

logger = logging.getLogger('AMRDetectionMLST')


class AMRDetectionMLST(unittest.TestCase):

    def setUp(self):
        blast_databases_repositories = AMRDatabasesManager.create_default_manager().get_database_repos()
        self.resfinder_dir = blast_databases_repositories.get_repo_dir(
            'resfinder')
        self.pointfinder_dir = blast_databases_repositories.get_repo_dir(
            'pointfinder')
        self.plasmidfinder_dir = blast_databases_repositories.get_repo_dir(
            'plasmidfinder')

        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_dir)
        self.resfinder_drug_table = ARGDrugTableResfinder()
        self.pointfinder_drug_table = ARGDrugTablePointfinder()
        self.plasmidfinder_database = PlasmidfinderBlastDatabase(
            self.plasmidfinder_dir)
        self.pointfinder_database = None
        self.blast_out = tempfile.TemporaryDirectory()
        self.blast_handler = JobHandler(
            {'resfinder': self.resfinder_database, 'pointfinder': self.pointfinder_database,
             'plasmidfinder': self.plasmidfinder_database}, 2, self.blast_out.name)

        self.outdir = tempfile.TemporaryDirectory()
        self.amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table,
                                                    self.blast_handler, self.pointfinder_drug_table,
                                                    self.pointfinder_database, output_dir=self.outdir.name)

        self.test_data_dir = path.join(path.dirname(__file__), '..', 'data')

    def tearDown(self):
        self.blast_out.cleanup()
        self.outdir.cleanup()

    def testMLSTResults(self):
        file = path.join(self.test_data_dir, "test-mlst-summary.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        mlst_results = self.amr_detection.get_mlst_results()

        self.assertEqual(len(mlst_results.index), 1, 'Wrong number of results detected')
        self.assertEqual(len(mlst_results.columns), 9, 'Wrong number of columns detected')

        self.assertTrue(mlst_results['Scheme'].iloc[0] in ['senterica', 'senterica_achtman_2'], msg='Wrong Scheme')
        self.assertTrue(mlst_results['Sequence Type'].iloc[0] in ['-', '1'], msg='Wrong Sequence Type')
        self.assertEqual(mlst_results['Locus 1'].iloc[0], 'aroC(1)', msg='Wrong Locus 1 Result')
        self.assertEqual(mlst_results['Locus 2'].iloc[0], 'dnaN(1)', msg='Wrong Locus 2 Result')
        self.assertEqual(mlst_results['Locus 3'].iloc[0], 'hemD(1)', msg='Wrong Locus 3 Result')
        self.assertEqual(mlst_results['Locus 4'].iloc[0], 'hisD(1)', msg='Wrong Locus 4 Result')
        self.assertEqual(mlst_results['Locus 5'].iloc[0], 'purE(1)', msg='Wrong Locus 5 Result')
        self.assertEqual(mlst_results['Locus 6'].iloc[0], 'sucA(1)', msg='Wrong Locus 6 Result')
        self.assertTrue(mlst_results['Locus 7'].iloc[0] in ['thrA(5)', 'thrA(5,783)'], msg='Wrong Locus 7 Result')

    def testNoMLSTResults(self):
        file = path.join(self.test_data_dir, "gyrA-S97N.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        mlst_results = self.amr_detection.get_mlst_results()

        self.assertEqual(len(mlst_results.index), 1, 'Wrong number of results detected')
        self.assertEqual(len(mlst_results.columns), 2, 'Wrong number of columns detected')

        self.assertEqual(mlst_results['Scheme'].iloc[0], '-', msg='Scheme is found, expected none')
        self.assertEqual(mlst_results['Sequence Type'].iloc[0], '-', msg='Sequence Type is found, expected none')
