import logging
import os
import tempfile
import unittest
from os import path

import pandas as pd
from Bio import SeqIO

from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager
from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.databases.resistance.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder
from staramr.detection.AMRDetection import AMRDetection
from staramr.detection.AMRDetectionResistance import AMRDetectionResistance

logger = logging.getLogger('AMRDetectionPlasmid')


class AMRDetectionPlasmid(unittest.TestCase):

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
        self.blast_handler = BlastHandler(
            {'resfinder': self.resfinder_database, 'pointfinder': self.pointfinder_database, 'plasmidfinder': self.plasmidfinder_database}, 2, self.blast_out.name)

        self.outdir = tempfile.TemporaryDirectory()
        self.amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table,
                                                    self.blast_handler, self.pointfinder_drug_table,
                                                    self.pointfinder_database, output_dir=self.outdir.name)

        self.test_data_dir = path.join(path.dirname(__file__), '..', 'data')

    def tearDown(self):
        self.blast_out.cleanup()
        self.outdir.cleanup()

    def testPlasmidfinderNameSuccess(self):
        file = path.join(self.test_data_dir, "test-plasmids-seq.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        plasmidfinder_results = self.amr_detection.get_plasmidfinder_results()
        logger.debug("results is %s", plasmidfinder_results)
        self.assertEqual(len(plasmidfinder_results.index), 1, 'Wrong number of rows in result')

        result = plasmidfinder_results[plasmidfinder_results['Gene'] == "IncW"]
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['Accession'].iloc[0], 'EF633507', msg='Wrong accession')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '243/243', msg='Wrong lengths')
