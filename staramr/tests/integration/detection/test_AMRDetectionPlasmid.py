import logging
import tempfile
import unittest
from os import path

from Bio import SeqIO

from staramr.blast.JobHandler import JobHandler
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager
from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.databases.resistance.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder
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

    def testPlasmidfinderNameSuccess(self):
        file = path.join(self.test_data_dir, "test-plasmids-seq.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        plasmidfinder_results = self.amr_detection.get_plasmidfinder_results()
        self.assertEqual(len(plasmidfinder_results.index), 1, 'Wrong number of rows in result')

        result = plasmidfinder_results[plasmidfinder_results['Plasmid'] == "IncW"]

        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['Accession'].iloc[0], 'EF633507', msg='Wrong accession')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '243/243', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'plasmidfinder_test-plasmids-seq.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['IncW_1__EF633507'].seq, records['IncW_1__EF633507'].seq,
                         "records don't match")

    def testDetailedSummary_ResPlasmid(self):
        file = path.join(self.test_data_dir, "test-detailed-summary.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        detailed_summary_results = self.amr_detection.get_detailed_summary_results()
        self.assertEqual(len(detailed_summary_results.index), 3, 'Wrong number of rows in result')

        plasmid_type = detailed_summary_results[detailed_summary_results['Data'] == "rep1"]
        self.assertEqual(len(plasmid_type.index), 1, 'Wrong number of results detected')
        self.assertEqual(plasmid_type['Predicted Phenotype'].iloc[0], '', msg='Wrong predicted phenotype')
        self.assertAlmostEqual(plasmid_type['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')
        self.assertAlmostEqual(plasmid_type['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(plasmid_type['Accession'].iloc[0], 'X64695', msg='Wrong accession')
        self.assertEqual(plasmid_type['HSP Length/Total Length'].iloc[0], '1491/1491', msg='Wrong lengths')
        self.assertEqual(plasmid_type['Data Type'].iloc[0], 'Plasmid', msg='Wrong data type')

        res_type = detailed_summary_results[detailed_summary_results['Data'] == "tet(47)"]
        self.assertEqual(len(res_type.index), 1, 'Wrong number of results detected')
        self.assertEqual(res_type['Predicted Phenotype'].iloc[0], 'tetracycline', msg='Wrong predicted phenotype')
        self.assertAlmostEqual(res_type['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')
        self.assertAlmostEqual(res_type['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(res_type['Accession'].iloc[0], 'KR857681', msg='Wrong accession')
        self.assertEqual(res_type['HSP Length/Total Length'].iloc[0], '1248/1248', msg='Wrong lengths')
        self.assertEqual(res_type['Data Type'].iloc[0], 'Resistance', msg='Wrong data type')

    def testResistancePlasmidGenesSummary(self):
        file = path.join(self.test_data_dir, "test-resistance-plasmid.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        summary_results = self.amr_detection.get_summary_results()

        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        result = summary_results[summary_results['Genotype'] == "blaIMP-42"]
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result['Predicted Phenotype'].iloc[0],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         msg='Wrong Predicted Phenotype')
        self.assertEqual(result['Plasmid'].iloc[0], 'IncW', msg='Wrong Plasmid Type')

    def testIndexRangePlasmids(self):
        file = path.join(self.test_data_dir, "test-index-range-plasmid.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        summary_results = self.amr_detection.get_summary_results()

        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows')

        self.assertEqual(summary_results['Genotype'].iloc[0], 'None', msg='Wrong Genotype value')
        self.assertEqual(summary_results['Predicted Phenotype'].iloc[0], 'Sensitive',
                         msg='Wrong Predicted Phenotype value')
        self.assertEqual(summary_results['Plasmid'].iloc[0], 'IncFII(pKPX1)', msg='Wrong Plasmid Type')
