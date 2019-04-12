import logging
import os
import tempfile
import unittest
from os import path

import pandas as pd
from Bio import SeqIO

from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager
from staramr.databases.resistance.pointfinder.ARGDrugTablePointfinder import ARGDrugTablePointfinder
from staramr.databases.resistance.resfinder.ARGDrugTableResfinder import ARGDrugTableResfinder
from staramr.detection.AMRDetection import AMRDetection
from staramr.detection.AMRDetectionResistance import AMRDetectionResistance

logger = logging.getLogger('AMRDetectionIT')


class AMRDetectionIT(unittest.TestCase):

    def setUp(self):
        blast_databases_repositories = AMRDatabasesManager.create_default_manager().get_database_repos()
        self.resfinder_dir = blast_databases_repositories.get_repo_dir('resfinder')
        self.pointfinder_dir = blast_databases_repositories.get_repo_dir('pointfinder')
        self.plasmidfinder_dir = blast_databases_repositories.get_repo_dir('plasmidfinder')

        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_dir)
        self.resfinder_drug_table = ARGDrugTableResfinder()
        self.pointfinder_drug_table = ARGDrugTablePointfinder()
        self.plasmidfinder_database = PlasmidfinderBlastDatabase(self.plasmidfinder_dir)
        self.pointfinder_database = None
        self.blast_out = tempfile.TemporaryDirectory()
        self.blast_handler = BlastHandler(
            {'resfinder': self.resfinder_database, 'pointfinder': self.pointfinder_database}, 2, self.blast_out.name)

        self.outdir = tempfile.TemporaryDirectory()
        self.amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table,
                                                    self.blast_handler, self.pointfinder_drug_table,
                                                    self.pointfinder_database, output_dir=self.outdir.name)

        self.test_data_dir = path.join(path.dirname(__file__), '..', 'data')
        self.drug_key_resfinder_invalid_file = path.join(self.test_data_dir, 'gene-drug-tables',
                                                         'drug_key_resfinder_invalid.tsv')
        self.drug_key_pointfinder_invalid_file = path.join(self.test_data_dir, 'gene-drug-tables',
                                                           'drug_key_pointfinder_invalid.tsv')

    def tearDown(self):
        self.blast_out.cleanup()
        self.outdir.cleanup()

    def testResfinderAminoglycosideNameSuccess(self):
        file = path.join(self.test_data_dir, "test-aminoglycoside.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == "aac(6')-Iaa"]
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['Accession'].iloc[0], 'NC_003197', msg='Wrong accession')

    def testResfinderExcludeGeneListSuccess(self):
        self.amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table,
                                                    self.blast_handler, self.pointfinder_drug_table,
                                                    self.pointfinder_database, output_dir=self.outdir.name,
                                                    genes_to_exclude=["aac(6')-Iaa_1_NC_003197"])

        file = path.join(self.test_data_dir, "test-aminoglycoside.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

    def testNumericalSequenceID(self):
        file = path.join(self.test_data_dir, "test-seq-id.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == "aac(6')-Iaa"]
        self.assertEqual(result['Contig'].iloc[0], "1", "Incorrect contig id")

    def testResfinderBetaLactam2MutationsSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '741/741', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactam2MutationsSuccessNoPredictedPhenotype(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '741/741', msg='Wrong lengths')
        self.assertFalse('Predicted Phenotype' in result.columns, 'Should not exist phenotype column')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactam2MutationsSuccessNoMatchDrugTable(self):
        resfinder_drug_table = ARGDrugTableResfinder(self.drug_key_resfinder_invalid_file)
        self.amr_detection = AMRDetectionResistance(self.resfinder_database, resfinder_drug_table,
                                                    self.blast_handler, self.pointfinder_drug_table,
                                                    self.pointfinder_database, output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[blaIMP-42_1_AB753456]', 'Wrong phenotype')

    def testResfinderBetaLactam2MutationsFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")]
        self.amr_detection.run_amr_detection(files, 99.8, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testResfinderBetaLactamDelStartSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 91.90, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '681/741', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-del-start.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamDelStartFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 92, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testResfinderBetaLactamInsStartSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-ins-start.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '741/741', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-ins-start.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'beta-lactam-blaIMP-42-mut-2.fsa'), 'fasta'))
        logger.debug("expected_seq=%s", expected_records['blaIMP-42_1_AB753456'].seq)
        logger.debug("actual_seq=%s", records['blaIMP-42_1_AB753456'].seq)
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamDelMiddleSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-middle.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.33, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '741/741', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-del-middle.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'resfinder_beta-lactam-blaIMP-42-del-middle.fsa'), 'fasta'))
        logger.debug("expected_seq=%s", expected_records['blaIMP-42_1_AB753456'].seq)
        logger.debug("actual_seq=%s", records['blaIMP-42_1_AB753456'].seq)
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamDelMiddleReverseComplementSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-middle-rc.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.33, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong percent overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '741/741', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-del-middle-rc.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'resfinder_beta-lactam-blaIMP-42-del-middle.fsa'), 'fasta'))
        logger.debug("expected_seq=%s", expected_records['blaIMP-42_1_AB753456'].seq)
        logger.debug("actual_seq=%s", records['blaIMP-42_1_AB753456'].seq)
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamInsMiddleSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-ins-middle.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 97, 99, 99, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 98.14, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 101.62, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '753/741', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-ins-middle.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'beta-lactam-blaIMP-42-ins-middle.fsa'), 'fasta'))
        logger.debug("expected_seq=%s", expected_records['blaIMP-42_1_AB753456'].seq)
        logger.debug("actual_seq=%s", records['blaIMP-42_1_AB753456'].seq)
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq.upper(), records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamTwoCopies(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2-two-copies.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 2, 'Wrong number of rows in result')

        result = resfinder_results.iloc[0]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'], 100.00, places=2, msg='Wrong percent overlap')
        self.assertEqual(result['HSP Length/Total Length'], '741/741', msg='Wrong lengths')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 61, msg='Wrong start')
        self.assertEqual(result['End'], 801, msg='Wrong end')
        self.assertEqual(result['Predicted Phenotype'],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        result = resfinder_results.iloc[1]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'], 100.00, places=2, msg='Wrong percent overlap')
        self.assertEqual(result['HSP Length/Total Length'], '741/741', msg='Wrong lengths')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 841, msg='Wrong start')
        self.assertEqual(result['End'], 1581, msg='Wrong end')
        self.assertEqual(result['Predicted Phenotype'],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2-two-copies.fsa')
        records = list(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 2, 'Wrong number of hit records')

        expected_file_data = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        expected_records = SeqIO.to_dict(SeqIO.parse(expected_file_data, 'fasta'))
        # verify I get both copies' sequence data
        for record in records:
            self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, record.seq, "records don't match")

    def testResfinderBetaLactamTwoCopiesOneReverseComplement(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2-two-copies-one-rev-complement.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 2, 'Wrong number of rows in result')

        result = resfinder_results.iloc[0]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'], 100.00, places=2, msg='Wrong percent overlap')
        self.assertEqual(result['HSP Length/Total Length'], '741/741', msg='Wrong lengths')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 61, msg='Wrong start')
        self.assertEqual(result['End'], 801, msg='Wrong end')
        self.assertEqual(result['Predicted Phenotype'],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        result = resfinder_results.iloc[1]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'], 100.00, places=2, msg='Wrong percent overlap')
        self.assertEqual(result['HSP Length/Total Length'], '741/741', msg='Wrong lengths')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 1581, msg='Wrong start')
        self.assertEqual(result['End'], 841, msg='Wrong end')
        self.assertEqual(result['Predicted Phenotype'],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name,
                             'resfinder_beta-lactam-blaIMP-42-mut-2-two-copies-one-rev-complement.fsa')
        records = list(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 2, 'Wrong number of hit records')

        expected_file_data = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        expected_records = SeqIO.to_dict(SeqIO.parse(expected_file_data, 'fasta'))
        # verify I get both copies' sequence data
        for record in records:
            self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, record.seq, "records don't match")

    def testPointfinderSalmonellaA67PSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R, nalidixic acid',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonellaA67PSuccessNoPhenotype(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')
        self.assertFalse('Predicted Phenotype' in result.columns, 'Should not exist Predicted Phenotype column')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonellaA67PDelEndSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-del-end.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 98, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P-del-end', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 98.22, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2590/2637', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R, nalidixic acid',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P-del-end.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonellaA67PDelEndFailPlength(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-del-end.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 99, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67PFailPID(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99.97, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67TFail(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        files = [path.join(self.test_data_dir, "gyrA-A67T.fsa")]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67PReverseComplementSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-rc.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P-rc', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R, nalidixic acid',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P-rc.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, "gyrA-A67P.fsa"), 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonella_16S_rrSD_C1065T_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S_rrsD-1T1065.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == '16S_rrsD (C1065T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S_rrsD-1T1065', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 1065, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'C -> T', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.94, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '1544/1544', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'spectinomycin',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_16S_rrsD-1T1065.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['16S_rrsD'].seq.upper(), records['16S_rrsD'].seq.upper(),
                         "records don't match")

    def testResfinderPointfinderSalmonella_16S_C1065T_gyrA_A67_beta_lactam_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S_gyrA_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        resfinder_results = amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')
        self.assertEqual(result['Predicted Phenotype'].iloc[0],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'resfinder_16S_gyrA_beta-lactam.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))
        self.assertEqual(len(records), 1, 'Wrong number of hit records')
        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'beta-lactam-blaIMP-42-mut-2.fsa'), 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq.upper(),
                         records['blaIMP-42_1_AB753456'].seq.upper(), "records don't match")

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 2, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == '16S_rrsD (C1065T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S_gyrA_beta-lactam', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 1065, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'C -> T', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.94, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '1544/1544', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'spectinomycin',
                         'Wrong phenotype')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S_gyrA_beta-lactam', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R, nalidixic acid',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_16S_gyrA_beta-lactam.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))
        self.assertEqual(len(records), 2, 'Wrong number of hit records')
        expected_records1 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, 'gyrA-A67P.fsa'), 'fasta'))
        self.assertEqual(expected_records1['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")
        expected_records2 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, '16S_rrsD-1T1065.fsa'), 'fasta'))
        self.assertEqual(expected_records2['16S_rrsD'].seq.upper(), records['16S_rrsD'].seq.upper(),
                         "records don't match")

    def testResfinderPointfinderSalmonellaExcludeGenesListSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name, genes_to_exclude=['gyrA'])

        file = path.join(self.test_data_dir, "16S_gyrA_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        resfinder_results = amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == '16S_rrsD (C1065T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')

    def testResfinderPointfinderSalmonella_16Src_C1065T_gyrArc_A67_beta_lactam_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S-rc_gyrA-rc_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        resfinder_results = amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')
        self.assertEqual(result['Predicted Phenotype'].iloc[0],
                         'ampicillin, amoxicillin/clavulanic acid, cefoxitin, ceftriaxone, meropenem',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'resfinder_16S-rc_gyrA-rc_beta-lactam.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))
        self.assertEqual(len(records), 1, 'Wrong number of hit records')
        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'beta-lactam-blaIMP-42-mut-2.fsa'), 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq.upper(),
                         records['blaIMP-42_1_AB753456'].seq.upper(), "records don't match")

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 2, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == '16S_rrsD (C1065T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S-rc_gyrA-rc_beta-lactam', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 1065, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'C -> T', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.94, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '1544/1544', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'spectinomycin',
                         'Wrong phenotype')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S-rc_gyrA-rc_beta-lactam', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R, nalidixic acid',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_16S-rc_gyrA-rc_beta-lactam.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))
        self.assertEqual(len(records), 2, 'Wrong number of hit records')
        expected_records1 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, 'gyrA-A67P.fsa'), 'fasta'))
        self.assertEqual(expected_records1['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")
        expected_records2 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, '16S_rrsD-1T1065.fsa'), 'fasta'))
        self.assertEqual(expected_records2['16S_rrsD'].seq.upper(), records['16S_rrsD'].seq.upper(),
                         "records don't match")

    def testResfinderExcludeNonMatches(self):
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, self.blast_handler,
                                               self.pointfinder_drug_table, self.pointfinder_database,
                                               include_negative_results=False, output_dir=self.outdir.name)
        file_beta_lactam = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        file_non_match = path.join(self.test_data_dir, "non-match.fsa")
        files = [file_beta_lactam, file_non_match]
        amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file_beta_lactam, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderIncludeNonMatches(self):
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, self.blast_handler,
                                               self.pointfinder_drug_table, self.pointfinder_database,
                                               include_negative_results=True, output_dir=self.outdir.name)
        file_beta_lactam = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        file_non_match = path.join(self.test_data_dir, "non-match.fsa")
        files = [file_beta_lactam, file_non_match]
        amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 2, 'Wrong number of rows in result')

        result_beta_lactam = summary_results.loc['beta-lactam-blaIMP-42-mut-2']
        self.assertTrue(isinstance(result_beta_lactam, pd.Series), 'Wrong type of results returned')
        self.assertEqual(result_beta_lactam['Genotype'], 'blaIMP-42', 'Wrong genotype')

        result_sensitive = summary_results.loc['non-match']
        self.assertTrue(isinstance(result_sensitive, pd.Series), 'Wrong type of results returned')
        self.assertEqual(result_sensitive['Genotype'], 'None', 'Wrong genotype')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file_beta_lactam, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testNonMatches(self):
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, self.blast_handler,
                                               self.pointfinder_drug_table, self.pointfinder_database,
                                               include_negative_results=True, output_dir=self.outdir.name)
        files = [path.join(self.test_data_dir, "non-match.fsa")]
        amr_detection.run_amr_detection(files, 99, 90, 90, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        result_sensitive = summary_results.loc['non-match']
        self.assertTrue(isinstance(result_sensitive, pd.Series), 'Wrong number of results detected')
        self.assertEqual(result_sensitive['Genotype'], 'None', 'Wrong genotype')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderCampylobacterA70TSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'campylobacter')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A70T.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A70T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A70T', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 70, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> ACC (A -> T)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2592/2592', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R', 'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A70T.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderCampylobacterA2075GSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'campylobacter')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "23S-A2075G.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == '23S (A2075G)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '23S-A2075G', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 2075, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'A -> G', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.97, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2912/2912', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0],
                         'erythromycin, azithromycin, telithromycin, clindamycin', 'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_23S-A2075G.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['23S'].seq.upper(), records['23S'].seq.upper(), "records don't match")

    @unittest.SkipTest  # type: ignore
    def testPointfinderEFaecalisS97NSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'enterococcus_faecalis')
        blast_handler = BlastHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-S97N.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (S97N)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-S97N', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 97, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'AGT -> AAT (S -> N)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2499/2499', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[gyrA (S97N)]', 'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-S97N.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")


if __name__ == '__main__':
    unittest.main()
