import logging
import os
import tempfile
import unittest
from os import path

import pandas as pd
from Bio import SeqIO

from staramr.blast.JobHandler import JobHandler
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

        self.columns_resfinder = ('Isolate ID', 'Gene', '%Identity', '%Overlap',
                                  'HSP Length/Total Length', 'Contig', 'Start', 'End', 'Accession')

        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_dir)
        self.resfinder_drug_table = ARGDrugTableResfinder()
        self.pointfinder_drug_table = ARGDrugTablePointfinder()
        self.plasmidfinder_database = PlasmidfinderBlastDatabase(self.plasmidfinder_dir)
        self.pointfinder_database = None
        self.blast_out = tempfile.TemporaryDirectory()
        self.blast_handler = JobHandler(
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
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

    def testResFinderCorrectSeq(self):
        self.amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table,
                                                    self.blast_handler, self.pointfinder_drug_table,
                                                    self.pointfinder_database, output_dir=self.outdir.name,)

        file = path.join(self.test_data_dir, "test-resfinder-correct-seq.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(results.index), 1, 'Wrong number of rows in result')
        self.assertEqual(results['Sequence'].iloc[0], "ATGAGCAAGTTATCTGCATTCTTTATATTTTTGTTTTGCAGCATTGATACCGCAGCAGAGTCTTTGCCAGATTTAAAAATTGAAAAGCTTGATGAAGGCGTTTATGTTCATACTTCGTTTGAAGAAGTTAACAGGTGGGGCGTTGTTCCTAAACATGGTTTGGTGGTTCTTGTAAATGCTGAGGCTTACCTAATTGACACTCCATTTACGGCTAAAGATACTGAAAAGTTAGTCACTTGGTTTGTGGAGCGTGGCTATAAAATAAAAGGCAGCATTTCCTCTCATTTTCATAGCGACAGCACGGGCGGAATAGAGTGGCTTAATTCTCGATCTATCCCCACGTATGCATCTGAATTAACAAATGAACTGCTTAAAAAAGACGGTAAGGTTCAAGCCACAAATTCATTTAGCGGAGTTAACTATTGGCTAGTTAAAAATAAAATTGAAGTTTTTTATCCAGGCCCGGGACACACTCCAGATAACGTAGTGGTTTGGTTGCCTGAAAGGAAAATATTATTCGGTGGTTGTTTTATTAAACCGTACGGTTTAGGCAATTTGGGTGACGCAAATATAGAAGCTTGGCCAAAGTCCGCCAAATTATTAAAGTCCAAATATGGTAAGGCAAAACTGGTTGTTCCAAGTCACAGTGAAGTTGGAGACGCATCACTCTTGAAACTTACATTAGAGCAGGCGGTTAAAGGGTTAAACGAAAGTAAAAAACCATCAAAACCAAGCAACTAA", "Incorrect sequence")

    def testNumericalSequenceID(self):
        file = path.join(self.test_data_dir, "test-seq-id.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == "aac(6')-Iaa"]
        self.assertEqual(result['Contig'].iloc[0], "1", "Incorrect contig id")

    def testResfinderBetaLactam2MutationsSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[blaIMP-42_1_AB753456]', 'Wrong phenotype')

    def testResfinderBetaLactam2MutationsFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")]
        self.amr_detection.run_amr_detection(files, 99.8, 90, 90, 90,0,0,0,0,0)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testResfinderBetaLactamDelStartSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 92, 90, 90,0,0,0,0,0)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testResfinderBetaLactamInsStartSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-ins-start.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 91, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 97, 99, 99, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        self.amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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

    def testPointfinderSalmonellaS83ISuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-S83I.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (S83I)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-S83I', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 83, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'TCC -> ATC (S -> I)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.92, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin I/R, nalidixic acid',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-S83I.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonellaA67PSuccessNoPhenotype(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-del-end.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 98, 90,0,0,0,0,0)

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-del-end.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 99, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67PFailPID(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99.97, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67TFail(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        files = [path.join(self.test_data_dir, "gyrA-A67T.fsa")]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67PReverseComplementSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-rc.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        #Make sure that the quality module fails when we pass 0,0,0,0,0 for the quality checking parameters
        summary_results = amr_detection.get_summary_results()
        self.assertEqual('Failed', summary_results['Quality Module'].iloc[0], 'Quality result not equal')

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S_rrsD-1T1065.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,1000,1000000,1000,300,1000)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        #Make sure that the quality module passes when we pass 1000,1000000,1000,300,1000 for the quality checking parameters
        summary_results = amr_detection.get_summary_results()
        self.assertEqual('Passed', summary_results['Quality Module'].iloc[0], 'Quality result not equal')

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S_gyrA_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name, genes_to_exclude=['gyrA'])

        file = path.join(self.test_data_dir, "16S_gyrA_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S-rc_gyrA-rc_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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
        amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

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
        amr_detection.run_amr_detection(files, 99, 90, 90, 90,0,0,0,0,0)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        result_sensitive = summary_results.loc['non-match']
        self.assertTrue(isinstance(result_sensitive, pd.Series), 'Wrong number of results detected')
        self.assertEqual(result_sensitive['Genotype'], 'None', 'Wrong genotype')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderCampylobacterA70TSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'campylobacter')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A70T.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ciprofloxacin, nalidixic acid', 'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A70T.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderCampylobacterA2075GSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'campylobacter')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "23S-A2075G.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

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

    """
    def testPointfinderEFaecalisS97NSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'enterococcus_faecalis')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
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

    def testPointfinderEFaeciumS83RSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'enterococcus_faecium')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-S83R.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (S83R)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-S83R', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 83, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'AGT -> CGT (S -> R)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.96, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2469/2469', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[gyrA (S83R)]', 'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-S83R.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderEFaeciumIns466Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'enterococcus_faecium')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "pbp5-ins-466.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()

        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'pbp5 (ins466S)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'pbp5-ins-466', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 466, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], '--T -> AGT (ins -> S)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.85, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.15, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2040/2037', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[pbp5 (ins466S)]', 'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_pbp5-ins-466.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['pbp5'].seq.upper(), records['pbp5'].seq.upper(), "records don't match")
    """

    def testPointfinderEscherichiaColiCn42TSuccess(self):
        # C-42T (negative coordinate)
        # This tests the ability to find a single point mutation in the nucleotide part of a promoter.
        # Verified CGE identifies the ampC-promoter (g.-42C>T) mutation.
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'escherichia_coli')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "ampC_promoter_size_53bp-Cn42T.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'ampC_promoter_size_53bp (C-42T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'ampC_promoter_size_53bp-Cn42T', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], -42, msg='Wrong nucleotide position')
        self.assertEqual(result['Mutation'].iloc[0], 'C -> T', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.31, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '145/145', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ampicillin, amoxicillin/clavulanic acid, cefoxitin',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_ampC_promoter_size_53bp-Cn42T.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['ampC_promoter_size_53bp'].seq.upper(), records['ampC_promoter_size_53bp'].seq.upper(), "records don't match")

    def testPointfinderMycobacteriumTuberculosisF10ISuccess(self):
        # F10I
        # This tests the ability to find a single mutation in the codon part of a promoter.
        # Verified CGE finds the ahpC-promoter (p.F10I) mutation.
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'mycobacterium_tuberculosis')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "ahpC_promoter_size_180bp-F10I.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'ahpC_promoter_size_180bp (F10I)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'ahpC_promoter_size_180bp-F10I', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 10, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'TTC -> ATC (F -> I)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.87, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '768/768', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[ahpC_promoter_size_180bp (F10I)]',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_ahpC_promoter_size_180bp-F10I.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['ahpC_promoter_size_180bp'].seq.upper(), records['ahpC_promoter_size_180bp'].seq.upper(), "records don't match")

    def testPointfinderEscherichiaColin13insGSuccess(self):
        # Insert a G at -13 in the promoter.
        # (ampC_promoter_size_53bp / ampC promoter / -13 / - / ins / G)
        # This tests the ability to find a single insertion in the nucleotide part of a promoter.
        # Verified CGE identifies the ampC-promoter (g.-13_-12insG) mutation.
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'escherichia_coli')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "ampC_promoter_size_53bp-n13insG.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'ampC_promoter_size_53bp (ins-13G)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'ampC_promoter_size_53bp-n13insG', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], -13, msg='Wrong nucleotide position')
        self.assertEqual(result['Mutation'].iloc[0], 'ins -> G', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.31, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.69, places=2, msg='Wrong overlap')  # Blast reports 100.69
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '146/145', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'ampicillin, amoxicillin/clavulanic acid, cefoxitin',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_ampC_promoter_size_53bp-n13insG.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['ampC_promoter_size_53bp'].seq.upper(), records['ampC_promoter_size_53bp'].seq.upper(), "records don't match")

    def testPointfinderEscherichiaColin16insGTSuccess(self):
        # Insert a GT at -16 in the promoter.
        # (ampC_promoter_size_53bp / ampC promoter / -16 / - / ins / GT)
        # This tests the ability to find a double insertion in the nucleotide part of a promoter.
        # Verified CGE identifies the ampC-promoter (g.-16_-15insGT) mutation.
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'escherichia_coli')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "ampC_promoter_size_53bp-n16insGT.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 90, 90, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'ampC_promoter_size_53bp (ins-16GT)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'ampC_promoter_size_53bp-n16insGT', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], -16, msg='Wrong nucleotide position')
        self.assertEqual(result['Mutation'].iloc[0], 'ins -> GT', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 98.639, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 101.38, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '147/145', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[ampC_promoter_size_53bp (ins-16GT)]',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_ampC_promoter_size_53bp-n16insGT.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['ampC_promoter_size_53bp'].seq.upper(), records['ampC_promoter_size_53bp'].seq.upper(), "records don't match")

    def testPointfinderNeisseriaGonorrhoeaen57delSuccess(self):
        # Delete the A nucleotide at -57 in the promoter.
        # (mtrR_promoter_size_66bp / mtrR promoter / -57 / del / A)
        # This tests the ability to find a single deletion in the nucleotide part of a promoter.
        # Verified CGE identifies the mtrR-promoter (g.-57_-57del) mutation.
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'neisseria_gonorrhoeae')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "mtrR_promoter_size_66bp-n57del.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 90, 90, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'mtrR_promoter_size_66bp (del-57A)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'mtrR_promoter_size_66bp-n57del', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], -57, msg='Wrong nucleotide position')
        self.assertEqual(result['Mutation'].iloc[0], 'del -> A', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.86, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.0, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '699/699', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[mtrR_promoter_size_66bp (del-57A)]',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_mtrR_promoter_size_66bp-n57del.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['mtrR_promoter_size_66bp'].seq.upper(), records['mtrR_promoter_size_66bp'].seq.upper().replace('-', ''), "records don't match")

    def testPointfinderNeisseriaGonorrhoeae78delSuccess(self):
        # Delete the GTT/V nucleotides/codon at pos 78 (nucleotide coords) / pos 27 (codon coords)
        # Reminder that the Pointfinder database uses 0-base for nucleotides and 1-base for codons.
        # DB entry:
        # rpsE / rpsE / 78 / - / del / GTT / Spectinomycin / 23183436
        # This tests the ability to find a single codon deletion.
        # Verified CGE identifies the rpsE codon deletion (rpsE:p.V27_None79del): 
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'neisseria_gonorrhoeae')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "rpsE-del78gtt.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 90, 90, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'rpsE (del27V)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'rpsE-del78gtt', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 27, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GTT -> --- (V -> del)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.42, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.0, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '519/519', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[rpsE (del27V)]',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_rpsE-del78gtt.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['rpsE'].seq.upper(), records['rpsE'].seq.upper().replace('-', ''), "records don't match")

    def testPointfinderEnterococcusFaecium466insSSuccess(self):
        # Insert a AGC (S) at position 466 (codon coords).
        # DB entry:
        # pbp5 / pbp5 / 466 / - / ins / S,D / Ampicillin / 25182648,25182648
        # This tests the ability to find a single codon insertion.
        # Verified CGE identifies the codon insertion pbp5 (ins466_None467insS). 
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'enterococcus_faecium')
        blast_handler = JobHandler({'resfinder': self.resfinder_database, 'pointfinder': pointfinder_database}, 2,
                                     self.blast_out.name)
        amr_detection = AMRDetectionResistance(self.resfinder_database, self.resfinder_drug_table, blast_handler,
                                               self.pointfinder_drug_table, pointfinder_database,
                                               output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "pbp5-ins465S.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 90, 90, 90, 90,0,0,0,0,0)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'pbp5 (ins465S)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'pbp5-ins465S', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 465, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], '--- -> AGC (ins -> S)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.85, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.15, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2040/2037', msg='Wrong lengths')
        self.assertEqual(result['Predicted Phenotype'].iloc[0], 'unknown[pbp5 (ins465S)]',
                         'Wrong phenotype')

        hit_file = path.join(self.outdir.name, 'pointfinder_pbp5-ins465S.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['pbp5'].seq.upper(), records['pbp5'].seq.upper().replace('-', ''), "records don't match")


if __name__ == '__main__':
    unittest.main()
