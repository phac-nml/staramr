import logging
import os
import tempfile
import unittest
from os import path

import pandas
from Bio import SeqIO

from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.AMRDatabaseHandlerFactory import AMRDatabaseHandlerFactory
from staramr.detection.AMRDetection import AMRDetection

logger = logging.getLogger('AMRDetectionIT')


class AMRDetectionIT(unittest.TestCase):

    def setUp(self):
        database_handler = AMRDatabaseHandlerFactory.create_default_factory().get_database_handler()
        self.resfinder_dir = database_handler.get_resfinder_dir()
        self.pointfinder_dir = database_handler.get_pointfinder_dir()

        self.resfinder_database = ResfinderBlastDatabase(self.resfinder_dir)
        self.pointfinder_database = None
        self.blast_handler = BlastHandler(self.resfinder_database, 2, self.pointfinder_database)

        self.outdir = tempfile.TemporaryDirectory()
        self.amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database,
                                          output_dir=self.outdir.name)

        self.test_data_dir = path.join(path.dirname(__file__), '..', 'data')

    def tearDown(self):
        self.outdir.cleanup()

    def testResfinderBetaLactam2MutationsSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactam2MutationsFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")]
        self.amr_detection.run_amr_detection(files, 99.8, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testResfinderBetaLactamDelStartSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 100.00, places=2, msg='Wrong pid')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-del-start.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamDelStartFail(self):
        files = [path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-start.fsa")]
        self.amr_detection.run_amr_detection(files, 99, 92, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testResfinderBetaLactamDelMiddleSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-middle.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.33, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong percent overlap')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-del-middle.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'resfinder_beta-lactam-blaIMP-42-del-middle.fsa'), 'fasta'))
        logger.debug('expected_seq=' + expected_records['blaIMP-42_1_AB753456'].seq)
        logger.debug('actual_seq=' + records['blaIMP-42_1_AB753456'].seq)
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamDelMiddleReverseComplementSuccess(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-del-middle-rc.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 91, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.33, places=2, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong percent overlap')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-del-middle-rc.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(
            SeqIO.parse(path.join(self.test_data_dir, 'resfinder_beta-lactam-blaIMP-42-del-middle.fsa'), 'fasta'))
        logger.debug('expected_seq=' + expected_records['blaIMP-42_1_AB753456'].seq)
        logger.debug('actual_seq=' + records['blaIMP-42_1_AB753456'].seq)
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderBetaLactamTwoCopies(self):
        file = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2-two-copies.fsa")
        files = [file]
        self.amr_detection.run_amr_detection(files, 99, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 2, 'Wrong number of rows in result')

        result = resfinder_results.iloc[0]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 61, msg='Wrong start')
        self.assertEqual(result['End'], 801, msg='Wrong end')

        result = resfinder_results.iloc[1]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 841, msg='Wrong start')
        self.assertEqual(result['End'], 1581, msg='Wrong end')

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
        self.amr_detection.run_amr_detection(files, 99, 90, 90)

        resfinder_results = self.amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 2, 'Wrong number of rows in result')

        result = resfinder_results.iloc[0]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 61, msg='Wrong start')
        self.assertEqual(result['End'], 801, msg='Wrong end')

        result = resfinder_results.iloc[1]
        self.assertEqual(result['Gene'], 'blaIMP-42', 'Wrong gene for result')
        self.assertAlmostEqual(result['%Identity'], 99.73, places=2, msg='Wrong pid')
        self.assertEqual(result['Contig'], 'blaIMP-42_1_AB753456', msg='Wrong contig name')
        self.assertEqual(result['Start'], 841, msg='Wrong start')
        self.assertEqual(result['End'], 1581, msg='Wrong end')

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
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonellaA67PDelEndSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-del-end.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 98)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P-del-end', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.961, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 98.22, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2590/2637', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P-del-end.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonellaA67PDelEndFailPlength(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-del-end.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 99)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67PFailPID(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        files = [path.join(self.test_data_dir, "gyrA-A67P.fsa")]
        amr_detection.run_amr_detection(files, 99.97, 99, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67TFail(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        files = [path.join(self.test_data_dir, "gyrA-A67T.fsa")]
        amr_detection.run_amr_detection(files, 99, 99, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 0, 'Wrong number of rows in result')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')

    def testPointfinderSalmonellaA67PReverseComplementSuccess(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "gyrA-A67P-rc.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], 'gyrA-A67P-rc', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'pointfinder_gyrA-A67P-rc.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, "gyrA-A67P.fsa"), 'fasta'))
        self.assertEqual(expected_records['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")

    def testPointfinderSalmonella_16S_rrSD_C1065T_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S_rrsD-1T1065.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90)

        pointfinder_results = amr_detection.get_pointfinder_results()
        self.assertEqual(len(pointfinder_results.index), 1, 'Wrong number of rows in result')

        result = pointfinder_results[pointfinder_results['Gene'] == '16S_rrsD (C1065T)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S_rrsD-1T1065', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'nucleotide', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 1065, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'C -> T', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.935, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '1544/1544', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'pointfinder_16S_rrsD-1T1065.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        self.assertEqual(expected_records['16S_rrsD'].seq.upper(), records['16S_rrsD'].seq.upper(),
                         "records don't match")

    def testResfinderPointfinderSalmonella_16S_C1065T_gyrA_A67_beta_lactam_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S_gyrA_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90)

        resfinder_results = amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')

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
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.935, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '1544/1544', msg='Wrong lengths')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S_gyrA_beta-lactam', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'pointfinder_16S_gyrA_beta-lactam.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))
        self.assertEqual(len(records), 2, 'Wrong number of hit records')
        expected_records1 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, 'gyrA-A67P.fsa'), 'fasta'))
        self.assertEqual(expected_records1['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")
        expected_records2 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, '16S_rrsD-1T1065.fsa'), 'fasta'))
        self.assertEqual(expected_records2['16S_rrsD'].seq.upper(), records['16S_rrsD'].seq.upper(),
                         "records don't match")

    def testResfinderPointfinderSalmonella_16Src_C1065T_gyrArc_A67_beta_lactam_Success(self):
        pointfinder_database = PointfinderBlastDatabase(self.pointfinder_dir, 'salmonella')
        blast_handler = BlastHandler(self.resfinder_database, 2, pointfinder_database)
        amr_detection = AMRDetection(self.resfinder_database, blast_handler, pointfinder_database,
                                     output_dir=self.outdir.name)

        file = path.join(self.test_data_dir, "16S-rc_gyrA-rc_beta-lactam.fsa")
        files = [file]
        amr_detection.run_amr_detection(files, 99, 99, 90)

        resfinder_results = amr_detection.get_resfinder_results()
        self.assertEqual(len(resfinder_results.index), 1, 'Wrong number of rows in result')

        result = resfinder_results[resfinder_results['Gene'] == 'blaIMP-42']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.73, places=2, msg='Wrong pid')

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
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.935, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '1544/1544', msg='Wrong lengths')

        result = pointfinder_results[pointfinder_results['Gene'] == 'gyrA (A67P)']
        self.assertEqual(len(result.index), 1, 'Wrong number of results detected')
        self.assertEqual(result.index[0], '16S-rc_gyrA-rc_beta-lactam', msg='Wrong file')
        self.assertEqual(result['Type'].iloc[0], 'codon', msg='Wrong type')
        self.assertEqual(result['Position'].iloc[0], 67, msg='Wrong codon position')
        self.assertEqual(result['Mutation'].iloc[0], 'GCC -> CCC (A -> P)', msg='Wrong mutation')
        self.assertAlmostEqual(result['%Identity'].iloc[0], 99.962, places=3, msg='Wrong pid')
        self.assertAlmostEqual(result['%Overlap'].iloc[0], 100.00, places=2, msg='Wrong overlap')
        self.assertEqual(result['HSP Length/Total Length'].iloc[0], '2637/2637', msg='Wrong lengths')

        hit_file = path.join(self.outdir.name, 'pointfinder_16S-rc_gyrA-rc_beta-lactam.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))
        self.assertEqual(len(records), 2, 'Wrong number of hit records')
        expected_records1 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, 'gyrA-A67P.fsa'), 'fasta'))
        self.assertEqual(expected_records1['gyrA'].seq.upper(), records['gyrA'].seq.upper(), "records don't match")
        expected_records2 = SeqIO.to_dict(SeqIO.parse(path.join(self.test_data_dir, '16S_rrsD-1T1065.fsa'), 'fasta'))
        self.assertEqual(expected_records2['16S_rrsD'].seq.upper(), records['16S_rrsD'].seq.upper(),
                         "records don't match")

    def testResfinderExcludeNonMatches(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database, False,
                                     output_dir=self.outdir.name)
        file_beta_lactam = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        file_non_match = path.join(self.test_data_dir, "non-match.fsa")
        files = [file_beta_lactam, file_non_match]
        amr_detection.run_amr_detection(files, 99, 90, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file_beta_lactam, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testResfinderIncludeNonMatches(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database, True,
                                     output_dir=self.outdir.name)
        file_beta_lactam = path.join(self.test_data_dir, "beta-lactam-blaIMP-42-mut-2.fsa")
        file_non_match = path.join(self.test_data_dir, "non-match.fsa")
        files = [file_beta_lactam, file_non_match]
        amr_detection.run_amr_detection(files, 99, 90, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 2, 'Wrong number of rows in result')

        result_beta_lactam = summary_results.loc['beta-lactam-blaIMP-42-mut-2']
        self.assertTrue(isinstance(result_beta_lactam, pandas.Series), 'Wrong type of results returned')
        self.assertEqual(result_beta_lactam['Genotype'], 'blaIMP-42', 'Wrong genotype')

        result_sensitive = summary_results.loc['non-match']
        self.assertTrue(isinstance(result_sensitive, pandas.Series), 'Wrong type of results returned')
        self.assertEqual(result_sensitive['Genotype'], 'None', 'Wrong genotype')

        hit_file = path.join(self.outdir.name, 'resfinder_beta-lactam-blaIMP-42-mut-2.fsa')
        records = SeqIO.to_dict(SeqIO.parse(hit_file, 'fasta'))

        self.assertEqual(len(records), 1, 'Wrong number of hit records')

        expected_records = SeqIO.to_dict(SeqIO.parse(file_beta_lactam, 'fasta'))
        self.assertEqual(expected_records['blaIMP-42_1_AB753456'].seq, records['blaIMP-42_1_AB753456'].seq,
                         "records don't match")

    def testNonMatches(self):
        amr_detection = AMRDetection(self.resfinder_database, self.blast_handler, self.pointfinder_database, True,
                                     output_dir=self.outdir.name)
        files = [path.join(self.test_data_dir, "non-match.fsa")]
        amr_detection.run_amr_detection(files, 99, 90, 90)

        summary_results = amr_detection.get_summary_results()
        self.assertEqual(len(summary_results.index), 1, 'Wrong number of rows in result')

        result_sensitive = summary_results.loc['non-match']
        self.assertTrue(isinstance(result_sensitive, pandas.Series), 'Wrong number of results detected')
        self.assertEqual(result_sensitive['Genotype'], 'None', 'Wrong genotype')

        self.assertEqual(len(os.listdir(self.outdir.name)), 0, 'File found where none should exist')


if __name__ == '__main__':
    unittest.main()
