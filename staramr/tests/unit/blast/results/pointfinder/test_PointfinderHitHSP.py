import unittest

import pandas as pd

from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.exceptions.InvalidPositionException import InvalidPositionException


class PointfinderHitHSPTest(unittest.TestCase):

    def testBuildPointfinderHitHSPSuccess(self):
        blast_record = pd.Series({'sstart': 1, 'send': 10, 'qstart': 1, 'qend': 10, 'sstrand': 'plus'})

        PointfinderHitHSP(file=None, blast_record=blast_record)

    def testBuildPointfinderHitHSPFailInvalidSubjectCoords(self):
        blast_record = pd.Series({'sstart': 10, 'send': 1, 'qstart': 1, 'qend': 10, 'sstrand': 'plus'})

        self.assertRaises(InvalidPositionException, PointfinderHitHSP, None, blast_record)

    def testBuildPointfinderHitHSPInvalidQueryCoords(self):
        blast_record = pd.Series({'sstart': 1, 'send': 10, 'qstart': 10, 'qend': 1, 'sstrand': 'plus'})

        self.assertRaises(InvalidPositionException, PointfinderHitHSP, None, blast_record)

    def testGetMutations_MultipleIndels(self):
        # This tests correctly handling the following sequence of mutations:
        # insertion, insertion, deletion, deletion

        blast_record = pd.Series({
            "qseqid":"pmrA",
            "sseqid":"pmrA",
            "pident":98.222,
            "length":675,
            "qstart":1,
            "qend":669,
            "sstart":1,
            "send":669,
            "slen":669,
            "qlen":669,
            "sstrand":"plus",
            "sseq":"ATGAAAATTCTGATTGTTGAAGACGATCCCACGCTGTTATTGCAGGGACTGATTCTGGCGGCGCAAACCGAAGGCTACGCGTGCGATGGCGTGACAACCGCGCGGATGGCGGAACAAAGCCTTGAGGCAGGTCATTACAGCCTGGTGGTACTGGATTTAGGGTTACCCGACGAAGATGGACTGCATTTTCTCGCCCGTATCCGGCAGAAAAAATATACCCTGCCGGTACTGATCCTCACCAAAGCTCGCGATACGCTGACCGACAAAATCGCCGGGCTGGATGTCGGTGCCGACGACTATCTGGTGAAGCCTTTTGCGCTGGAAGAGTTACATGCCCGTATCCGCGCCCTGCTACGACGCCATAATAATCAGGGCGAAAGTGAGCTGATTGTTGGCAATCTGACGCTGAACATGGGTCGCCGTCAGGTATGGATGGGCGGTGAAGAGTTGATT---ACGCCCAAAGAATATGCTCTGCTGTCACGGTTAATGCTCAAAGCAGGCAGTCCGGTGCATCGGGAAATTCTCTACAACGACATCTATAACTGGGACAATGAACCCTCGACCAACACCCTGGAAGTGCATATCCACAAT---CGCGACAAAGTGGGCAAAGCCCGTATCCGCACCGTGCGCGGCTTTGGCTATATGCTGGTCGCGAATGAGGAAAACTAA",
            "qseq":"ATGAAAATTCTGATTGTTGAAGACGAT---ACGCTGTTATTGCAGGGACTGATTCTGGCGGCGCAAACCGAAGGCTACGCGTGCGATGGCGTGACAACCGCGCGGATGGCGGAACAAAGCCTTGAGGCAGGTCATTACAGCCTGGTGGTACTGGATTTAGGGTTACCCGACGAAGATGGACTGCATTTTCTCGCCCGTATCCGGCAGAAAAAATATACCCTGCCGGTACTGATCCTCACC---GCTCGCGATACGCTGACCGACAAAATCGCCGGGCTGGATGTCGGTGCCGACGACTATCTGGTGAAGCCTTTTGCGCTGGAAGAGTTACATGCCCGTATCCGCGCCCTGCTACGACGCCATAATAATCAGGGCGAAAGTGAGCTGATTGTTGGCAATCTGACGCTGAACATGGGTCGCCGTCAGGTATGGATGGGCGGTGAAGAGTTGATTCTGACGCCCAAAGAATATGCTCTGCTGTCACGGTTAATGCTCAAAGCAGGCAGTCCGGTGCATCGGGAAATTCTCTACAACGACATCTATAACTGGGACAATGAACCCTCGACCAACACCCTGGAAGTGCATATCCACAATCTGCGCGACAAAGTGGGCAAAGCCCGTATCCGCACCGTGCGCGGCTTTGGCTATATGCTGGTCGCGAATGAGGAAAACTAA",
            "plength":100.89686098654708})
        
        hit = PointfinderHitHSP("staramr/tests/unit/data/pmrA-multi-indel.fsa", blast_record)
        mutations = hit.get_mutations()

        self.assertEqual(mutations[0]._codon_start, 9, msg='Wrong codon position.')
        self.assertEqual(mutations[0]._database_amr_gene_codon, "---", msg='Wrong AMR codon.')
        self.assertEqual(mutations[0]._input_genome_codon, "CCC", msg='Wrong input codon.')

        self.assertEqual(mutations[1]._codon_start, 79, msg='Wrong codon position.')
        self.assertEqual(mutations[1]._database_amr_gene_codon, "---", msg='Wrong AMR codon.')
        self.assertEqual(mutations[1]._input_genome_codon, "AAA", msg='Wrong input codon.')

        self.assertEqual(mutations[2]._codon_start, 150, msg='Wrong codon position.')
        self.assertEqual(mutations[2]._database_amr_gene_codon, "CTG", msg='Wrong AMR codon.')
        self.assertEqual(mutations[2]._input_genome_codon, "---", msg='Wrong input codon.')

        self.assertEqual(mutations[3]._codon_start, 197, msg='Wrong codon position.')
        self.assertEqual(mutations[3]._database_amr_gene_codon, "CTG", msg='Wrong AMR codon.')
        self.assertEqual(mutations[3]._input_genome_codon, "---", msg='Wrong input codon.')
