import unittest
import pandas

from amr.blast.pointfinder.PointfinderDatabaseInfo import PointfinderDatabaseInfo
from amr.blast.results.pointfinder.NucleotideMutationPosition import NucleotideMutationPosition

class PointfinderDatabaseInfoTest(unittest.TestCase):

    def setUp(self):
        pandas_pointfinder_table = pandas.DataFrame([
                ['gyrA', 'gyrA', 1, 1, 'ATC', 'I', 'F', 'Quinolones', 15848289],
                ['gyrA', 'gyrA', 1, 2, 'GAT', 'D', 'N,H', 'Quinolones', 15848289],
            ],
            columns=('#Gene_ID', 'Gene_name', 'No of mutations needed', 'Codon_pos', 'Ref_nuc', 'Ref_codon', 'Res_codon', 'Resistance', 'PMID'))

        self.database = PointfinderDatabaseInfo.from_pandas_table(pandas_pointfinder_table)

        mutation_position  = 0
        database_string = "ATCGATCGA"
        query_string    = "TTCGATCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        self.mutation1 = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        mutation_position  = 3
        database_string = "ATCGATCGA"
        query_string    = "ATCAATCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        self.mutation2 = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        mutation_position  = 8
        database_string = "ATCGATCGA"
        query_string    = "ATCGATCGT"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        self.mutation_missing = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        mutation_position  = 3
        database_string = "ATCGATCGA"
        query_string    = "ATCGAACGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        self.mutation_aa_not_match = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)


    def testGetResistanceCodons1Mutation1Codon(self):
        resistance_mutations = self.database.get_resistance_codons('gyrA', [self.mutation1])

        self.assertEqual(resistance_mutations, [self.mutation1], "Did not pick up correct mutations")


    def testGetResistanceCodons1Mutation2Codons(self):
        resistance_mutations = self.database.get_resistance_codons('gyrA', [self.mutation2])

        self.assertEqual(resistance_mutations, [self.mutation2], "Did not pick up correct mutations")


    def testGetResistanceCodons2Mutation2Codons(self):
        resistance_mutations = self.database.get_resistance_codons('gyrA', [self.mutation1, self.mutation2])

        self.assertEqual(resistance_mutations, [self.mutation1, self.mutation2], "Did not pick up correct mutations")


    def testGetResistanceCodons2Mutation1Missing(self):
        resistance_mutations = self.database.get_resistance_codons('gyrA', [self.mutation1, self.mutation_missing])

        self.assertEqual(resistance_mutations, [self.mutation1], "Did not pick up correct mutations")


    def testGetResistanceCodons1Mutation1Missing(self):
        resistance_mutations = self.database.get_resistance_codons('gyrA', [self.mutation_missing])

        self.assertEqual(resistance_mutations, [], "Did not pick up correct mutations")


    def testGetResistanceCodons1Mutation1MissingAANotMatch(self):
        resistance_mutations = self.database.get_resistance_codons('gyrA', [self.mutation_aa_not_match])

        self.assertEqual(resistance_mutations, [], "Did not pick up correct mutations")


    def testGetResfinderPhenotype(self):
        phenotype = self.database.get_phenotype('gyrA', self.mutation1)

        self.assertEqual(phenotype, 'Quinolones')


    def testGetResfinderPhenotypeMissingFail(self):
        self.assertRaises(Exception, self.database.get_phenotype, 'gyrA', self.mutation_missing)