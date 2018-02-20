import unittest

from amr.blast.results.pointfinder.NucleotideMutationPosition import NucleotideMutationPosition

class AMRDetectionIT(unittest.TestCase):


    def testMutationPositionStartCodon1(self):
        mutation_position  = 0
        database_string = "ATCGATCGA"
        query_string    = "TTCGATCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TTC', 'Incorrect query codon')


    def testMutationPositionMiddleCodon1(self):
        mutation_position  = 1
        database_string = "ATCGATCGA"
        query_string    = "AGCGATCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AGC', 'Incorrect query codon')


    def testMutationPositionEndCodon1(self):
        mutation_position  = 2
        database_string = "ATCGATCGA"
        query_string    = "ATGGATCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATG', 'Incorrect query codon')


    def testMutationPositionStartCodon2(self):
        mutation_position  = 3
        database_string = "ATCGATCGA"
        query_string    = "ATCAATCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AAT', 'Incorrect query codon')


    def testMutationPositionEndCodon2(self):
        mutation_position  = 5
        database_string = "ATCGATCGA"
        query_string    = "ATCGACCGA"
        database_start  = 1
        database_frame  = 1
        query_frame     = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start, database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 6, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'GAC', 'Incorrect query codon')