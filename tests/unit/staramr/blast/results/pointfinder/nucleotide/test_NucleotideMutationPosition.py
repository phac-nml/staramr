import unittest

from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


class NucleotideMutationPositionTest(unittest.TestCase):

    def testMutationPositionNucleotideStart(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "TTCGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'T', 'Incorrect query mutation')

    def testMutationPositionNucleotideMiddle(self):
        mutation_position = 4
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCGCTCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 5, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'C', 'Incorrect query mutation')

    def testMutationPositionNucleotideEnd(self):
        mutation_position = 8
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCGATCGG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 9, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'G', 'Incorrect query mutation')

    def testMutationPositionGapStart(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCG"
        query_string    = "-TCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), '-', 'Incorrect query mutation')

    def testMutationPositionGapEnd(self):
        mutation_position = 3
        # @formatter:off
        database_string = "ATCG"
        query_string    = "ATC-"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'G', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), '-', 'Incorrect query mutation')

    def testMutationPositionGapReference(self):
        mutation_position = 0
        # @formatter:off
        database_string = "-TCG"
        query_string    = "ATCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), '-', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'A', 'Incorrect query mutation')

    def testMutationPositionStartDBNegative(self):
        mutation_position = 8
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "TCGATCGAA"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'T', 'Incorrect query mutation')

    def testMutationPositionMiddleDBNegative(self):
        mutation_position = 7
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "TCGATCGGT"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'T', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'C', 'Incorrect query mutation')

    def testMutationPositionEndDBNegative(self):
        mutation_position = 0
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "CCGATCGAT"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 9, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_query_mutation(), 'G', 'Incorrect query mutation')

