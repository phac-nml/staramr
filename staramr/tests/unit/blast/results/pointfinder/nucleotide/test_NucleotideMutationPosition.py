import unittest

from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


class NucleotideMutationPositionTest(unittest.TestCase):

    def testMutationPositionNucleotideStart(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "TTCGATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = NucleotideMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                              amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_input_genome_mutation(), 'T', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A1T', 'Incorrect string')

    def testMutationPositionNucleotideMiddle(self):
        mutation_position = 4
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATCGCTCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = NucleotideMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                              amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 5, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_input_genome_mutation(), 'C', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A5C', 'Incorrect string')

    def testMutationPositionNucleotideEnd(self):
        mutation_position = 8
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATCGATCGG"
        # @formatter:on
        amr_gene_start = 1

        mutation = NucleotideMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                              amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 9, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_input_genome_mutation(), 'G', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A9G', 'Incorrect string')

    def testMutationPositionGapStart(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "ATCG"
        input_genome_string = "-TCG"
        # @formatter:on
        amr_gene_start = 1

        mutation = NucleotideMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                              amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_input_genome_mutation(), '-', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A1-', 'Incorrect string')

    def testMutationPositionGapEnd(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "ATCG"
        input_genome_string = "ATC-"
        # @formatter:on
        amr_gene_start = 1

        mutation = NucleotideMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                              amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'G', 'Incorrect database mutation')
        self.assertEqual(mutation.get_input_genome_mutation(), '-', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'G4-', 'Incorrect string')

    def testMutationPositionGapReference(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "-TCG"
        input_genome_string = "ATCG"
        # @formatter:on
        amr_gene_start = 1

        mutation = NucleotideMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                              amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), '-', 'Incorrect database mutation')
        self.assertEqual(mutation.get_input_genome_mutation(), 'A', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), '-1A', 'Incorrect string')
