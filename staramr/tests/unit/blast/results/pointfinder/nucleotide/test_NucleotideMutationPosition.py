import unittest

from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


class NucleotideMutationPositionTest(unittest.TestCase):

    def testMutationPositionNucleotideStart(self):
        mutation_position = 0
        # @formatter:off
        amr_gene_string = "ATCGATCGA"
        genome_string   = "TTCGATCGA"
        #@formatter:on
        amr_gene_start = 1
        amr_gene_strand = 'plus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'T', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A1T', 'Incorrect string')

    def testMutationPositionNucleotideMiddle(self):
        mutation_position = 4
        # @formatter:off
        amr_gene_string = "ATCGATCGA"
        genome_string   = "ATCGCTCGA"
        #@formatter:on
        amr_gene_start = 1
        amr_gene_strand = 'plus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 5, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'C', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A5C', 'Incorrect string')

    def testMutationPositionNucleotideEnd(self):
        mutation_position = 8
        # @formatter:off
        amr_gene_string = "ATCGATCGA"
        genome_string   = "ATCGATCGG"
        #@formatter:on
        amr_gene_start = 1
        amr_gene_strand = 'plus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 9, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'G', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A9G', 'Incorrect string')

    def testMutationPositionGapStart(self):
        mutation_position = 0
        # @formatter:off
        amr_gene_string = "ATCG"
        genome_string   = "-TCG"
        #@formatter:on
        amr_gene_start = 1
        amr_gene_strand = 'plus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), '-', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A1-', 'Incorrect string')

    def testMutationPositionGapEnd(self):
        mutation_position = 3
        # @formatter:off
        amr_gene_string = "ATCG"
        genome_string   = "ATC-"
        #@formatter:on
        amr_gene_start = 1
        amr_gene_strand = 'plus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'G', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), '-', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'G4-', 'Incorrect string')

    def testMutationPositionGapReference(self):
        mutation_position = 0
        # @formatter:off
        amr_gene_string = "-TCG"
        genome_string   = "ATCG"
        #@formatter:on
        amr_gene_start = 1
        amr_gene_strand = 'plus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), '-', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'A', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), '-1A', 'Incorrect string')

    def testMutationPositionStartDBNegative(self):
        mutation_position = 8
        # @formatter:off
        amr_gene_string = "TCGATCGAT" # rc("ATCGATCGA")
        genome_string   = "TCGATCGAA"
        #@formatter:on
        amr_gene_start = 9
        amr_gene_strand = 'minus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'T', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A1T', 'Incorrect string')

    def testMutationPositionMiddleDBNegative(self):
        mutation_position = 7
        # @formatter:off
        amr_gene_string = "TCGATCGAT" # rc("ATCGATCGA")
        genome_string   = "TCGATCGGT"
        #@formatter:on
        amr_gene_start = 9
        amr_gene_strand = 'minus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'T', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'C', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'T2C', 'Incorrect string')

    def testMutationPositionEndDBNegative(self):
        mutation_position = 0
        # @formatter:off
        amr_gene_string = "TCGATCGAT" # rc("ATCGATCGA")
        genome_string   = "CCGATCGAT"
        #@formatter:on
        amr_gene_start = 9
        amr_gene_strand = 'minus'

        mutation = NucleotideMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start,
                                              amr_gene_strand)

        self.assertEqual(mutation.get_mutation_position(), 9, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_amr_gene_mutation(), 'A', 'Incorrect database mutation')
        self.assertEqual(mutation.get_genome_mutation(), 'G', 'Incorrect query mutation')
        self.assertEqual(mutation.get_mutation_string_short(), 'A9G', 'Incorrect string')
