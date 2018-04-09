import unittest

from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition


class CodonMutationPositionTest(unittest.TestCase):

    def testMutationPositionStartCodon1(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "TTCGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TTC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'F', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1F', 'Incorrect string')

    def testMutationPositionMiddleCodon1(self):
        mutation_position = 1
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "AGCGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AGC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'S', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1S', 'Incorrect string')

    def testMutationPositionEndCodon1(self):
        mutation_position = 2
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATGGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATG', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'M', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1M', 'Incorrect string')

    def testMutationPositionStartCodon2(self):
        mutation_position = 3
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCAATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AAT', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'N', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'D2N', 'Incorrect string')

    def testMutationPositionEndCodon2(self):
        mutation_position = 5
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCGACCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 6, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'GAC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'D', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'D2D', 'Incorrect string')

    def testMutationPositionStartCodon3(self):
        mutation_position = 6
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCGATGGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 7, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 3, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 3, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'CGA', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'GGA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'R', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'G', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'R3G', 'Incorrect string')

    def testMutationPositionStartCodon1StartMethionine(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATGGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATG', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'M', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1M', 'Incorrect string')

    def testMutationPositionStartCodon1Stop(self):
        mutation_position = 2
        # @formatter:off
        database_string = "TACGATCGA"
        query_string    = "TAAGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'TAC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TAA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'Y', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), '*', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'Y1*', 'Incorrect string')

    def testMutationPositionGapStart(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCG"
        query_string    = "-TCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '-TC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapMiddle(self):
        mutation_position = 1
        # @formatter:off
        database_string = "ATCG"
        query_string    = "A-CG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'A-C', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapEnd(self):
        mutation_position = 2
        # @formatter:off
        database_string = "ATCG"
        query_string    = "AT-G"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AT-', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapMiddleEnd(self):
        mutation_position = 2
        # @formatter:off
        database_string = "ATCGG"
        query_string    = "AT--G"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AT-', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapStartMiddleEnd(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGAC"
        query_string    = "CCC---GAC"
        #@formatter:on
        # @formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapPreviousCodon(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGACT"
        query_string    = "CC----GACT"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapLargerPreviousCodon(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGACTT"
        query_string    = "C-----GACTT"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapBefore(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGAC"
        query_string    = "-CCA--GAC"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'A--', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapBeforeAfter(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGACT"
        query_string    = "-CCA--GA-T"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'A--', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapReferenceStart(self):
        mutation_position = 0
        # @formatter:off
        database_string = "-TCG"
        query_string    = "ATCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), '-TC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'X1I', 'Incorrect string')

    def testMutationPositionGapReferenceMiddle(self):
        mutation_position = 1
        # @formatter:off
        database_string = "A-CG"
        query_string    = "ATCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'A-C', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'X1I', 'Incorrect string')

    def testMutationPositionGapReferenceEnd(self):
        mutation_position = 2
        # @formatter:off
        database_string = "AT-G"
        query_string    = "ATCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'AT-', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'X1I', 'Incorrect string')

    def testMutationPositionStartCodon1DBNegative(self):
        mutation_position = 8
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "TCGATCGAA"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TTC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'F', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1F', 'Incorrect string')

    def testMutationPositionMiddleCodon1DBNegative(self):
        mutation_position = 7
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "TCGATCGGT"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ACC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'T', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1T', 'Incorrect string')

    def testMutationPositionEndCodon1DBNegative(self):
        mutation_position = 6
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "TCGATCTAT"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1I', 'Incorrect string')

    def testMutationPositionStartCodon2DBNegative(self):
        mutation_position = 5
        # @formatter:off
        database_string = "TCGATCGAT" # rc("ATCGATCGA")
        query_string    = "TCGATTGAT"
        #@formatter:on
        database_start = 9
        database_frame = -1
        query_frame = 1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AAT', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'N', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'D2N', 'Incorrect string')

    def testMutationPositionStartCodon1QNegative(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "TTCGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = -1
        mutation = CodonMutationPosition(mutation_position, database_string, query_string, database_start,
                                         database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TTC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_mutation(), 'F', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1F', 'Incorrect string')
