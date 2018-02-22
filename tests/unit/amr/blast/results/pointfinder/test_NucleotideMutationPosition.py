import unittest

from amr.blast.results.pointfinder.NucleotideMutationPosition import NucleotideMutationPosition


class AMRDetectionIT(unittest.TestCase):

    def testMutationPositionStartCodon1(self):
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

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TTC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'F', 'Incorrect query amino acid')

    def testMutationPositionMiddleCodon1(self):
        mutation_position = 1
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "AGCGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AGC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'S', 'Incorrect query amino acid')

    def testMutationPositionEndCodon1(self):
        mutation_position = 2
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATGGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATG', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'M', 'Incorrect query amino acid')

    def testMutationPositionStartCodon2(self):
        mutation_position = 3
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCAATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AAT', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'N', 'Incorrect query amino acid')

    def testMutationPositionEndCodon2(self):
        mutation_position = 5
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCGACCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 6, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'GAC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'D', 'Incorrect query amino acid')

    def testMutationPositionStartCodon3(self):
        mutation_position = 6
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATCGATGGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 7, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 3, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'CGA', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'GGA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'R', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'G', 'Incorrect query amino acid')

    def testMutationPositionStartCodon1StartMethionine(self):
        mutation_position = 0
        # @formatter:off
        database_string = "ATCGATCGA"
        query_string    = "ATGGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATG', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'M', 'Incorrect query amino acid')

    def testMutationPositionStartCodon1Stop(self):
        mutation_position = 2
        # @formatter:off
        database_string = "TACGATCGA"
        query_string    = "TAAGATCGA"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'TAC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'TAA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'Y', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), '*', 'Incorrect query amino acid')

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

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '-TC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

    def testMutationPositionGapMiddle(self):
        mutation_position = 1
        # @formatter:off
        database_string = "ATCG"
        query_string    = "A-CG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'A-C', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

    def testMutationPositionGapEnd(self):
        mutation_position = 2
        # @formatter:off
        database_string = "ATCG"
        query_string    = "AT-G"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AT-', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

    def testMutationPositionGapMiddleEnd(self):
        mutation_position = 2
        # @formatter:off
        database_string = "ATCGG"
        query_string    = "AT--G"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'AT-', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

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
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

    def testMutationPositionGapPreviousCodon(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGACT"
        query_string    = "CC----GACT"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

    def testMutationPositionGapLargerPreviousCodon(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGACTT"
        query_string    = "C-----GACTT"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')


    def testMutationPositionGapBefore(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGAC"
        query_string    = "-CCA--GAC"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'A--', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')

    def testMutationPositionGapBeforeAfter(self):
        mutation_position = 3
        # @formatter:off
        database_string = "CCCATCGACT"
        query_string    = "-CCA--GA-T"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'A--', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'X', 'Incorrect query amino acid')


    def testMutationPositionGapReferenceStart(self):
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

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), '-TC', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'I', 'Incorrect query amino acid')


    def testMutationPositionGapReferenceMiddle(self):
        mutation_position = 1
        # @formatter:off
        database_string = "A-CG"
        query_string    = "ATCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'A-C', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'I', 'Incorrect query amino acid')


    def testMutationPositionGapReferenceEnd(self):
        mutation_position = 2
        # @formatter:off
        database_string = "AT-G"
        query_string    = "ATCG"
        #@formatter:on
        database_start = 1
        database_frame = 1
        query_frame = 1
        mutation = NucleotideMutationPosition(mutation_position, database_string, query_string, database_start,
                                              database_frame, query_frame)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_database_codon(), 'AT-', 'Incorrect database codon')
        self.assertEqual(mutation.get_query_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amino_acid(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_query_amino_acid(), 'I', 'Incorrect query amino acid')