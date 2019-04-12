import unittest

from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition


class CodonMutationPositionTest(unittest.TestCase):

    def testMutationPositionStartCodon1(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "TTCGATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'TTC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'F', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1F', 'Incorrect string')

    def testMutationPositionMiddleCodon1(self):
        mutation_position = 1
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "AGCGATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'AGC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'S', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1S', 'Incorrect string')

    def testMutationPositionEndCodon1(self):
        mutation_position = 2
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATGGATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'ATG', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'M', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1M', 'Incorrect string')

    def testMutationPositionStartCodon2(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATCAATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'AAT', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'N', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'D2N', 'Incorrect string')

    def testMutationPositionEndCodon2(self):
        mutation_position = 5
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATCGACCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 6, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'GAT', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'GAC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'D', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'D', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'D2D', 'Incorrect string')

    def testMutationPositionStartCodon3(self):
        mutation_position = 6
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATCGATGGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 7, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 3, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 3, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'CGA', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'GGA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'R', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'G', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'R3G', 'Incorrect string')

    def testMutationPositionStartCodon1StartMethionine(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "ATCGATCGA"
        input_genome_string = "ATGGATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'ATG', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'M', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1M', 'Incorrect string')

    def testMutationPositionStartCodon1Stop(self):
        mutation_position = 2
        # @formatter:off
        database_amr_gene_string = "TACGATCGA"
        input_genome_string = "TAAGATCGA"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'TAC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'TAA', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'Y', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), '*', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'Y1*', 'Incorrect string')

    def testMutationPositionGapStart(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "ATCG"
        input_genome_string = "-TCG"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), '-TC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapMiddle(self):
        mutation_position = 1
        # @formatter:off
        database_amr_gene_string = "ATCG"
        input_genome_string = "A-CG"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'A-C', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapEnd(self):
        mutation_position = 2
        # @formatter:off
        database_amr_gene_string = "ATCG"
        input_genome_string = "AT-G"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'AT-', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapMiddleEnd(self):
        mutation_position = 2
        # @formatter:off
        database_amr_gene_string = "ATCGG"
        input_genome_string = "AT--G"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'AT-', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I1X', 'Incorrect string')

    def testMutationPositionGapStartMiddleEnd(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "CCCATCGAC"
        input_genome_string = "CCC---GAC"
        # @formatter:on
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapPreviousCodon(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "CCCATCGACT"
        input_genome_string = "CC----GACT"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapLargerPreviousCodon(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "CCCATCGACTT"
        input_genome_string = "C-----GACTT"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), '---', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapBefore(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "CCCATCGAC"
        input_genome_string = "-CCA--GAC"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'A--', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapBeforeAfter(self):
        mutation_position = 3
        # @formatter:off
        database_amr_gene_string = "CCCATCGACT"
        input_genome_string = "-CCA--GA-T"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 4, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 2, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 2, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'ATC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'A--', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'I', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'X', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'I2X', 'Incorrect string')

    def testMutationPositionGapReferenceStart(self):
        mutation_position = 0
        # @formatter:off
        database_amr_gene_string = "-TCG"
        input_genome_string = "ATCG"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 1, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), '-TC', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'X1I', 'Incorrect string')

    def testMutationPositionGapReferenceMiddle(self):
        mutation_position = 1
        # @formatter:off
        database_amr_gene_string = "A-CG"
        input_genome_string = "ATCG"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 2, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'A-C', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'X1I', 'Incorrect string')

    def testMutationPositionGapReferenceEnd(self):
        mutation_position = 2
        # @formatter:off
        database_amr_gene_string = "AT-G"
        input_genome_string = "ATCG"
        # @formatter:on
        amr_gene_start = 1

        mutation = CodonMutationPosition(mutation_position, database_amr_gene_string, input_genome_string,
                                         amr_gene_start)

        self.assertEqual(mutation.get_nucleotide_position(), 3, 'Incorrect nucleotide position')
        self.assertEqual(mutation.get_codon_start(), 1, 'Incorrect codon start')
        self.assertEqual(mutation.get_mutation_position(), 1, 'Incorrect mutation start')
        self.assertEqual(mutation.get_database_amr_gene_codon(), 'AT-', 'Incorrect database codon')
        self.assertEqual(mutation.get_input_genome_codon(), 'ATC', 'Incorrect query codon')
        self.assertEqual(mutation.get_database_amr_gene_mutation(), 'X', 'Incorrect database amino acid')
        self.assertEqual(mutation.get_input_genome_mutation(), 'I', 'Incorrect query amino acid')
        self.assertEqual(mutation.get_mutation_string_short(), 'X1I', 'Incorrect string')
