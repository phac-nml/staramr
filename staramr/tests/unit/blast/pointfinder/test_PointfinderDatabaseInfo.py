import unittest

import pandas as pd

from staramr.blast.pointfinder.PointfinderDatabaseInfo import PointfinderDatabaseInfo
from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition


class PointfinderDatabaseInfoTest(unittest.TestCase):

    def setUp(self):
        pandas_pointfinder_table = pd.DataFrame([
            ['gyrA', 'gyrA', 1, 1, 'ATC', 'I', 'F', 'Quinolones', 15848289],
            ['gyrA', 'gyrA', 1, 2, 'GAT', 'D', 'N,H', 'Quinolones', 15848289],
        ],
            columns=(
                '#Gene_ID', 'Gene_name', 'No of mutations needed', 'Codon_pos', 'Ref_nuc', 'Ref_codon', 'Res_codon',
                'Resistance', 'PMID'))

        self.database = PointfinderDatabaseInfo.from_pandas_table(pandas_pointfinder_table)

        mutation_position = 0
        amr_gene_string = "ATCGATCGA"
        genome_string = "TTCGATCGA"
        amr_gene_start = 1
        self.mutation1 = CodonMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start)

        mutation_position = 3
        amr_gene_string = "ATCGATCGA"
        genome_string = "ATCAATCGA"
        amr_gene_start = 1
        self.mutation2 = CodonMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start)

        mutation_position = 8
        amr_gene_string = "ATCGATCGA"
        genome_string = "ATCGATCGT"
        amr_gene_start = 1
        self.mutation_missing = CodonMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start)

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

    def testGetResistanceCodons1MutationAANotMatch(self):
        mutation_position = 3
        amr_gene_string = "ATCGATCGA"
        genome_string = "ATCGAACGA"
        amr_gene_start = 1
        mutation_aa_not_match = CodonMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start)
        resistance_mutations = self.database.get_resistance_codons('gyrA', [mutation_aa_not_match])

        self.assertEqual(resistance_mutations, [], "Did not pick up correct mutations")

    def testGetResistanceCodons1MutationStartCodon(self):
        mutation_position = 0
        amr_gene_string = "ATCGATCGA"
        genome_string = "ATGGATCGA"
        amr_gene_start = 1
        mutation_start_methionine = CodonMutationPosition(mutation_position, amr_gene_string, genome_string,
                                                          amr_gene_start)
        resistance_mutations = self.database.get_resistance_codons('gyrA', [mutation_start_methionine])

        self.assertEqual(resistance_mutations, [], "Did not pick up correct mutations")

    def testGetResistanceCodons1MutationStopCodon(self):
        mutation_position = 2
        amr_gene_string = "TACGATCGA"
        genome_string = "TAAGATCGA"
        amr_gene_start = 1
        mutation_stop = CodonMutationPosition(mutation_position, amr_gene_string, genome_string, amr_gene_start)
        resistance_mutations = self.database.get_resistance_codons('gyrA', [mutation_stop])

        self.assertEqual(resistance_mutations, [], "Did not pick up correct mutations")

    def testGetResfinderPhenotype(self):
        phenotype = self.database.get_phenotype('gyrA', self.mutation1)

        self.assertEqual(phenotype, 'Quinolones')

    def testGetResfinderPhenotypeMissingFail(self):
        self.assertRaises(Exception, self.database.get_phenotype, 'gyrA', self.mutation_missing)
