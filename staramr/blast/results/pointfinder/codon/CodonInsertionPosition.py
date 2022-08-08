import math

from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition

"""
A class defining a codon-based insertion mutation for PointFinder.

This class is needed to handle the special case of codon insertions needing to report
the codon to the LEFT of the insertion, not the sequence position of the insertion itself.
"""


class CodonInsertionPosition(CodonMutationPosition):

    def __init__(self, match_position, database_amr_gene_string, input_genome_blast_string, database_amr_gene_start, offset=0):
        """
        Creates a new CodonMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_amr_gene_string: The database amr gene string from BLAST.
        :param input_genome_blast_string: The genome BLAST string from the input genome.
        :param database_amr_gene_start: The start coordinates of the BLAST amr gene hit.
        :param offset: The amount to offset the mutation by (important for promoter mutations).
        """
        super().__init__(match_position, database_amr_gene_string, input_genome_blast_string, database_amr_gene_start, offset)

        # We need to correctly set the codon position to be 1 to the left (-1) of the observed insertion codon.
        self._codon_start = math.ceil(self._nucleotide_position_amr_gene / 3) - 1

    def __repr__(self):
        return (
            'CodonInsertionPosition(_database_amr_gene_start={_database_amr_gene_start}, _nucleotide_position_amr_gene={_nucleotide_position_amr_gene}, '
            '_codon_start={_codon_start}, _database_amr_gene_codon={_database_amr_gene_codon}, _input_genome_codon={_input_genome_codon})').format(
            **self.__dict__)

    def get_pointfinder_mutation_string(self):
        # Since the position codon insertions in the Pointfinder database is off by one, we need
        # to adjust the position in the report to show what Pointfinder shows.
        return self.get_database_amr_gene_mutation() + str(
            self.get_mutation_position() + 1) + self.get_input_genome_mutation()

