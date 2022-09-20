from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class NucleotideMutationPosition(MutationPosition):

    def __init__(self, match_position, database_amr_gene_string, input_genome_string, database_amr_gene_start, offset=0):
        """
        Creates a new NucleotideMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_amr_gene_string: The amr gene string.
        :param input_genome_string: The input genome string.
        :param database_amr_gene_start: The start coordinates of the BLAST database hit.
        :param offset: The amount to offset the mutation by (particularly for negative-coordinate promoter mutations).
        """
        super().__init__(match_position - offset, database_amr_gene_start)

        self._database_amr_gene_mutation = database_amr_gene_string[match_position].upper()
        self._input_genome_mutation = input_genome_string[match_position].upper()

    def get_type(self):
        return 'nucleotide'

    def get_mutation_position(self):
        return self.get_nucleotide_position()

    def get_database_amr_gene_mutation(self):
        # If there's an insertion ('-'), then return 'ins' instead of '-':
        if '-' in self._database_amr_gene_mutation:
            return 'ins'
        # If there's a deletion in the INPUT sequence, then we need to return 'del'
        # instead of the actual AMR reference sequence
        # (which would normally be some sequence)
        if '-' in self._input_genome_mutation:
            return 'del'
        # If neither an insertion or deletion, return the actual AMR gene sequence:
        else:
            return self._database_amr_gene_mutation

    def get_input_genome_mutation(self):
        # If there's a deletion in the input genome, the database format actually
        # requires us to return the reference sequence that was deleted,
        # not the '-' in the sequence that you would expect:
        if '-' in self._input_genome_mutation:
            return self._database_amr_gene_mutation
        # If not a deletion, return the actual input genome sequence:
        else:
            return self._input_genome_mutation

    def get_mutation_string(self):
        return self.get_database_amr_gene_mutation() + ' -> ' + self.get_input_genome_mutation()

    def __repr__(self):
        return (
            'NucleotideMutationPosition(_database_amr_gene_start={_database_amr_gene_start}, _nucleotide_position_amr_gene={_nucleotide_position_amr_gene}, '
            '_database_amr_gene_mutation={_database_amr_gene_mutation}, _input_genome_mutation={_input_genome_mutation})').format(
            **self.__dict__)
