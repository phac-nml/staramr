from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class NucleotideMutationPosition(MutationPosition):

    def __init__(self, match_position, amr_gene_string, genome_string, amr_gene_start):
        """
        Creates a new NucleotideMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param amr_gene_string: The amr gene string.
        :param genome_string: The genome string.
        :param amr_gene_start: The start coordinates of the BLAST database hit.
        """
        super().__init__(match_position, amr_gene_start)

        self._amr_gene_nucleotide = amr_gene_string[match_position].upper()
        self._genome_nucleotide = genome_string[match_position].upper()

    def get_type(self):
        return 'nucleotide'

    def get_mutation_position(self):
        return self.get_nucleotide_position()

    def get_amr_gene_mutation(self):
        return self._amr_gene_nucleotide

    def get_genome_mutation(self):
        return self._genome_nucleotide

    def get_mutation_string(self):
        return self.get_amr_gene_mutation() + ' -> ' + self.get_genome_mutation()

    def __repr__(self):
        return ('NucleotideMutationPosition(_amr_gene_start={_amr_gene_start}, _nucleotide_position_amr_gene={_nucleotide_position_amr_gene}, '
            '_amr_gene_nucleotide={_amr_gene_nucleotide}, _genome_nucleotide={_genome_nucleotide})').format(**self.__dict__)
