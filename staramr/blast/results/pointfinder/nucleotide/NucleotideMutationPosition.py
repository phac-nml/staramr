import Bio.Seq

from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class NucleotideMutationPosition(MutationPosition):

    def __init__(self, match_position, amr_gene_string, genome_string, amr_gene_start, amr_gene_strand):
        """
        Creates a new NucleotideMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param amr_gene_string: The amr gene string.
        :param genome_string: The genome string.
        :param amr_gene_start: The start coordinates of the BLAST database hit.
        :param amr_gene_strand: The strand of the BLAST database.
        """
        super().__init__(match_position, amr_gene_start, amr_gene_strand)

        self._amr_gene_nucleotide = amr_gene_string[match_position].upper()
        self._genome_nucleotide = genome_string[match_position].upper()

        if self._amr_gene_strand == 'minus':
            self._amr_gene_nucleotide = Bio.Seq.reverse_complement(self._amr_gene_nucleotide)
            self._genome_nucleotide = Bio.Seq.reverse_complement(self._genome_nucleotide)

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
        return "[amr_gene_start=" + str(self._amr_gene_start) + ", amr_gene_strand=" + str(
            self._amr_gene_strand) + ", nucleotide_position=" \
               + str(self._nucleotide_position_amr_gene) + ", mutation_start=" + str(self.get_mutation_position()) \
               + ", mutation=" + self.get_mutation_string() + "]"
