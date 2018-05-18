import math

import Bio.Seq

from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a codon-based mutation for PointFinder.
"""


class CodonMutationPosition(MutationPosition):

    def __init__(self, match_position, amr_gene_string, genome_string, amr_gene_start, amr_gene_strand):
        """
        Creates a new CodonMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param amr_gene_string: The amr gene string from BLAST.
        :param genome_string: The genome BLAST string.
        :param amr_gene_start: The start coordinates of the BLAST amr gene hit.
        :param amr_gene_strand: The strand of the amr gene.
        """
        super().__init__(match_position, amr_gene_start, amr_gene_strand)

        self._codon_start = math.ceil(self._nucleotide_position_amr_gene / 3)
        frame_shift = (self._nucleotide_position_amr_gene - 1) % 3

        self._amr_gene_codon = self._find_codon(amr_gene_string, match_position, amr_gene_strand, frame_shift)
        self._genome_codon = self._find_codon(genome_string, match_position, amr_gene_strand, frame_shift)

    def _find_codon(self, nucleotides, match_position, strand, frame_shift):
        if strand == 'plus':
            codon_start_index = match_position - frame_shift
            return nucleotides[codon_start_index:(codon_start_index + 3)].upper()
        else:
            codon_end_index = match_position + frame_shift
            return Bio.Seq.reverse_complement(nucleotides[(codon_end_index - 3 + 1):(codon_end_index + 1)].upper())

    def get_codon_start(self):
        """
        Gets the codon start for PointFinder (1-based).
        :return: The codon start.
        """
        return self._codon_start

    def get_amr_gene_codon(self):
        """
        Gets the particular codon from the amr gene.
        :return: The codon.
        """
        return self._amr_gene_codon

    def get_amr_gene_amino_acid(self):
        """
        Gets the corresponding amino acid from the amr gene. If there is an indel, returns 'X'.
        :return: The amino acid from the amr gene.
        """
        if '-' in self.get_amr_gene_codon():
            return 'X'
        else:
            return Bio.Seq.translate(self.get_amr_gene_codon(), table='Standard')

    def get_genome_amino_acid(self):
        """
        Gets the corresponding amino acid from the genome.  If there is an indel returns 'X'.
        :return: The amino acid from the genome.
        """
        if '-' in self.get_genome_codon():
            return 'X'
        else:
            return Bio.Seq.translate(self.get_genome_codon(), table='Standard')

    def get_genome_codon(self):
        """
        Gets the codon from the genome.
        :return: The codon from the genome.
        """
        return self._genome_codon

    def get_mutation_position(self):
        return self.get_codon_start()

    def get_mutation_string(self):
        return self.get_amr_gene_codon() + ' -> ' + self.get_genome_codon() + ' (' + self.get_amr_gene_amino_acid() \
               + ' -> ' + self.get_genome_amino_acid() + ')'

    def get_amr_gene_mutation(self):
        return self.get_amr_gene_amino_acid().upper()

    def get_genome_mutation(self):
        return self.get_genome_amino_acid().upper()

    def get_type(self):
        return 'codon'

    def __repr__(self):
        return "[amr_gene_start=" + str(self._amr_gene_start) + ", amr_gene_strand=" + str(
            self._amr_gene_strand) + ", nucleotide_position=" \
               + str(self._nucleotide_position_amr_gene) + ", codon_start=" + str(self._codon_start) \
               + ", codon=" + self._amr_gene_codon + "]"
