import math

import Bio.Seq

from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a codon-based mutation for PointFinder.
"""


class CodonMutationPosition(MutationPosition):

    def __init__(self, match_position, database_amr_gene_string, input_genome_blast_string, database_amr_gene_start):
        """
        Creates a new CodonMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_amr_gene_string: The database amr gene string from BLAST.
        :param input_genome_blast_string: The genome BLAST string from the input genome.
        :param database_amr_gene_start: The start coordinates of the BLAST amr gene hit.
        """
        super().__init__(match_position, database_amr_gene_start)

        self._codon_start = math.ceil(self._nucleotide_position_amr_gene / 3)
        frame_shift = (self._nucleotide_position_amr_gene - 1) % 3

        self._database_amr_gene_codon = self._find_codon(database_amr_gene_string, match_position, frame_shift)
        self._input_genome_codon = self._find_codon(input_genome_blast_string, match_position, frame_shift)

    def _find_codon(self, nucleotides, match_position, frame_shift):
        codon_start_index = match_position - frame_shift
        return nucleotides[codon_start_index:(codon_start_index + 3)].upper()

    def get_codon_start(self):
        """
        Gets the codon start for PointFinder (1-based).
        :return: The codon start.
        """
        return self._codon_start

    def get_database_amr_gene_codon(self):
        """
        Gets the particular codon from the amr gene.
        :return: The codon.
        """
        return self._database_amr_gene_codon

    def get_database_amr_gene_amino_acid(self):
        """
        Gets the corresponding amino acid from the amr gene. If there is an insertion, returns 'ins'.
        :return: The amino acid from the amr gene.
        """
        if '-' in self.get_database_amr_gene_codon():
            return 'ins'
        else:
            return Bio.Seq.translate(self.get_database_amr_gene_codon(), table='Standard')

    def get_input_genome_amino_acid(self):
        """
        Gets the corresponding amino acid from the genome.  If there is a deletion returns 'del'.
        :return: The amino acid from the genome.
        """
        if '-' in self.get_input_genome_codon():
            return 'del'
        else:
            return Bio.Seq.translate(self.get_input_genome_codon(), table='Standard')

    def get_input_genome_codon(self):
        """
        Gets the codon from the input genome.
        :return: The codon from the input genome.
        """
        return self._input_genome_codon

    def get_mutation_position(self):
        return self.get_codon_start()

    def get_mutation_string(self):
        return self.get_database_amr_gene_codon() + ' -> ' + self.get_input_genome_codon() + ' (' + self.get_database_amr_gene_amino_acid() \
               + ' -> ' + self.get_input_genome_amino_acid() + ')'

    def get_database_amr_gene_mutation(self):
        if '-' in self.get_database_amr_gene_codon():
            return self.get_database_amr_gene_amino_acid()
        else:
            return self.get_database_amr_gene_amino_acid().upper()

    def get_input_genome_mutation(self):
        # Keep 'ins' or 'del' lowercase 
        if '-' in self.get_input_genome_codon():
            return self.get_input_genome_amino_acid()
        else:
            return self.get_input_genome_amino_acid().upper()

    def get_type(self):
        return 'codon'

    def __repr__(self):
        return (
            'CodonMutationPosition(_database_amr_gene_start={_database_amr_gene_start}, _nucleotide_position_amr_gene={_nucleotide_position_amr_gene}, '
            '_codon_start={_codon_start}, _database_amr_gene_codon={_database_amr_gene_codon}, _input_genome_codon={_input_genome_codon})').format(
            **self.__dict__)
