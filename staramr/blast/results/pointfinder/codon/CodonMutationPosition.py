import math

import Bio.Seq

from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a codon-based mutation for PointFinder.
"""


class CodonMutationPosition(MutationPosition):

    def __init__(self, match_position, database_string, query_string, database_start, database_frame, query_frame):
        """
        Creates a new CodonMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_string: The BLAST database string.
        :param query_string: The BLAST query string.
        :param database_start: The start coordinates of the BLAST database hit.
        :param database_frame: The frame (strand) of the BLAST database.
        :param query_frame: The frame (strand) of the BLAST query.
        """
        super().__init__(match_position, database_start, database_frame, query_frame)

        self._codon_start_database = math.ceil(self._nucleotide_position_database / 3)
        frame_shift = (self._nucleotide_position_database - 1) % 3

        self._database_codon = self._find_codon(database_string, match_position, database_frame, frame_shift)
        self._query_codon = self._find_codon(query_string, match_position, database_frame, frame_shift)

    def _find_codon(self, string, match_position, frame, frame_shift):
        if frame == 1:
            codon_start_index = match_position - frame_shift
            return string[codon_start_index:(codon_start_index + 3)].upper()
        else:
            codon_end_index = match_position + frame_shift
            return Bio.Seq.reverse_complement(string[(codon_end_index - 3 + 1):(codon_end_index + 1)].upper())

    def get_codon_start(self):
        """
        Gets the codon start in the BLAST database for PointFinder (1-based).
        :return: The codon start.
        """
        return self._codon_start_database

    def get_database_codon(self):
        """
        Gets the particular codon from the BLAST database.
        :return: The codon.
        """
        return self._database_codon

    def get_database_amino_acid(self):
        """
        Gets the corresponding amino acid from the BLAST database. If there is an indel, returns 'X'.
        :return: The amino acid from the BLAST database.
        """
        if '-' in self.get_database_codon():
            return 'X'
        else:
            return Bio.Seq.translate(self.get_database_codon(), table='Standard')

    def get_query_amino_acid(self):
        """
        Gets the corresponding amino acid from the BLAST query.  If there is an indel returns 'X'.
        :return: The amino acid from the BLAST query.
        """
        if '-' in self.get_query_codon():
            return 'X'
        else:
            return Bio.Seq.translate(self.get_query_codon(), table='Standard')

    def get_query_codon(self):
        """
        Gets the codon from the BLAST query.
        :return: The codon from the BLAST query.
        """
        return self._query_codon

    def get_mutation_position(self):
        return self.get_codon_start()

    def get_mutation_string(self):
        return self.get_database_codon() + ' -> ' + self.get_query_codon() + ' (' + self.get_database_amino_acid() \
               + ' -> ' + self.get_query_amino_acid() + ')'

    def get_database_mutation(self):
        return self.get_database_amino_acid().upper()

    def get_query_mutation(self):
        return self.get_query_amino_acid().upper()

    def get_type(self):
        return 'codon'

    def __repr__(self):
        return "[database_start=" + str(self._database_start) + ", database_frame=" + str(
            self._database_frame) + ", query_frame=" + str(self._query_frame) + ", nucleotide_position=" \
               + str(self._nucleotide_position_database) + ", codon_start=" + str(self._codon_start_database) \
               + ", codon=" + self._database_codon + "]"
