import Bio.Seq

from staramr.blast.results.pointfinder.MutationPosition import MutationPosition

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class NucleotideMutationPosition(MutationPosition):

    def __init__(self, match_position, database_string, query_string, database_start, database_frame, query_frame):
        """
        Creates a new NucleotideMutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_string: The BLAST database string.
        :param query_string: The BLAST query string.
        :param database_start: The start coordinates of the BLAST database hit.
        :param database_frame: The frame (strand) of the BLAST database.
        :param query_frame: The frame (strand) of the BLAST query.
        """
        super().__init__(match_position, database_start, database_frame, query_frame)

        self._database_nucleotide = database_string[match_position].upper()
        self._query_nucleotide = query_string[match_position].upper()

        if self._database_frame == -1:
            self._database_nucleotide = Bio.Seq.reverse_complement(self._database_nucleotide)
            self._query_nucleotide = Bio.Seq.reverse_complement(self._query_nucleotide)

    def get_type(self):
        return 'nucleotide'

    def get_mutation_position(self):
        return self.get_nucleotide_position()

    def get_database_mutation(self):
        return self._database_nucleotide

    def get_query_mutation(self):
        return self._query_nucleotide

    def get_mutation_string(self):
        return self.get_database_mutation() + ' -> ' + self.get_query_mutation()

    def __repr__(self):
        return "[database_start=" + str(self._database_start) + ", database_frame=" + str(
            self._database_frame) + ", query_frame=" + str(self._query_frame) + ", nucleotide_position=" \
               + str(self._nucleotide_position_database) + ", mutation_start=" + str(self.get_mutation_position()) \
               + ", mutation=" + self.get_mutation_string() + "]"
