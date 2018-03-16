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
        super().__init__(match_position, database_string, query_string, database_start, database_frame, query_frame)

        self._database_nucleotide = database_string[match_position].upper()
        self._query_nucleotide = query_string[match_position].upper()

    def get_type(self):
        return 'nucleotide'

    def get_mutation_position(self):
        return self.get_nucleotide_position()

    def get_database_mutation(self):
        self._database_nucleotide

    def get_query_mutation(self):
        self._query_nucleotide

    def get_mutation_string(self):
        return self.get_database_mutation() + ' -> ' + self.get_query_mutation()