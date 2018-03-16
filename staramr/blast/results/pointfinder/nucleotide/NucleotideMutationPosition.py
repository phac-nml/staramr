from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class NucleotideMutationPosition(CodonMutationPosition):

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
