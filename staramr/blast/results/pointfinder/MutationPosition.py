import abc

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class MutationPosition:

    def __init__(self, match_position, database_start, database_strand):
        """
        Creates a new MutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_start: The start coordinates of the BLAST database hit.
        :param database_strand: The strand of the BLAST database.
        :param query_strand: The strand of the BLAST query.
        """
        __metaclass__ = abc.ABCMeta

        self._check_strand(database_strand)

        self._database_start = database_start
        self._database_strand = database_strand

        if database_strand == 'plus':
            self._nucleotide_position_database = database_start + match_position
        else:
            self._nucleotide_position_database = database_start - match_position

    def _check_strand(self, strand):
        if strand not in ['plus', 'minus']:
            raise Exception("Error, strand=" + strand + " not in [plus, minus].")

    def get_nucleotide_position(self):
        """
        Gets the nucleotide position in the BLAST database (1-based coords).
        :return: The nucleotide position.
        """
        return self._nucleotide_position_database

    def get_mutation_string_short(self):
        return self.get_database_mutation() + str(self.get_mutation_position()) + self.get_query_mutation()

    @abc.abstractmethod
    def get_type(self):
        """
        Gets the type of this mutation.
        :return: The type of this mutation.
        """
        pass

    @abc.abstractmethod
    def get_mutation_position(self):
        """
        Gets the position of this mutation.
        :return: The position of this mutation.
        """
        pass

    @abc.abstractmethod
    def get_database_mutation(self):
        """
        Gets the database characters corresponding to the mutation.
        :return: The database characters..
        """
        pass

    @abc.abstractmethod
    def get_query_mutation(self):
        """
        Gets the query characters corresponding to the mutation.
        :return: The query characters..
        """
        pass

    @abc.abstractmethod
    def get_mutation_string(self):
        """
        Gets the mutation as a string.
        :return: The mutation as a string.
        """
        pass
