import abc

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class MutationPosition:

    def __init__(self, match_position, database_amr_gene_start):
        """
        Creates a new MutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_amr_gene_start: The start coordinates of the amr gene from the BLAST hit.
        """
        __metaclass__ = abc.ABCMeta

        self._database_amr_gene_start = database_amr_gene_start
        self._nucleotide_position_amr_gene = database_amr_gene_start + match_position

    def get_nucleotide_position(self):
        """
        Gets the nucleotide position in the amr gene (1-based coords).
        :return: The nucleotide position.
        """
        return self._nucleotide_position_amr_gene

    def get_mutation_string_short(self):
        return self.get_database_amr_gene_mutation() + str(
            self.get_mutation_position()) + self.get_input_genome_mutation()

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
    def get_database_amr_gene_mutation(self):
        """
        Gets the database amr gene characters corresponding to the mutation.
        :return: The database amr gene characters.
        """
        pass

    @abc.abstractmethod
    def get_input_genome_mutation(self):
        """
        Gets the input genome characters corresponding to the mutation.
        :return: The input genome characters.
        """
        pass

    @abc.abstractmethod
    def get_mutation_string(self):
        """
        Gets the mutation as a string.
        :return: The mutation as a string.
        """
        pass
