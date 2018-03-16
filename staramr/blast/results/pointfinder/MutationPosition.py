import abc

"""
A Class defining a nucleotide-based mutation for PointFinder.
"""


class MutationPosition:

    def __init__(self, match_position, database_start, database_frame, query_frame):
        """
        Creates a new MutationPosition.
        :param match_position: The particular position (0-based index) of the BLAST match string for this mutation.
        :param database_start: The start coordinates of the BLAST database hit.
        :param database_frame: The frame (strand) of the BLAST database.
        :param query_frame: The frame (strand) of the BLAST query.
        """
        __metaclass__ = abc.ABCMeta

        # TODO: I realise that frame and strand are different, and that I want to account for the strand here since I'm
        #  parsing BLASTN results <http://biopython.org/DIST/docs/api/Bio.Blast.Record.HSP-class.html>.
        #  But from all my testing with BioPython/BLAST parsing, it looks like even though the documentation refers
        #  to a "strand", the values I want to use as "strand" are accessible in BioPython under the "frame" labels.
        #  I need to look more into this.
        self._check_frames(database_frame, query_frame)

        self._database_start = database_start
        self._database_frame = database_frame
        self._query_frame = query_frame

        if database_frame == 1:
            self._nucleotide_position_database = database_start + match_position
        else:
            self._nucleotide_position_database = database_start - match_position

    def _check_frames(self, database_frame, query_frame):
        self._check_frame(database_frame)
        self._check_frame(query_frame)

    def _check_frame(self, frame):
        if frame not in [1, -1]:
            raise Exception("Error, frame=" + frame + " not in [1, -1].")

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
