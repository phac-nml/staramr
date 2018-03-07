import abc

"""
An Abstract Class for interacting with particular ResFinder/PointFinder databases (fasta files, etc).
"""


class AbstractBlastDatabase:
    fasta_suffix = ".fsa"

    def __init__(self, database_dir):
        """
        Creates a new AbstractBlastDatabase
        :param database_dir: The directory containing the database files.
        """
        __metaclass__ = abc.ABCMeta
        self.database_dir = database_dir

    @abc.abstractmethod
    def get_database_names(self):
        """
        Get the names of the databases (fasta files) used for BLAST.
        :return: The names of the databases.
        """
        pass

    @abc.abstractmethod
    def get_path(self, database_name):
        """
        Gets the path to a particular database with the given name.
        :param database_name: The name of the database.
        :return: The path to the database (fasta) file.
        """
        pass

    def get_database_paths(self):
        """
        Gets a list of all database paths.
        :return: A list of all database (fasta file) paths.
        """
        return [self.get_path(x) for x in self.get_database_names()]
