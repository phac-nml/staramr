import logging
from os import path

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.databases.AMRDatabaseHandlerStripGitDir import AMRDatabaseHandlerStripGitDir

logger = logging.getLogger('AMRDatabaseHandlerFactory')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files.
"""


class AMRDatabaseHandlerFactory:

    def __init__(self, database_dir):
        """
        Builds a new AMRDatabaseHandlerFactory with the passed directory.
        :param database_dir: The directory containing the ResFinder/PointFinder databases.
        """
        self._database_dir = database_dir
        self._git_database_dir = path.join(database_dir, 'git')
        self._git_strip_database_dir = path.join(database_dir, 'strip')

    def get_database_handler(self, strip_git=False):
        """
        Gets the appropriate database handler.
        :param strip_git: Whether or not to strip out the .git directory from the databases.
        :return: The database handler.
        """
        if strip_git:
            return AMRDatabaseHandlerStripGitDir(self._git_strip_database_dir)
        else:
            return AMRDatabaseHandler(self._git_database_dir)

    @classmethod
    def get_default_database_directory(cls):
        """
        Class method for getting the default database root directory.
        :return: The default database root directory.
        """
        return path.join(path.dirname(__file__), 'data')

    @classmethod
    def create_default_factory(cls):
        """
        Class method for getting the default database handler factory.
        :return: The default database handler factory.
        """
        return cls(cls.get_default_database_directory())
