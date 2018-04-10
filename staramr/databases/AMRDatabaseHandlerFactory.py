import logging
from os import path

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.databases.AMRDatabaseHandlerStripGitDir import AMRDatabaseHandlerStripGitDir

logger = logging.getLogger('AMRDatabaseHandlerFactory')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files.
"""


class AMRDatabaseHandlerFactory:

    def __init__(self, database_dir, sub_dirs=False):
        """
        Builds a new AMRDatabaseHandlerFactory with the passed directory.
        :param database_dir: The directory containing the ResFinder/PointFinder databases.
        :param sub_dirs: If True, assumes we are using subdirectories to store databases
                            and searching for stripped git directories.
        """
        self._database_dir = database_dir
        self._git_database_dir = path.join(database_dir, 'update')
        self._git_strip_database_dir = path.join(database_dir, 'dist')
        self._sub_dirs = sub_dirs

    def get_database_handler(self, force_use_git=False):
        """
        Gets the appropriate database handler.
        :param force_use_git: Force use of git database handler.
        :return: The database handler.
        """
        if self._sub_dirs:
            if force_use_git or path.exists(self._git_database_dir):
                return AMRDatabaseHandler(self._git_database_dir)
            else:
                return AMRDatabaseHandlerStripGitDir(self._git_strip_database_dir)
        else:
            return AMRDatabaseHandler(self._database_dir)

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
        return cls(cls.get_default_database_directory(), sub_dirs=True)
