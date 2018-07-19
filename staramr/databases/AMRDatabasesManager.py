import logging
from os import path

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.databases.AMRDatabaseHandlerStripGitDir import AMRDatabaseHandlerStripGitDir

logger = logging.getLogger('AMRDatabasesManager')

"""
A Class used to manage interactions with default and updatable ResFinder/PointFinder database installations.
"""


class AMRDatabasesManager:
    DEFAULT_RESFINDER_COMMIT = 'e8f1eb2585cd9610c4034a54ce7fc4f93aa95535'
    DEFAULT_POINTFINDER_COMMIT = '8706a6363bb29e47e0e398c53043b037c24b99a7'

    def __init__(self, database_dir, sub_dirs=False):
        """
        Builds a new AMRDatabasesManager with the passed directory.
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
        if self._sub_dirs and (force_use_git or path.exists(self._git_database_dir)):
            return AMRDatabaseHandler(self._git_database_dir)
        elif self._sub_dirs:
            return AMRDatabaseHandlerStripGitDir(self._git_strip_database_dir)
        else:
            return AMRDatabaseHandler(self._database_dir)

    def setup_default(self):
        """
        Sets up a default database.
        :return: None
        """

        if path.exists(self._git_strip_database_dir):
            logger.warning("Default database already exists in [%s]", self._git_strip_database_dir)
        else:
            logger.info("Setting up default database in [%s]", self._git_strip_database_dir)
            database_handler = AMRDatabaseHandlerStripGitDir(self._git_strip_database_dir)
            database_handler.build(resfinder_commit=self.DEFAULT_RESFINDER_COMMIT,
                                   pointfinder_commit=self.DEFAULT_POINTFINDER_COMMIT)

    def restore_default(self):
        """
        Restores the default database.
        :return: None
        """

        if path.exists(self._git_database_dir):
            logger.info("Removing database in [%s]", self._git_database_dir)
            database_handler = AMRDatabaseHandler(self._git_database_dir)
            database_handler.remove()

            if not path.exists(self._git_strip_database_dir):
                self.setup_default()

            logger.info("Restored default database to [%s]", self._git_strip_database_dir)
        else:
            if not path.exists(self._git_strip_database_dir):
                self.setup_default()
            else:
                logger.info("Default database already in use under directory [%s]", self._git_strip_database_dir)

    @classmethod
    def is_handler_default_commits(self, amr_database_handler: AMRDatabaseHandler) -> bool:
        """
        Checks whether the past database handler is linked to default commits of the database.
        :param amr_database_handler: The database handler.
        :return: True if it's setup with default commit versions, false otherwise.
        """
        database_info = amr_database_handler.info()

        return database_info['resfinder_db_commit'] == self.DEFAULT_RESFINDER_COMMIT and database_info[
            'pointfinder_db_commit'] == self.DEFAULT_POINTFINDER_COMMIT

    @classmethod
    def get_default_database_directory(cls):
        """
        Class method for getting the default database root directory.
        :return: The default database root directory.
        """
        return path.join(path.dirname(__file__), 'data')

    @classmethod
    def create_default_manager(cls):
        """
        Class method for getting the default database manager.
        :return: The default database manager.
        """
        return cls(cls.get_default_database_directory(), sub_dirs=True)
