import logging
from os import path

from staramr.databases.BlastDatabaseRepositories import BlastDatabaseRepositories

logger = logging.getLogger('AMRDatabasesManager')

"""
A Class used to manage interactions with default and updatable ResFinder/PointFinder/PlasmidFinder database installations.
"""


class AMRDatabasesManager:
    # Update to commits corresponding to dates listed on <https://cge.cbs.dtu.dk/services/ResFinder/> (and PlasmidFinder)
    # As of May 26, 2022
    DEFAULT_COMMITS = {
        'resfinder': 'd1e607b8989260c7b6a3fbce8fa3204ecfc09022', # 2024-08-06
        'pointfinder': '694919f59a38980204009e7ade76bf319cb7df0b', # 2024-08-08
        'plasmidfinder': 'c18e08c17a5988d4f075fc1171636e47546a323d', # 2023-01-18
    }

    def __init__(self, database_dir: str, sub_dirs: bool = False) -> None:
        """
        Builds a new AMRDatabasesManager with the passed directory.
        :param database_dir: The directory containing the ResFinder/PointFinder/PlasmidFinder databases.
        :param sub_dirs: If True, assumes we are using subdirectories to store databases
                            and searching for stripped git directories.
        """
        self._database_dir = database_dir
        self._git_database_dir = path.join(database_dir, 'update')
        self._git_strip_database_dir = path.join(database_dir, 'dist')
        self._sub_dirs = sub_dirs

    def get_database_repos(self, force_use_git: bool = False) -> BlastDatabaseRepositories:
        """
        Gets the appropriate database repositories.
        :param force_use_git: Force use of git database repos.
        :return: The database repos object.
        """
        if self._sub_dirs and (force_use_git or path.exists(self._git_database_dir)):
            return BlastDatabaseRepositories.create_default_repositories(self._git_database_dir)
        elif self._sub_dirs:
            return BlastDatabaseRepositories.create_default_repositories(self._git_strip_database_dir, is_dist=True)
        else:
            return BlastDatabaseRepositories.create_default_repositories(self._git_database_dir)

    def setup_default(self):
        """
        Sets up a default database.
        :return: None
        """

        if path.exists(self._git_strip_database_dir):
            logger.warning("Default database already exists in [%s]", self._git_strip_database_dir)
        else:
            logger.info("Setting up default database in [%s]", self._git_strip_database_dir)
            database_repos = BlastDatabaseRepositories.create_default_repositories(self._git_strip_database_dir,
                                                                                   is_dist=True)
            database_repos.build(self.DEFAULT_COMMITS)

    def restore_default(self):
        """
        Restores the default database.
        :return: None
        """

        if path.exists(self._git_database_dir):
            logger.info("Removing database in [%s]", self._git_database_dir)
            database_repos = BlastDatabaseRepositories.create_default_repositories(self._git_database_dir)
            database_repos.remove()

            if not path.exists(self._git_strip_database_dir):
                self.setup_default()

            logger.info("Restored default database to [%s]", self._git_strip_database_dir)
        else:
            if not path.exists(self._git_strip_database_dir):
                self.setup_default()
            else:
                logger.info("Default database already in use under directory [%s]", self._git_strip_database_dir)

    @classmethod
    def is_database_repos_default_commits(self, database_repos: BlastDatabaseRepositories) -> bool:
        """
        Checks whether the past database repos handler is linked to default commits of the database.
        :param database_repos: The database repos handler.
        :return: True if it's setup with default commit versions, false otherwise.
        """
        return database_repos.is_at_commits(self.DEFAULT_COMMITS)

    @classmethod
    def get_default_database_directory(cls) -> str:
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
