import logging
import shutil
from collections import OrderedDict
from typing import Dict

from staramr.databases.BlastDatabaseRepository import BlastDatabaseRepository, BlastDatabaseRepositoryStripGitDir

logger = logging.getLogger('BlastDatabaseRepositories')

"""
A Class used to handle interactions with blast database repository files.
"""


class BlastDatabaseRepositories:

    def __init__(self, database_dir: str):
        """
        Creates a new AMRDatabaseHandler.
        :param database_dir: The root directory for the databases.
        """
        self._database_dir = database_dir
        self._database_repositories = {}

    def register_database_repository(self, database_name: str, git_repository_url: str, is_dist: bool = False):
        """
        Registers a new database repository.
        :param database_name: The name of the database.
        :param git_repository_url: The git repository url.
        :param is_dist: True if this database should be interpreted as the distributable version (no .git directory).
        :return: None
        """
        if is_dist:
            database_repository = BlastDatabaseRepositoryStripGitDir(self._database_dir, database_name,
                                                                     git_repository_url)
        else:
            database_repository = BlastDatabaseRepository(self._database_dir, database_name, git_repository_url)

        if database_name in self._database_repositories:
            raise Exception("A database with name [{}] already exists", database_name)
        else:
            self._database_repositories[database_name] = database_repository

    def build(self, commits: Dict[str, str] = {}):
        """
        Downloads and builds new databases.
        :param commits: A map of {'database_name' : 'commit'} defining the particular commits to build.
        :return: None
        """
        for database_name in self._database_repositories:
            commit = commits.get(database_name)
            self._database_repositories[database_name].build(commit)

    def update(self, commits: Dict[str, str] = {}):
        """
        Updates an existing ResFinder/PointFinder database to the latest revisions (or passed specific revisions).
        :param commits: A map of {'database_name' : 'commit'} defining the particular commits to update to.
        :return: None
        """
        for database_name in self._database_repositories:
            commit = commits.get(database_name)
            self._database_repositories[database_name].update(commit)

    def remove(self):
        """
        Removes the databases stored in this directory.
        :return: None
        """
        for name, repo in self._database_repositories.items():
            repo.remove()

        shutil.rmtree(self._database_dir)

    def info(self) -> OrderedDict:
        """
        Gets information on the ResFinder/PointFinder databases.
        :return: Database information as a OrderedDict of key/value pairs.
        """
        info = OrderedDict()

        for name, repo in self._database_repositories.items():
            info.update(repo.info())

        return info

    def get_database_dir(self) -> str:
        """
        Gets the root database dir.
        :return: The root database dir.
        """
        return self._database_dir

    def get_repo_dir(self, name: str) -> str:
        """
        Gets database repo directory for the given database name.
        :param name: The database name.
        :return: The database dir for the given database name.
        """
        return self._database_repositories[name].get_git_dir()

    @classmethod
    def create_default_repositories(cls, root_database_dir: str, is_dist: bool = False):
        """
        Class method for creating a BlastDatabaseRepositories object configured with the default repositories.
        :param database_dir: The root database directory.
        :param is_dist: Whether or not we are building distributable versions of the blast database repositories
            (that is, should we strip out the .git directories).
        :return: The BlastDatabaseRepositories.
        """
        repos = cls(root_database_dir)
        repos.register_database_repository('resfinder', 'https://bitbucket.org/genomicepidemiology/resfinder_db.git',
                                           is_dist)
        repos.register_database_repository('pointfinder',
                                           'https://bitbucket.org/genomicepidemiology/pointfinder_db.git',
                                           is_dist)

        return repos
