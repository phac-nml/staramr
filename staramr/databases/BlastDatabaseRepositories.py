import logging
import shutil
from collections import OrderedDict
from typing import Dict

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.BlastDatabaseRepository import BlastDatabaseRepository, BlastDatabaseRepositoryStripGitDir

logger = logging.getLogger('BlastDatabaseRepositories')

"""
A Class used to handle interactions with blast database repository files.
"""


class BlastDatabaseRepositories:

    def __init__(self, database_dir: str, is_dist: bool = False) -> None:
        """
        Creates a new AMRDatabaseHandler.
        :param database_dir: The root directory for the databases.
        :param is_dist: Whether or not we are building distributable versions of the blast database repositories
            (that is, should we strip out the .git directories).
        """
        self._database_dir = database_dir
        self._database_repositories = {}  # type: Dict[str,BlastDatabaseRepository]
        self._is_dist = is_dist

    def register_database_repository(self, database_name: str, git_repository_url: str) -> None:
        """
        Registers a new database repository.
        :param database_name: The name of the database.
        :param git_repository_url: The git repository url.
        :param is_dist: True if this database should be interpreted as the distributable version (no .git directory).
        :return: None
        """
        database_repository = BlastDatabaseRepository(self._database_dir, database_name,
                                                      git_repository_url)  # type: BlastDatabaseRepository
        if self._is_dist:
            database_repository = BlastDatabaseRepositoryStripGitDir(self._database_dir, database_name,
                                                                     git_repository_url)

        if database_name in self._database_repositories:
            raise Exception("A database with name [{}] already exists", database_name)
        else:
            self._database_repositories[database_name] = database_repository

    def build(self, commits: Dict[str, str] = None):
        """
        Downloads and builds new databases.
        :param commits: A map of {'database_name' : 'commit'} defining the particular commits to build.
        :return: None
        """
        for database_name in self._database_repositories:
            commit = commits.get(database_name) if commits else None
            self._database_repositories[database_name].build(commit)

    def update(self, commits: Dict[str, str] = None):
        """
        Updates an existing database to the latest revisions (or passed specific revisions).
        :param commits: A map of {'database_name' : 'commit'} defining the particular commits to update to.
        :return: None
        """
        for database_name in self._database_repositories:
            commit = commits.get(database_name) if commits else None
            self._database_repositories[database_name].update(commit)

    def remove(self):
        """
        Removes the databases stored in this directory.
        :return: None
        """
        for name, repo in self._database_repositories.items():
            repo.remove()

        shutil.rmtree(self._database_dir)

    def info(self) -> Dict[str, str]:
        """
        Gets information on the ResFinder/PointFinder databases.
        :return: Database information as a OrderedDict of key/value pairs.
        """
        info = OrderedDict()  # type: Dict[str,str]

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

    def is_at_commits(self, commits: Dict[str, str]):
        """
        Are the database repositories at the passed commits?
        :param commits: A dict of the commits {'database_name': 'commit'}.
        :return: True if the database repositories are at the passed commits (ignores repos not passed in dict). False otherwise.
        """
        for name, repo in self._database_repositories.items():
            if name in commits and not repo.is_at_commit(commits[name]):
                return False

        return True

    def is_dist(self):
        """
        Whether or not we are building distributable versions of the blast database repositories (that is, should we strip out the .git directories).
        :return: True if is_dist, False otherwise.
        """
        return self._is_dist

    @classmethod
    def create_default_repositories(cls, root_database_dir: str, is_dist: bool = False):
        """
        Class method for creating a BlastDatabaseRepositories object configured with the default repositories.
        :param database_dir: The root database directory.
        :param is_dist: Whether or not we are building distributable versions of the blast database repositories
            (that is, should we strip out the .git directories).
        :return: The BlastDatabaseRepositories.
        """
        repos = cls(root_database_dir, is_dist)
        repos.register_database_repository('resfinder', 'https://bitbucket.org/genomicepidemiology/resfinder_db.git')
        repos.register_database_repository('pointfinder',
                                           'https://bitbucket.org/genomicepidemiology/pointfinder_db.git')
        repos.register_database_repository('plasmidfinder',
                                           'https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git')

        return repos

    def build_blast_database(self, database_name: str, options: Dict[str, str] = {}) -> AbstractBlastDatabase:
        """
        Builds a staramr.blast.AbstractBlastDatabase from the given parameters.
        :param database_name: The name of the database to build.
        :param options: Options for the particular database in the form of a map {'key': 'value'}
        :return: A new staramr.blast.AbstractBlastDatabase.
        """
        if database_name not in self._database_repositories:
            raise Exception("database_name={} not registered", database_name)

        if database_name == 'resfinder':
            return ResfinderBlastDatabase(self.get_repo_dir(database_name))
        elif database_name == 'pointfinder':
            return PointfinderBlastDatabase(self.get_repo_dir(database_name), options['organism'])
        elif database_name == 'plasmidfinder':
            if options:
                return PlasmidfinderBlastDatabase(self.get_repo_dir(database_name), options['database_type'])
            else:
                return PlasmidfinderBlastDatabase(self.get_repo_dir(database_name))
        else:
            raise Exception("Unknown database name [{}]", database_name)
