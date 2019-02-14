import configparser
import logging
import shutil
import time
from collections import OrderedDict
from os import path
from typing import Dict

import git

from staramr.exceptions.DatabaseErrorException import DatabaseErrorException
from staramr.exceptions.DatabaseNotFoundException import DatabaseNotFoundException

"""
A Class used to handle interactions with the BLAST database repositories.
"""


class BlastDatabaseRepository:
    TIME_FORMAT = "%a, %d %b %Y %H:%M"
    LOGGER = logging.getLogger('BlastDatabaseRepository')

    def __init__(self, database_root_dir: str, database_name: str, git_repository_url: str) -> None:
        """
        Creates a new BlastDatabaseRepository.
        :param database_root_dir: The root directory for both the Blast databases.
        :param database_name: A name for this database.
        :param git_repository_url: A URL to the git repository managing the database files.
        """
        self._database_dir = database_root_dir
        self._database_name = database_name
        self._git_repository_url = git_repository_url

        self._git_dir = path.join(database_root_dir, database_name)

    def build(self, commit: str = None) -> None:
        """
        Downloads and builds a new Blast database.
        :param commit: The specific git commit to download. Defaults to latest commit.
        :return: None
        """

        try:
            self.LOGGER.info("Cloning %s db [%s] to [%s]", self._database_name, self._git_repository_url, self._git_dir)
            repo = git.repo.base.Repo.clone_from(self._git_repository_url, self._git_dir)

            if commit is not None:
                self.LOGGER.info("Checking out %s commit %s", self._database_name, commit)
                repo.git.checkout(commit)
        except Exception as e:
            raise DatabaseErrorException("Could not build database in [" + self._database_dir + "]") from e

    def update(self, commit: str = None) -> None:
        """
        Updates an existing Blast database to the latest revisions (or passed specific revisions).
        :param commit: The specific git commit to update to. Defaults to latest commit.
        :return: None
        """

        if not path.exists(self._git_dir):
            self.build(commit=commit)
        else:
            try:
                repo = git.Repo(self._git_dir)

                self.LOGGER.info("Updating %s", self._git_dir)
                repo.heads.master.checkout()
                repo.remotes.origin.pull()

                if commit is not None:
                    self.LOGGER.info("Checking out %s commit %s", self._database_name, commit)
                    repo.git.checkout(commit)

                    repo.git.reset('--hard')
            except Exception as e:
                raise DatabaseErrorException("Could not build database in [" + self._database_dir + "]") from e

    def remove(self):
        """
        Removes the databases stored in this directory.
        :return: None
        """
        shutil.rmtree(self._git_dir)

    def is_at_commit(self, commit: str = None) -> bool:
        """
        Determines whether this database repo is at the specified commit.
        :param commit: The commit to check.
        :return: True if the database is at the specified commit, otherwise False.
        """
        return self.info()[self._get_info_name('commit')] == commit

    def info(self) -> Dict[str, str]:
        """
        Gets information on the Blast databases.
        :return: Database information as a OrderedDict of key/value pairs.
        """
        data = OrderedDict()  # type: Dict[str,str]

        try:
            repo = git.Repo(self._git_dir)
            repo_head = repo.commit('HEAD')

            data[self._get_info_name('dir')] = self._git_dir
            data[self._get_info_name('url')] = self._git_repository_url
            data[self._get_info_name('commit')] = str(repo_head)
            data[self._get_info_name('date')] = time.strftime(self.TIME_FORMAT,
                                                              time.gmtime(repo_head.committed_date))

        except git.exc.NoSuchPathError as e:
            raise DatabaseNotFoundException('Invalid database in [' + self._database_dir + ']') from e

        return data

    def _get_info_name(self, info_type):
        return self._database_name + '_db_' + info_type

    def get_database_dir(self) -> str:
        """
        Gets the root database dir.
        :return: The root database dir.
        """
        return self._database_dir

    def get_git_dir(self) -> str:
        """
        Gets the database git directory.
        :return: The database git directory.
        """
        return self._git_dir


"""
A Class used to handle interactions with the BLAST database repositories, stripping out the .git directory.
"""


class BlastDatabaseRepositoryStripGitDir(BlastDatabaseRepository):
    GIT_INFO_SECTION = 'GitInfo'
    LOGGER = logging.getLogger('BlastDatabaseRepositoryStripGitDir')

    def __init__(self, database_root_dir: str, database_name: str, git_repository_url: str) -> None:
        """
        Creates a new BlastDatabaseRepositoryStripGitDir.
        :param database_root_dir: The root directory for both the Blast databases.
        :param database_name: A name for this database.
        :param git_repository_url: A URL to the git repository managing the database files.
        """
        super().__init__(database_root_dir, database_name, git_repository_url)

        self._git_dot_git_dir = path.join(self._git_dir, '.git')
        self._info_file = path.join(database_root_dir, database_name + '-info.ini')

    def build(self, commit: str = None):
        """
        Downloads and builds a new Blast database.
        :param commit: The specific git commit to download. Defaults to latest commit.
        :return: None
        """
        super().build(commit=commit)

        database_info = super().info()

        # remove directories from info as they are unimportant here
        database_info_stripped = OrderedDict(database_info)
        del database_info_stripped[self._get_info_name('dir')]

        self._write_database_info_to_file(database_info_stripped, self._info_file)

        self.LOGGER.info("Removing %s", self._git_dot_git_dir)
        shutil.rmtree(self._git_dot_git_dir)

    def _write_database_info_to_file(self, database_info, file):
        config = configparser.ConfigParser()
        config[self.GIT_INFO_SECTION] = database_info

        with open(file, 'w') as file_handle:
            config.write(file_handle)

    def _read_database_info_from_file(self, file):
        config = configparser.ConfigParser()
        config.read(file)
        return OrderedDict(config[self.GIT_INFO_SECTION])

    def update(self, commit: str = None) -> None:
        """
        Updates an existing Blast database to the latest revisions (or passed specific revisions).
        :param commit: The commit to update to.
        :return: None
        """
        raise Exception("Cannot update when .git directory has been removed")

    def info(self) -> Dict[str, str]:
        """
        Gets information on the ResFinder/PointFinder databases.
        :return: Database information as a list containing key/value pairs.
        """

        try:
            data = self._read_database_info_from_file(self._info_file)
            data[self._get_info_name('dir')] = self._git_dir

            # re-order all fields
            data.move_to_end(self._get_info_name('dir'), last=True)
            data.move_to_end(self._get_info_name('url'), last=True)
            data.move_to_end(self._get_info_name('commit'), last=True)
            data.move_to_end(self._get_info_name('date'), last=True)

            return data
        except FileNotFoundError as e:
            raise DatabaseNotFoundException('Database could not be found in [' + self.get_database_dir() + ']') from e
