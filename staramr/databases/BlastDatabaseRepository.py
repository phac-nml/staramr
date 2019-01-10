import logging
import shutil
import time
from collections import OrderedDict
from os import path

import git

from staramr.exceptions.DatabaseErrorException import DatabaseErrorException
from staramr.exceptions.DatabaseNotFoundException import DatabaseNotFoundException

logger = logging.getLogger('BlastDatabaseRepository')

"""
A Class used to handle interactions with the BLAST database repositories.
"""


class BlastDatabaseRepository:
    TIME_FORMAT = "%a, %d %b %Y %H:%M"

    def __init__(self, database_root_dir, database_name, git_repository_url):
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

    def build(self, commit=None):
        """
        Downloads and builds a new Blast database.
        :param commit: The specific git commit to download (if unspecified defaults to latest commit).
        :return: None
        """

        try:
            logger.info("Cloning %s db [%s] to [%s]", self._database_name, self._git_repository_url, self._git_dir)
            repo = git.repo.base.Repo.clone_from(self._git_repository_url, self._git_dir)

            if commit is not None:
                logger.info("Checking out %s commit %s", self._database_name, commit)
                repo.git.checkout(commit)
        except Exception as e:
            raise DatabaseErrorException("Could not build database in [" + self._database_dir + "]") from e

    def update(self, commit=None):
        """
        Updates an existing Blast database to the latest revisions (or passed specific revisions).
        :param commit: The specific git commit (if unspecified defaults to latest commit).
        :return: None
        """

        if not path.exists(self._git_dir):
            self.build(commit=commit)
        else:
            try:
                repo = git.Repo(self._git_dir)

                logger.info("Updating %s", self._git_dir)
                repo.heads.master.checkout()
                repo.remotes.origin.pull()

                if commit is not None:
                    logger.info("Checking out %s commit %s", self._database_name, commit)
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

    def info(self):
        """
        Gets information on the Blast databases.
        :return: Database information as a OrderedDict of key/value pairs.
        """
        data = OrderedDict()

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

    def get_database_dir(self):
        """
        Gets the root database dir.
        :return: The root database dir.
        """
        return self._database_dir

    def get_git_dir(self):
        """
        Gets the database git directory.
        :return: The database git directory.
        """
        return self._git_dir
