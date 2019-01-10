import configparser
import logging
import shutil
from collections import OrderedDict
from os import path

from staramr.databases.BlastDatabaseRepository import BlastDatabaseRepository
from staramr.exceptions.DatabaseNotFoundException import DatabaseNotFoundException

logger = logging.getLogger('BlastDatabaseRepositoryStripGitDir')

"""
A Class used to handle interactions with the BLAST database repositories, stripping out the .git directory.
"""


class BlastDatabaseRepositoryStripGitDir(BlastDatabaseRepository):
    GIT_INFO_SECTION = 'GitInfo'

    def __init__(self, database_root_dir, database_name, git_repository_url):
        """
        Creates a new BlastDatabaseRepositoryStripGitDir.
        :param database_root_dir: The root directory for both the Blast databases.
        :param database_name: A name for this database.
        :param git_repository_url: A URL to the git repository managing the database files.
        """
        super().__init__(database_root_dir, database_name, git_repository_url)

        self._git_dot_git_dir = path.join(self._git_dir, '.git')
        self._info_file = path.join(database_root_dir, database_name + '-info.ini')

    def build(self, commit=None):
        """
        Downloads and builds a new Blast database.
        :param commit: The specific git commit to download (if unspecified defaults to latest commit).
        :return: None
        """
        super().build(commit=commit)

        database_info = super().info()

        # remove directories from info as they are unimportant here
        database_info_stripped = OrderedDict(database_info)
        del database_info_stripped[self._get_info_name('dir')]

        self._write_database_info_to_file(database_info_stripped, self._info_file)

        logger.info("Removing %s", self._git_dot_git_dir)
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

    def update(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Updates an existing ResFinder/PointFinder database to the latest revisions (or passed specific revisions).
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """
        raise Exception("Cannot update when .git directory has been removed")

    def info(self):
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
