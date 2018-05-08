import configparser
import logging
import shutil
from os import path

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.exceptions.DatabaseNotFoundException import DatabaseNotFoundException
from staramr.exceptions.DatabaseErrorException import DatabaseErrorException

logger = logging.getLogger('AMRDatabaseHandlerStripGitDir')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files, stripping out the .git directory.
"""


class AMRDatabaseHandlerStripGitDir(AMRDatabaseHandler):
    GIT_INFO_SECTION = 'GitInfo'

    def __init__(self, database_dir):
        """
        Creates a new AMRDatabaseHandlerStripGitDir.
        :param database_dir: The root directory for both the ResFinder/PointFinder databases.
        """
        super().__init__(database_dir)

        self._resfinder_dir_git = path.join(self._resfinder_dir, '.git')
        self._pointfinder_dir_git = path.join(self._pointfinder_dir, '.git')
        self._info_file = path.join(database_dir, 'info.ini')

    def build(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Downloads and builds a new ResFinder/PointFinder database.
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """
        super().build(resfinder_commit=resfinder_commit, pointfinder_commit=pointfinder_commit)

        database_info = super().info()

        # remove directories from info as they are unimportant here
        database_info_stripped = []
        for i in database_info:
            if i[0] != 'resfinder_db_dir' and i[0] != 'pointfinder_db_dir':
                database_info_stripped.append(i)

        self._write_database_info_to_file(database_info_stripped, self._info_file)

        logger.info("Removing " + self._resfinder_dir_git)
        shutil.rmtree(self._resfinder_dir_git)
        logger.info("Removing " + self._pointfinder_dir_git)
        shutil.rmtree(self._pointfinder_dir_git)

    def _write_database_info_to_file(self, database_info, file):
        config = configparser.ConfigParser()
        config[self.GIT_INFO_SECTION] = {k: v for k, v in database_info}

        with open(file, 'w') as file_handle:
            config.write(file_handle)

    def _read_database_info_from_file(self, file):
        config = configparser.ConfigParser()
        config.read(file)
        git_info = config[self.GIT_INFO_SECTION]
        return [[k, git_info[k]] for k in git_info]

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

        if self.is_error():
            raise DatabaseErrorException('Database [' + self._database_dir + '] is in an error state')

        try:
            data = self._read_database_info_from_file(self._info_file)
            data.insert(0, ['resfinder_db_dir', self._resfinder_dir])
            data.insert(4, ['pointfinder_db_dir', self._pointfinder_dir])

            return data
        except FileNotFoundError as e:
            raise DatabaseNotFoundException('Database could not be found in [' + self._database_dir + ']') from e
