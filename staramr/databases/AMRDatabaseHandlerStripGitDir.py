import configparser
import logging
import shutil
from collections import OrderedDict
from os import path

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.exceptions.DatabaseNotFoundException import DatabaseNotFoundException

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
        database_info_stripped = OrderedDict(database_info)
        del database_info_stripped['resfinder_db_dir']
        del database_info_stripped['pointfinder_db_dir']

        self._write_database_info_to_file(database_info_stripped, self._info_file)

        logger.info("Removing %s", self._resfinder_dir_git)
        shutil.rmtree(self._resfinder_dir_git)
        logger.info("Removing %s", self._pointfinder_dir_git)
        shutil.rmtree(self._pointfinder_dir_git)

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
            data['resfinder_db_dir'] = self._resfinder_dir
            data['pointfinder_db_dir'] = self._pointfinder_dir

            # re-order all fields
            data.move_to_end('resfinder_db_dir', last=True)
            data.move_to_end('resfinder_db_url', last=True)
            data.move_to_end('resfinder_db_commit', last=True)
            data.move_to_end('resfinder_db_date', last=True)

            data.move_to_end('pointfinder_db_dir', last=True)
            data.move_to_end('pointfinder_db_url', last=True)
            data.move_to_end('pointfinder_db_commit', last=True)
            data.move_to_end('pointfinder_db_date', last=True)

            return data
        except FileNotFoundError as e:
            raise DatabaseNotFoundException('Database could not be found in [' + self._database_dir + ']') from e
