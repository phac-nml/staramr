import logging
import shutil
from os import path

import pandas

import staramr.Utils as Utils
from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler

logger = logging.getLogger('AMRDatabaseHandlerStripGitDir')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files, stripping out the .git directory.
"""


class AMRDatabaseHandlerStripGitDir(AMRDatabaseHandler):

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
        file_handle = open(file, 'w')
        file_handle.write(Utils.get_string_with_spacing(database_info))
        file_handle.close()

    def _read_database_info_from_file(self, file):
        return pandas.read_csv(file, sep="=", index_col=False, header=None)

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
        data = self._read_database_info_from_file(self._info_file)
        data_matrix = data.as_matrix().tolist()
        data_matrix.insert(0, ['resfinder_db_dir', self._resfinder_dir])
        data_matrix.insert(3, ['pointfinder_db_dir', self._pointfinder_dir])

        return data_matrix
